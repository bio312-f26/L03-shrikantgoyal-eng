#!/usr/bin/env python3
"""Validate selected sequences and prepare the Lab 4 FASTA and summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import tempfile
from collections import OrderedDict
from pathlib import Path


BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qcovs"]
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWYU")


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        temporary = handle.name
        handle.write(content)
    os.replace(temporary, path)
    path.chmod(0o644)


def tsv_text(fieldnames: list[str], rows: list[dict[str, str]]) -> str:
    from io import StringIO

    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def load_species(path: Path) -> OrderedDict[str, dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty species key: {path}")
    species: OrderedDict[str, dict[str, str]] = OrderedDict()
    taxon_ids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            code = row.get("species_short_code", "")
            taxon = row.get("NCBI_taxon_id", "")
            if not code or code in species:
                raise ValueError(f"duplicate or empty species code: {code!r}")
            if not taxon or taxon in taxon_ids:
                raise ValueError(f"duplicate or empty taxon ID: {taxon!r}")
            species[code] = row
            taxon_ids.add(taxon)
    if not species:
        raise ValueError(f"no species rows in {path}")
    return species


def read_candidates(path: Path, evalue_cutoff: float, identity: float, coverage: float) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty filtered BLAST table: {path}")
    rows: list[dict[str, str]] = []
    query_ids: set[str] = set()
    seen_subjects: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            values = raw.rstrip("\n\r").split("\t")
            if len(values) != len(BLAST_FIELDS):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(BLAST_FIELDS)} tab-separated fields, found {len(values)}"
                )
            row = dict(zip(BLAST_FIELDS, values))
            if not row["qseqid"] or not row["sseqid"]:
                raise ValueError(f"{path}:{line_number}: empty query or subject ID")
            if row["sseqid"] in seen_subjects:
                raise ValueError(f"{path}:{line_number}: duplicate subject ID {row['sseqid']!r}")
            query_ids.add(row["qseqid"])
            seen_subjects.add(row["sseqid"])
            try:
                numeric = {field: float(row[field]) for field in BLAST_FIELDS[2:]}
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: nonnumeric BLAST value") from error
            if any(not math.isfinite(value) for value in numeric.values()):
                raise ValueError(f"{path}:{line_number}: non-finite BLAST value")
            if (
                numeric["evalue"] > evalue_cutoff
                or numeric["pident"] < identity
                or numeric["qcovs"] < coverage
            ):
                raise ValueError(
                    f"{path}:{line_number}: row does not meet E-value <= {evalue_cutoff:g}, "
                    f"identity >= {identity:g}, and query coverage >= {coverage:g}"
                )
            rows.append(row)
    if len(query_ids) != 1:
        raise ValueError(f"filtered table contains {len(query_ids)} query IDs; expected one")
    return rows


def read_candidate_ids(path: Path) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty candidate-ID list: {path}")
    identifiers = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"duplicate ID in {path}")
    return identifiers


def read_selected_fasta(path: Path) -> OrderedDict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty selected-sequence FASTA: {path}")
    records: OrderedDict[str, str] = OrderedDict()
    header: str | None = None
    sequence: list[str] = []

    def store() -> None:
        if header is None:
            return
        if header in records:
            raise ValueError(f"duplicate ID in selected-sequence FASTA: {header}")
        joined = "".join(sequence).upper()
        if not joined:
            raise ValueError(f"empty sequence for {header}")
        unexpected = sorted(set(joined) - STANDARD_AA)
        if unexpected:
            raise ValueError(f"unexpected amino-acid characters for {header}: {''.join(unexpected)}")
        records[header] = joined

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:].split()[0]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
                sequence.append(line)
        store()
    return records


def load_lookup(path: Path, selected_keys: set[tuple[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty protein lookup: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    found: dict[tuple[str, str], dict[str, str]] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "species_short_code",
            "uniprot_accession",
            "gene_symbol",
            "protein_name",
            "sequence_length",
            "sequence_sha256",
            "reviewed_status",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"protein lookup is missing columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            key = (row["species_short_code"], row["uniprot_accession"])
            if key not in selected_keys:
                continue
            if key in found:
                raise ValueError(f"duplicate species/accession in lookup: {key[0]} {key[1]}")
            found[key] = row
    missing = selected_keys - found.keys()
    if missing:
        detail = ", ".join(f"{species}|{accession}" for species, accession in sorted(missing))
        raise ValueError(f"selected records absent from protein lookup: {detail}")
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a Lab 4 FASTA from sequences retrieved with seqkit.")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--candidate-ids", required=True, type=Path)
    parser.add_argument("--source-fasta", required=True, type=Path)
    parser.add_argument("--species-key", required=True, type=Path)
    parser.add_argument("--protein-lookup", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", required=True)
    parser.add_argument("--evalue", type=float, default=1e-10)
    parser.add_argument("--min-identity", type=float, default=35.0)
    parser.add_argument("--min-query-coverage", type=float, default=70.0)
    args = parser.parse_args()

    project_id = args.project_id.strip().upper()
    if not project_id or "|" in project_id or any(character.isspace() for character in project_id):
        raise SystemExit("ERROR: project ID must be a nonempty, pipe-free word")
    if not args.basename or "/" in args.basename or "|" in args.basename:
        raise SystemExit("ERROR: basename must be a simple pipe-free filename stem")

    try:
        species = load_species(args.species_key)
        candidates = read_candidates(
            args.candidate_table, args.evalue, args.min_identity, args.min_query_coverage
        )
        candidate_ids = read_candidate_ids(args.candidate_ids)
        table_ids = [row["sseqid"] for row in candidates]
        if set(candidate_ids) != set(table_ids):
            raise ValueError("candidate-ID list does not match subject IDs in the filtered BLAST table")

        selected_sequences = read_selected_fasta(args.source_fasta)
        if set(selected_sequences) != set(candidate_ids):
            missing = set(candidate_ids) - set(selected_sequences)
            extra = set(selected_sequences) - set(candidate_ids)
            raise ValueError(
                f"selected-sequence FASTA/ID-list mismatch: {len(missing)} missing and {len(extra)} extra"
            )

        parsed: dict[str, tuple[str, str, str]] = {}
        selected_keys: set[tuple[str, str]] = set()
        for source_id in selected_sequences:
            fields = source_id.split("|")
            if len(fields) != 3:
                raise ValueError(
                    f"wrong source-ID schema for {source_id!r}; expected SpeciesCode|UniProtAccession|GeneSymbol"
                )
            species_code, accession, gene_symbol = fields
            if species_code not in species:
                raise ValueError(f"unknown species code in selected ID {source_id!r}: {species_code}")
            key = (species_code, accession)
            if key in selected_keys:
                raise ValueError(f"duplicate selected accession within species: {species_code} {accession}")
            selected_keys.add(key)
            parsed[source_id] = (species_code, accession, gene_symbol)

        lookup = load_lookup(args.protein_lookup, selected_keys)
        counts = {code: 0 for code in species}
        mapping_rows: list[dict[str, str]] = []
        fasta_parts: list[str] = []

        for source_id, sequence in selected_sequences.items():
            species_code, accession, source_gene_symbol = parsed[source_id]
            row = lookup[(species_code, accession)]
            if source_gene_symbol != (row["gene_symbol"] or "NA"):
                raise ValueError(f"gene-symbol mismatch for {source_id}")
            observed_hash = hashlib.sha256(sequence.encode()).hexdigest()
            if int(row["sequence_length"]) != len(sequence) or row["sequence_sha256"] != observed_hash:
                raise ValueError(f"sequence length/hash mismatch for {source_id}")
            taxon_id = species[species_code]["NCBI_taxon_id"]
            project_fasta_id = f"{project_id}|{accession}|{taxon_id}|{species_code}"
            fasta_parts.append(f">{project_fasta_id}\n")
            fasta_parts.extend(sequence[start : start + 80] + "\n" for start in range(0, len(sequence), 80))
            counts[species_code] += 1
            mapping_rows.append(
                {
                    "project_id": project_id,
                    "project_fasta_id": project_fasta_id,
                    "source_id": source_id,
                    "species_short_code": species_code,
                    "uniprot_accession": accession,
                    "NCBI_taxon_id": taxon_id,
                    "original_gene_symbol": source_gene_symbol,
                    "protein_name": row["protein_name"],
                    "reviewed_status": row["reviewed_status"],
                    "sequence_length": str(len(sequence)),
                    "sequence_sha256": observed_hash,
                }
            )

        count_rows = [
            {
                "species_short_code": code,
                "common_name": species[code]["common_name"],
                "candidate_protein_records": str(counts[code]),
            }
            for code in species
        ]
        mapping_fields = [
            "project_id",
            "project_fasta_id",
            "source_id",
            "species_short_code",
            "uniprot_accession",
            "NCBI_taxon_id",
            "original_gene_symbol",
            "protein_name",
            "reviewed_status",
            "sequence_length",
            "sequence_sha256",
        ]
        args.output_dir.mkdir(parents=True, exist_ok=True)
        atomic_text(
            args.output_dir / f"{args.basename}.species_counts.tsv",
            tsv_text(["species_short_code", "common_name", "candidate_protein_records"], count_rows),
        )
        atomic_text(
            args.output_dir / f"{args.basename}.identifier_mapping.tsv",
            tsv_text(mapping_fields, mapping_rows),
        )
        atomic_text(args.output_dir / f"{args.basename}.homologs.fas", "".join(fasta_parts))
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error

    print(f"PASS: {len(selected_sequences)} selected proteins prepared for Lab 4")
    print(f"Counts: {args.output_dir / f'{args.basename}.species_counts.tsv'}")
    print(f"Mapping: {args.output_dir / f'{args.basename}.identifier_mapping.tsv'}")
    print(f"FASTA: {args.output_dir / f'{args.basename}.homologs.fas'}")


if __name__ == "__main__":
    main()
