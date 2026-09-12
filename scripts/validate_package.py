#!/usr/bin/env python3
"""Validate the frozen Lab 3 package without modifying source data."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWYU")
ACCESSION_RE = re.compile(r"^[A-Z0-9]+(?:-[0-9]+)?$")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty checksum manifest: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "packaged_file",
        "packaged_sha256",
        "decompressed_sha256",
        "decompressed_bytes",
        "record_count",
        "record_type",
    }
    fields = set(rows[0]) if rows else set()
    if not rows or not required.issubset(fields):
        raise ValueError(f"checksum manifest lacks required columns: {sorted(required - fields)}")
    return rows


def validate_fasta(data: bytes, expected_species: str) -> tuple[int, set[tuple[str, str]]]:
    count = 0
    seen_ids: set[str] = set()
    seen_keys: set[tuple[str, str]] = set()
    header = None
    sequence: list[str] = []

    def check_record() -> None:
        nonlocal count
        if header is None:
            return
        fields = header.split("|")
        if len(fields) != 3:
            raise ValueError(f"wrong FASTA schema for {header!r}")
        species, accession, _gene = fields
        if species != expected_species:
            raise ValueError(f"species prefix {species!r} does not match file code {expected_species!r}")
        if not ACCESSION_RE.fullmatch(accession):
            raise ValueError(f"invalid UniProt accession syntax in {header!r}")
        if header in seen_ids or (species, accession) in seen_keys:
            raise ValueError(f"duplicate FASTA ID or species/accession: {header}")
        joined = "".join(sequence).upper()
        if not joined:
            raise ValueError(f"empty sequence for {header}")
        unexpected = sorted(set(joined) - STANDARD_AA)
        if unexpected:
            raise ValueError(f"unexpected amino-acid characters for {header}: {''.join(unexpected)}")
        seen_ids.add(header)
        seen_keys.add((species, accession))
        count += 1

    for line_number, raw in enumerate(io.StringIO(data.decode("utf-8")), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            check_record()
            header = line[1:].split()[0]
            sequence = []
        else:
            if header is None:
                raise ValueError(f"sequence before first FASTA header on line {line_number}")
            sequence.append(line)
    check_record()
    return count, seen_keys


def validate_lookup(data: bytes) -> tuple[int, set[tuple[str, str]]]:
    text = data.decode("utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {
        "species_short_code",
        "uniprot_accession",
        "gene_symbol",
        "sequence_length",
        "sequence_sha256",
    }
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"protein lookup lacks columns: {sorted(required - set(reader.fieldnames or []))}")
    keys: set[tuple[str, str]] = set()
    count = 0
    for row in reader:
        key = (row["species_short_code"], row["uniprot_accession"])
        if key in keys:
            raise ValueError(f"duplicate species/accession in protein lookup: {key}")
        keys.add(key)
        count += 1
    return count, keys


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify Lab 3 data checksums, IDs, and sequence content.")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    repo = args.repo.resolve()

    try:
        manifest = load_manifest(repo / "checksums" / "package_checksums.tsv")
        all_fasta_keys: set[tuple[str, str]] = set()
        lookup_keys: set[tuple[str, str]] | None = None
        total_proteins = 0

        for row in manifest:
            path = repo / row["packaged_file"]
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"missing or empty packaged file: {path}")
            packaged = path.read_bytes()
            if sha256_bytes(packaged) != row["packaged_sha256"]:
                raise ValueError(f"packaged checksum mismatch: {path}")
            data = gzip.decompress(packaged) if path.suffix == ".gz" else packaged
            if sha256_bytes(data) != row["decompressed_sha256"]:
                raise ValueError(f"content checksum mismatch: {path}")
            if len(data) != int(row["decompressed_bytes"]):
                raise ValueError(f"decompressed byte count mismatch: {path}")

            if row["record_type"] == "protein_fasta":
                species_code = path.name.split(".")[0]
                count, keys = validate_fasta(data, species_code)
                if keys & all_fasta_keys:
                    raise ValueError(f"duplicate accession across packaged proteomes: {path}")
                all_fasta_keys.update(keys)
                total_proteins += count
            elif row["record_type"] == "protein_lookup_tsv":
                count, lookup_keys = validate_lookup(data)
            else:
                raise ValueError(f"unknown record type in checksum manifest: {row['record_type']}")
            if count != int(row["record_count"]):
                raise ValueError(f"record-count mismatch: {path}")

        if lookup_keys != all_fasta_keys:
            missing_from_lookup = all_fasta_keys - (lookup_keys or set())
            missing_from_fasta = (lookup_keys or set()) - all_fasta_keys
            raise ValueError(
                f"FASTA/lookup key mismatch: {len(missing_from_lookup)} missing from lookup; "
                f"{len(missing_from_fasta)} missing from FASTA"
            )

        species_path = repo / "metadata" / "species_key.tsv"
        roster_path = repo / "metadata" / "starting_proteins.tsv"
        with species_path.open(encoding="utf-8", newline="") as handle:
            species_rows = list(csv.DictReader(handle, delimiter="\t"))
        species_codes = [row["species_short_code"] for row in species_rows]
        if len(species_codes) != len(set(species_codes)):
            raise ValueError("duplicate species code in species key")
        packaged_species = {path.name.split(".")[0] for path in (repo / "proteomes").glob("*.canonical.faa.gz")}
        if set(species_codes) != packaged_species:
            raise ValueError("species key does not match packaged proteomes")

        with roster_path.open(encoding="utf-8", newline="") as handle:
            roster_rows = list(csv.DictReader(handle, delimiter="\t"))
        projects = [row["project_id"] for row in roster_rows]
        accessions = [row["human_uniprot_accession"] for row in roster_rows]
        if len(projects) != len(set(projects)) or len(accessions) != len(set(accessions)):
            raise ValueError("duplicate project ID or human accession in starting-protein roster")
        missing_queries = {("Hsap", accession) for accession in accessions} - all_fasta_keys
        if missing_queries:
            raise ValueError(f"approved query accessions missing from human proteome: {sorted(missing_queries)}")
    except (KeyError, OSError, UnicodeError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error

    print(f"PASS: {len(species_codes)} proteomes; {total_proteins} protein records; {len(roster_rows)} approved projects")
    print("PASS: packaged and decompressed checksums, identifier schemas, metadata joins, and amino-acid alphabet")


if __name__ == "__main__":
    main()
