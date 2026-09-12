# BIO 312 Lab 3: Finding Protein Homologs with BLAST

**BIO 312-01 — Bioinformatics and Computational Biology — Fall 2026 (course 84182)**

Copyright Joshua Rest. Do not post or redistribute this lab outside the course.

## Learning objectives

By the end of this lab, you should be able to:

- connect a UniProt accession to an amino-acid FASTA sequence;
- explain what a proteome and a local protein BLAST database contain;
- combine proteomes and build a protein BLAST database;
- retrieve one protein sequence with `samtools faidx`;
- use `blastp` to find proteins similar to a query sequence;
- examine and filter BLAST results with standard command-line tools;
- adapt the myoglobin workflow to your assigned protein; and
- prepare the protein sequences that you will align in Lab 4.

## The semester project question

When mammals independently evolved aquatic lifestyles, were changes in protein sequence or predicted structure associated with these transitions, or were important features conserved?

The species panel includes three independent aquatic origins: cetaceans, pinnipeds, and sirenians. You will first investigate myoglobin and then begin the same investigation for your assigned protein.

## From nucleotide records to protein comparisons

In Labs 1 and 2, we used NCBI RefSeq records displayed in GenBank format to study genomic DNA, transcripts, coding sequences, and the resulting protein. Beginning in Lab 3, we will use UniProtKB and organismal proteomes because our analyses now focus on amino-acid sequences and protein evolution.

GenBank is an archival sequence database that is especially useful for nucleotide, mRNA, and protein records. UniProtKB is protein-centered: its records organize amino-acid sequences, protein names, functional annotations, domains, structures, locations, evidence, and links to other databases. Reviewed UniProtKB/Swiss-Prot entries have been manually curated; unreviewed UniProtKB/TrEMBL entries are computationally annotated.

The human myoglobin protein used in Labs 1 and 2 has RefSeq accession `NP_001349775.1`; its corresponding reviewed UniProtKB accession is `P02144`. We will use the UniProt accession for Labs 3–7. **BLAST is an analysis tool** that compares sequences.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A01. (2 points)** Which statements accurately describe the transition from Labs 1–2 to Lab 3? Select all that apply.

## What to complete

- **Brightspace assignment: 15 points.** Answer `L3-A01` through `L3-A10`. You have two attempts; Brightspace records the higher score.
- **GitHub repository check: 10 points.** Complete every [GitHub] answer block and commit the requested result files. Answer-block completion and generated-file completion are scored separately, so partial credit is possible. A detailed rubric is posted with Lab 3 in Brightspace.
- **Exit ticket: 10 points.** Complete the separate paper exit ticket during lab. Devices, notes, AI tools, and collaboration are not permitted during the exit ticket.
- **Project literature search.** Begin the separate Brightspace assignment during lab. It is due at the start of your next laboratory meeting; see Part VII for how to begin.

Lab work is due 24 hours before the next laboratory meeting.

## Part I. Reconnect and clone Lab 3

### Reconnect to the existing course instance

1. Start the AWS Academy Learner Lab session.
2. Wait while the existing EC2 instance starts. In the AWS console, confirm that it is **Running** and that both status checks pass.
3. Open VS Code and reconnect using the saved Remote–SSH entry from Lab 2.
4. If prompted, enter the course Linux password `bio312isCOOL!` for username `bio312-user`.
5. In the remote VS Code window, select **Terminal → New Terminal**.

Copy each command from a `bash` block into the VS Code terminal unless these instructions tell you to replace a value.

Confirm the Linux account:

```bash
whoami
```

Confirm the current directory:

```bash
pwd
```

**Expected output:** `bio312-user` and `/home/bio312-user`.

If the connection fails, first confirm that the AWS Academy Learner Lab has started and that the instance and both status checks have finished starting. After connecting, use `pwd` to show your location and `ls` to inspect its contents before repeating a command.

### Clone through VS Code

1. Open your assigned Lab 3 repository on GitHub. Its name follows the pattern `bio312-lab03-YOUR-GITHUB-USERNAME`.
2. Copy the repository URL from your browser.
3. In the remote VS Code window, open **View → Command Palette → Git: Clone**.
4. Paste the URL, complete VS Code's GitHub browser sign-in if requested, and choose `/home/bio312-user` as the parent folder.
5. Select **Open** when cloning finishes.
6. In the newly opened repository window, select **Terminal → New Terminal** again.

Make all edits in the remote VS Code copy. Editing the same file on GitHub's website and in VS Code at the same time can create competing versions that Git must reconcile.

### Save your GitHub username for this and later labs

A **shell variable** gives a name to a piece of text. In the commands below, `$MYGIT` tells the shell to insert your GitHub username. Adding an `export` line to `.bashrc` sets the variable each time you open a terminal on this instance.

Open `/home/bio312-user/.bashrc` in VS Code. Add the following line at the end, replacing `your-github-username` with your own GitHub username. If a line beginning with `export MYGIT=` is already present, update that line instead.

```text
export MYGIT="your-github-username"
```

Save `.bashrc`, close the terminal with the trash-can button, and select **Terminal → New Terminal**. Check the saved value:

```bash
echo "$MYGIT"
```

### Move around the directory system

The Explorer panel in VS Code lets you open files, but it does not change the terminal's location. Three commands help you orient yourself in the terminal: `cd` changes the current directory, `pwd` prints the current directory, and `ls` lists the current directory's contents.

Move into the Lab 3 repository with `cd`:

```bash
cd /home/bio312-user/bio312-lab03-$MYGIT
```

Confirm the location with `pwd`:

```bash
pwd
```

**Expected output:** `/home/bio312-user/bio312-lab03-YOUR-GITHUB-USERNAME`

List the files and directories in your current location with `ls`:

```bash
ls
```

**Expected output:** the listing should include `README.md`, `checksums`, `img`, `metadata`, `proteomes`, and `scripts`.

Open `README.md` in VS Code. Enter requested responses only inside the labeled `notes/answer` blocks, and be sure to save the file often.

## Part II. Explore the protein dataset

### Check the supplied files

Run the package check:

```bash
python3 /home/bio312-user/bio312-lab03-$MYGIT/scripts/validate_package.py --repo /home/bio312-user/bio312-lab03-$MYGIT
```

Read both lines beginning with `PASS`. They report how many proteomes, protein records, and assigned projects are included.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A02. (1 point)** According to the package check, how many proteomes and protein sequences will be included in the BLAST database? What does each number count?

List the compressed proteome files:

```bash
ls /home/bio312-user/bio312-lab03-$MYGIT/proteomes
```

Look inside the compressed human proteome. `gzip -cd` displays decompressed text without replacing the `.gz` file. The `|` symbol is a **pipe**: it sends the output from the first command into the next command. Here, `less` lets you move through that output.

```bash
gzip -cd /home/bio312-user/bio312-lab03-$MYGIT/proteomes/Hsap.canonical.faa.gz | less
```

Use the arrow keys or spacebar to move through the file. Press `q` to leave `less`.

### Species panel and aquatic origins

The complete species information is in `metadata/species_key.tsv`.

| Code | Species | Common name | Category | Aquatic origin or comparison |
|---|---|---|---|---|
| `Hsap` | *Homo sapiens* | Human | terrestrial | contextual taxon |
| `Pmac` | *Physeter macrocephalus* | Sperm whale | aquatic | cetacean |
| `Bmus` | *Balaenoptera musculus* | Blue whale | aquatic | cetacean |
| `Btau` | *Bos taurus* | Cattle | terrestrial | cetartiodactyl comparison |
| `Sscr` | *Sus scrofa* | Pig | terrestrial | cetartiodactyl comparison |
| `Oros` | *Odobenus rosmarus divergens* | Pacific walrus | aquatic | pinniped |
| `Zcal` | *Zalophus californianus* | California sea lion | aquatic | pinniped |
| `Clup` | *Canis lupus familiaris* | Dog | terrestrial | carnivoran comparison |
| `Mfur` | *Mustela putorius furo* | Domestic ferret | terrestrial | carnivoran comparison |
| `Tman` | *Trichechus manatus latirostris* | Florida manatee | aquatic | sirenian |
| `Lafr` | *Loxodonta africana* | African elephant | terrestrial | paenungulate comparison |

The Mirceta et al. reading introduced comparisons among independently evolved diving mammals. For this project, use the supplied `independent_aquatic_transition` and `terrestrial_comparison_lineage` columns rather than trying to reconstruct all of mammalian taxonomy yourself.

View those columns and the complete species information:

```bash
column -t -s $'\t' /home/bio312-user/bio312-lab03-$MYGIT/metadata/species_key.tsv | less -S
```

Walruses and sea lions forage in water and haul out on land. For this introductory comparison, the supplied table places these diving pinnipeds in the aquatic category.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A03. (1 point)** Use the supplied species table and the Mirceta et al. reading: which three independent aquatic origins are represented?

### Protein identifiers

A **proteome** is the set of proteins represented for an organism. These files contain canonical sequences from 11 UniProt reference proteomes. A canonical sequence is UniProt's main representative sequence for a protein entry; other isoforms may exist.

The supplied FASTA headers have three pipe-separated fields:

```text
>SpeciesCode|UniProtAccession|GeneSymbol
```

For example:

```text
>Hsap|P02144|MB
```

- `Hsap` is the course species code.
- `P02144` is the UniProt accession for the protein record.
- `MB` is the gene symbol associated with that record.

A **key** is a field used to identify one record uniquely. The UniProt accession is the protein-record key in this project. A gene symbol is useful biological information, but the same symbol appears in multiple species.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A04. (1 point)** In `Hsap|P02144|MB`, which field is the protein-record key, and why is the gene symbol not unique in this dataset?

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A05. (1 point)** Which statement correctly describes a proteome, a canonical sequence, and an isoform?

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G01. (0.5 point)** Record the human myoglobin RefSeq accession, UniProt accession, UniProt review status, and amino-acid sequence length.

```notes/answer
myoglobin Ref: NP_001349775.1
Uni prot: P02144
reviewed and swiss verified
154 amino acids
```

## Part III. Combine the proteomes and build a BLAST database

BLAST will search all 11 proteomes together. Because the supplied files are compressed, `gzip -cd` reads and decompresses them in the order shown inside the braces. Unlike the pipe you just used, the `>` symbol saves the combined output in a new FASTA file.

Combine the proteomes:

```bash
gzip -cd /home/bio312-user/bio312-lab03-$MYGIT/proteomes/{Hsap,Pmac,Bmus,Btau,Sscr,Oros,Zcal,Clup,Mfur,Tman,Lafr}.canonical.faa.gz > /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta
```

Count the protein sequences in the combined file:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta
```

Look at the beginning of the combined FASTA:

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta
```

Press `q` to leave `less`.

Create a directory for the BLAST database files:

```bash
mkdir -p /home/bio312-user/bio312-lab03-$MYGIT/blast_db
```

`makeblastdb` turns the combined FASTA into a database that BLAST can search. `-in` names the input FASTA, `-dbtype prot` specifies proteins, `-parse_seqids` retains the FASTA identifiers, and `-out` sets the database prefix.

Build the protein database:

```bash
makeblastdb -in /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta -dbtype prot -parse_seqids -out /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein
```

List the database files created by BLAST+:

```bash
ls /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein.*
```

`allprotein.fasta` is the sequence collection used to build the database. `blast_db/allprotein` is the shared prefix of the database files that `blastp` searches.

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G02. (0.5 point)** Record the combined FASTA sequence count, its first FASTA identifier, and the BLAST database filename extensions.

```notes/answer
218233
Hsap|A0A075B6H9|IGLV4-69
extensions: .pdb,.phr,.pin,.pjs,.pog,.pos,.pot,.psq,.ptf,.pto
```

## Part IV. Use BLAST to find proteins similar to human myoglobin

In BLAST terminology, the sequence you start with is the **query**, and the sequence collection being searched is the **database**. For our shared example, the query is the human myoglobin sequence in `myoglobin/P02144.fasta`. BLASTP will compare it with every protein in the database built from the 11 proteomes and report similar sequences. These matches are the starting set of possible myoglobin homologs that we will examine further.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A06. (1 point)** Which file contains the query sequence, and what collection does BLASTP search as the database?

### Retrieve the human myoglobin query with samtools

`samtools faidx` is a standard bioinformatics tool for quickly retrieving named sequences from a FASTA file. First, index the combined FASTA:

```bash
samtools faidx /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta
```

Create the myoglobin directory:

```bash
mkdir -p /home/bio312-user/bio312-lab03-$MYGIT/myoglobin
```

Retrieve the exact human myoglobin identifier. The quotes keep the pipe characters inside the identifier.

```bash
samtools faidx /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta 'Hsap|P02144|MB' > /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/P02144.fasta
```

Check that the new file contains one FASTA sequence:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/P02144.fasta
```

Look at the query sequence:

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/P02144.fasta
```

### Search for similar proteins and read the alignments

`blastp` compares a protein query with a protein database. This first search produces the familiar, readable BLAST report with a ranked hit list followed by pairwise alignments.

Important options are:

- `-query`: the human myoglobin FASTA;
- `-db`: the database prefix;
- `-evalue 1e-10`: report matches meeting this search threshold;
- `-max_target_seqs 5000`: allow a large result list;
- `-max_hsps 1`: report the best aligned segment for each query–subject pair;
- `-num_threads 2`: use two CPU threads;
- `-outfmt 0`: request the readable report; and
- `-out`: name the output file.

Run BLASTP:

```bash
blastp -query /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/P02144.fasta -db /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein -evalue 1e-10 -max_target_seqs 5000 -max_hsps 1 -num_threads 2 -outfmt 0 -out /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.txt
```

Examine the ranked hits and several alignments:

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.txt
```

The first result is the **self-hit**: the human query finds the identical human record in the database. The next result is the highest-scoring non-self hit. Compare their identifiers, percent identities, alignment lengths, and E-values.

An E-value estimates how many matches with this score or better would be expected by chance in a search of this size. Smaller values are stronger evidence that a match is not due to chance. BLAST may display an extremely small E-value as `0.0` because of output precision.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A07. (2 points)** Using the readable myoglobin BLAST report, identify the self-hit and the highest-scoring non-self hit. Compare their displayed percent identities and E-values, and use the query coordinates to decide whether each alignment covers the full 154-amino-acid query.

### Request tabular output from the same search

The readable report is useful for examining alignments. A tabular report is easier to count and filter. We will repeat the same BLAST search but change `-outfmt` and the output filename.

The tabular output has these nine columns:

| Column | Field | Meaning |
|---:|---|---|
| 1 | `qseqid` | query sequence ID |
| 2 | `sseqid` | matching database sequence ID |
| 3 | `pident` | percent identical amino acids in the alignment |
| 4 | `length` | alignment length |
| 5 | `qlen` | complete query length |
| 6 | `slen` | complete subject length |
| 7 | `evalue` | expected number of chance matches at least this strong |
| 8 | `bitscore` | normalized alignment score; larger is stronger |
| 9 | `qcovs` | percentage of the query covered by the alignment |

Run BLASTP and request the tabular output:

```bash
blastp -query /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/P02144.fasta -db /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein -evalue 1e-10 -max_target_seqs 5000 -max_hsps 1 -num_threads 2 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore qcovs" -out /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.tsv
```

Inspect the first five rows:

```bash
head -n 5 /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.tsv
```

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G03. (0.5 point)** Paste the complete tabular rows for the human self-hit and the highest-scoring non-self hit.

```notes/answer
Hsap|P02144|MB  Hsap|P02144|MB  100.000 154     154     154     1.50e-110       312     100
Hsap|P02144|MB  Sscr|P02189|MB  93.506  154     154     154     4.04e-103       294     100
Hsap|P02144|MB  Tman|A0A2Y9DG99|MB      90.260  154     154     154     1.93e-98        282     100
Hsap|P02144|MB  Oros|A0A2U3X238|MB      88.312  154     154     154     6.54e-98        281     100
Hsap|P02144|MB  Zcal|P02161|MB  87.662  154     154     154     3.21e-97        279     100
(base) [bio312-user@ip-172-31-82-182 bio312-lab03-shrikantgoyal-eng]$ 
```

### Filter the BLAST hits

First, count how many hits BLAST reported before filtering:

```bash
wc -l /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.tsv
```

We are going to filter the hits to meet the following conditions:

- E-value at or below `1e-10`;
- percent identity at least `35`; and
- query coverage at least `70` percent.

In the tabular output, these are columns 7, 3, and 9. `awk` tests each row and writes the rows that meet all three conditions to a new file.

```bash
awk '$7 <= 1e-10 && $3 >= 35 && $9 >= 70' /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.blastp.tsv > /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidates.tsv
```

Count the hits after filtering:

```bash
wc -l /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidates.tsv
```

Column 2 contains the exact identifier for each matching protein. Extract those identifiers and use `sort -u` to make a unique list:

```bash
awk '{print $2}' /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidates.tsv | sort -u > /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidate_ids.txt
```

Use `seqkit grep`, a standard FASTA-processing tool, to retrieve those proteins from the combined sequence file:

```bash
seqkit grep -f /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidate_ids.txt /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta > /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.selected_source.fasta
```

Check that the number of retrieved sequences equals the number of identifiers:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.selected_source.fasta
```

Look at the retrieved sequences:

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.selected_source.fasta
```

The next script performs the course-specific step that standard FASTA tools do not know how to do: it adds the project and taxonomy fields needed by later labs, checks the metadata, and creates the species summary.

```bash
python3 /home/bio312-user/bio312-lab03-$MYGIT/scripts/prepare_lab4_fasta.py --project-id MB --candidate-table /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidates.tsv --candidate-ids /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.candidate_ids.txt --source-fasta /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.selected_source.fasta --species-key /home/bio312-user/bio312-lab03-$MYGIT/metadata/species_key.tsv --protein-lookup /home/bio312-user/bio312-lab03-$MYGIT/metadata/protein_lookup.tsv.gz --output-dir /home/bio312-user/bio312-lab03-$MYGIT/myoglobin --basename myoglobin
```

Inspect the number of selected proteins from each species:

```bash
column -t -s $'\t' /home/bio312-user/bio312-lab03-$MYGIT/myoglobin/myoglobin.species_counts.tsv
```

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G04. (0.5 point)** Record the number of myoglobin hits before and after filtering and the count for each species.

```notes/answer
41->12
species_short_code  common_name          candidate_protein_records
Hsap                Human                1
Pmac                Sperm whale          1
Bmus                Blue whale           1
Btau                Cattle               1
Sscr                Pig                  1
Oros                Pacific walrus       1
Zcal                California sea lion  1
Clup                Dog                  1
Mfur                Domestic ferret      1
Tman                Florida manatee      1
Lafr                African elephant     2
```

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A08. (2 points)** How many myoglobin hits were present before and after filtering? How many were removed, and did every species retain at least one hit?

## Part V. Begin your assigned protein investigation

Your assigned **project ID** appears in Brightspace. For this course it is the gene symbol of your human starting protein and the short name used for your project directory and files. The UniProt accession identifies the particular protein record.

The shell variables `$PROJECT_ID` and `$ACCESSION` work like `$MYGIT`: whenever the shell sees one in a command, it substitutes the value you saved. In the commands you adapt below, you may use these variables or type your actual project ID and accession. If you type an actual value, omit the `$` and use that value consistently. The variables simply reduce repeated typing; they do not change what a command does.

### Find and save your project information

Open `/home/bio312-user/.bashrc` in VS Code. Add the following line, replacing `YOUR_PROJECT_ID` with the project ID assigned to you. Preserve capitalization.

```text
export PROJECT_ID="YOUR_PROJECT_ID"
```

Save `.bashrc`, close the terminal with the trash-can button, and select **Terminal → New Terminal**. Check the value:

```bash
echo "$PROJECT_ID"
```

Find the information for your assigned protein. This `awk` command displays the header plus the row whose first column exactly matches your project ID.

```bash
awk -F $'\t' -v project="$PROJECT_ID" 'NR==1 || $1==project' /home/bio312-user/bio312-lab03-$MYGIT/metadata/starting_proteins.tsv
```

Read the UniProt accession in the third column. Return to `/home/bio312-user/.bashrc` and add the following line, replacing the example with your assigned accession:

```text
export ACCESSION="YOUR_UNIPROT_ACCESSION"
```

Save `.bashrc`, open a new terminal, and check all three variables:

```bash
echo "$MYGIT $PROJECT_ID $ACCESSION"
```

Search for the accession at [UniProtKB](https://www.uniprot.org/uniprotkb/). Read the protein name and function sections and look for any family or related-protein names mentioned in the record.

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G05. (0.5 point)** Record your project ID, UniProt accession, protein name, review status, and sequence length.

```notes/answer
ID: ASS1
Acession: P00966
Name: Argininosuccinate synthase
Reviewed(Swiss prot)
length: 412 amino acids
```

### Retrieve your query protein

Create a directory named with your project ID:

```bash
mkdir -p /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID
```

Use `samtools faidx` to retrieve the exact human protein. Read the identifier inside the quotes and make sure you understand where each variable will be inserted.

```bash
samtools faidx /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta "Hsap|$ACCESSION|$PROJECT_ID" > /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/query.fasta
```

Confirm that you retrieved one sequence, then examine it with `less`:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/query.fasta
```

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/query.fasta
```

### Follow the files through your project

Each main output becomes an input to a later step:

`query.fasta` → `.blastp.tsv` → `.candidates.tsv` → `.candidate_ids.txt` → `.selected_source.fasta` → `.homologs.fas`

The readable `.blastp.txt` report is for examining alignments. As you adapt each command, make sure its input filename matches the file created in the preceding step.

### Decide what must change for your BLAST search

Return to the readable and tabular myoglobin `blastp` commands in Part IV. For your assigned protein:

- the query path must point to your new `query.fasta`;
- the two output paths must use your project directory and meaningful filenames; and
- the database and other search settings stay the same.

Later labs will ask you to make more of these command adaptations yourself. For now, use the myoglobin commands as a model and pay attention to which parts name inputs, outputs, and search settings.

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G06, part 1. (0.5 point)** Before running the commands, explain which parts of the two myoglobin BLAST commands you need to change, which parts remain the same, and how your filenames avoid overwriting the myoglobin results.

```notes/answer
The query now points to the ASS1 fasta. The output command will be different to accomidate different rows to be readable and the tabluar report. The database used is the same because its searching through all 11 proteomes, the e value and command output format is the same.Saving the new results in projects/ASS1 with ASS1 filenames keeps them separate from the results in myoglobin and avoids overwriting them.
```

Copy both myoglobin BLAST commands into the next answer block and adapt them. Use these two output filenames:

- `projects/$PROJECT_ID/$PROJECT_ID.blastp.txt`
- `projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv`

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G06, part 2. (0.5 point)** Paste your two adapted `blastp` commands here. Save `README.md`, then copy each command from your answer block into the terminal and run it.

```notes/answer
blastp -query /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/query.fasta -db /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein -evalue 1e-10 -max_target_seqs 5000 -max_hsps 1 -num_threads 2 -outfmt 0 -out /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.txt

blastp -query /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/query.fasta -db /home/bio312-user/bio312-lab03-$MYGIT/blast_db/allprotein -evalue 1e-10 -max_target_seqs 5000 -max_hsps 1 -num_threads 2 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore qcovs" -out /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv
```

Confirm that both output files exist and contain results:

```bash
ls -lh /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.txt /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv
```

Count the hits before filtering:

```bash
wc -l /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv
```

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G07, part 1. (0.5 point)** Record the number of BLAST hits before filtering.

```notes/answer
13
```

Examine the readable report. Look at the self-hit, the highest-scoring match from another mammal, and how the matches change farther down the list.

```bash
less /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.txt
```

### Adapt the filtering and sequence-retrieval commands

Return to the two `awk` commands used for myoglobin. Copy them into the next answer block and change the input and output paths for your project. They should create:

- `projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv`; and
- `projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt`.

The filtering values and column numbers remain the same.

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G07, part 2. (0.5 point)** Paste your two adapted `awk` commands here. Save `README.md`, then copy each command from your answer block into the terminal and run it.

```notes/answer
awk '$7 <= 1e-10 && $3 >= 35 && $9 >= 70' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv > /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv

awk '{print $2}' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv | sort -u > /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt
```

Count the hits after filtering:

```bash
wc -l /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv
```

Use `seqkit grep` to retrieve the proteins named in your ID list:

```bash
seqkit grep -f /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt /home/bio312-user/bio312-lab03-$MYGIT/allprotein.fasta > /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.selected_source.fasta
```

Check that the number of retrieved sequences matches the filtered hit count:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.selected_source.fasta
```

Prepare the Lab 4 FASTA and species summary:

```bash
python3 /home/bio312-user/bio312-lab03-$MYGIT/scripts/prepare_lab4_fasta.py --project-id "$PROJECT_ID" --candidate-table /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv --candidate-ids /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt --source-fasta /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.selected_source.fasta --species-key /home/bio312-user/bio312-lab03-$MYGIT/metadata/species_key.tsv --protein-lookup /home/bio312-user/bio312-lab03-$MYGIT/metadata/protein_lookup.tsv.gz --output-dir /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID --basename "$PROJECT_ID"
```

Inspect the count for each species:

```bash
column -t -s $'\t' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.species_counts.tsv
```

Summarize the gene symbols in the filtered set:

```bash
cut -d '|' -f 3 /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt | sort | uniq -c | sort -nr
```

Choose the highest-scoring non-self hit in your tabular output. Also consider what the UniProt record and your BLAST results have shown you about the protein or its family.

> <img src="img/github.png" alt="GitHub notes" width="20" height="20"> [GitHub] **L3-G07, part 3. (0.5 point)** Record the number of hits after filtering, how many were removed, all 11 species counts, and the gene-symbol summary. Paste the complete tabular row for the highest-scoring non-self hit. Then, in one or two sentences, describe something you learned about the protein or its family from the UniProt record and your BLAST results.

```notes/answer
11 hits after filter
species_short_code  common_name          candidate_protein_records
Hsap                Human                1
Pmac                Sperm whale          0
Bmus                Blue whale           1
Btau                Cattle               2
Sscr                Pig                  1
Oros                Pacific walrus       1
Zcal                California sea lion  1
Clup                Dog                  1
Mfur                Domestic ferret      1
Tman                Florida manatee      1
Lafr                African elephant     1

11 ASS1
Highest-scoring non-self hit:
Hsap|P00966|ASS1	Mfur|M3YMP6|ASS1	98.301	412	412	412	0.0	847	100

My filtered BLAST results contain 11 candidate proteins across 10 species, all labeled ASS1. No sperm whale record passed the screening, but this alone does not show that the species has lost the gene.
```

### TA checkpoint: confirm that your project is ready for Lab 4

Count the filtered rows again:

```bash
wc -l /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv
```

Count the sequences in the FASTA prepared for Lab 4:

```bash
grep -c '^>' /home/bio312-user/bio312-lab03-$MYGIT/projects/$PROJECT_ID/$PROJECT_ID.homologs.fas
```

Before beginning the exit ticket, show a TA these two files and the two counts. Both files should be nonempty, and the number of filtered rows should equal the number of FASTA sequences. This brief checkpoint confirms that your assigned-protein files are ready for Lab 4; it is not a separate submission.

### What has BLAST told you?

Your filtered result is the set of candidate homologs that you will examine in the next labs. In Lab 4, you will align their amino-acid sequences. In Lab 5, you will build a protein-family tree and investigate how the sequences are related.

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A09. (2 points)** How many candidate homologs remain in your assigned-protein results after filtering? Enter the number of rows in your `$PROJECT_ID.candidates.tsv` file.

The final Brightspace question asks you to apply the same three filtering conditions to a short example.

| Hit | E-value | Percent identity | Query coverage |
|---|---:|---:|---:|
| A | `6e-18` | 42% | 82% |
| B | `3e-50` | 72% | 65% |
| C | `2e-40` | 33% | 95% |
| D | `4e-8` | 49% | 90% |

> <img src="img/brightspace.png" alt="Brightspace question" width="20" height="20"> [Brightspace] **L3-A10. (2 points)** We retain hits with E-value ≤ `1e-10`, percent identity ≥ 35%, and query coverage ≥ 70%. Which hit in the table would be retained?

If your query file or either BLAST output is empty, show the TA the command and files at that step before continuing.

## Part VI. Files to keep for Lab 4

The table below lists the files created for myoglobin and your assigned protein. The `.homologs.fas` files are the sequence files that Lab 4 will align.

| Product | Myoglobin path | Assigned-protein path |
|---|---|---|
| Query FASTA | `myoglobin/P02144.fasta` | `projects/$PROJECT_ID/query.fasta` |
| Readable BLAST report | `myoglobin/myoglobin.blastp.txt` | `projects/$PROJECT_ID/$PROJECT_ID.blastp.txt` |
| Original tabular BLAST results | `myoglobin/myoglobin.blastp.tsv` | `projects/$PROJECT_ID/$PROJECT_ID.blastp.tsv` |
| Filtered BLAST rows | `myoglobin/myoglobin.candidates.tsv` | `projects/$PROJECT_ID/$PROJECT_ID.candidates.tsv` |
| Filtered exact IDs | `myoglobin/myoglobin.candidate_ids.txt` | `projects/$PROJECT_ID/$PROJECT_ID.candidate_ids.txt` |
| Retrieved source sequences | `myoglobin/myoglobin.selected_source.fasta` | `projects/$PROJECT_ID/$PROJECT_ID.selected_source.fasta` |
| Species counts | `myoglobin/myoglobin.species_counts.tsv` | `projects/$PROJECT_ID/$PROJECT_ID.species_counts.tsv` |
| Identifier mapping | `myoglobin/myoglobin.identifier_mapping.tsv` | `projects/$PROJECT_ID/$PROJECT_ID.identifier_mapping.tsv` |
| **Lab 4 sequence input** | `myoglobin/myoglobin.homologs.fas` | `projects/$PROJECT_ID/$PROJECT_ID.homologs.fas` |

`allprotein.fasta`, its `.fai` index, and the `blast_db/` directory can be regenerated, so they do not need to be committed.

## Part VII. Start your project literature search

Open the separate project literature-search assignment in Brightspace and follow its complete directions. Unlike the Lab 3 work above, that assignment is due at the start of your next laboratory meeting; verify the exact time in Brightspace.

Before leaving lab:

1. Search [Google Scholar](https://scholar.google.com/) using your project ID, the full protein name from UniProt, and one or two terms related to its function or evolution.
2. Identify at least three promising scholarly sources and save their complete citation information.
3. Read each abstract—and inspect the paper when available—well enough to decide how it could help you understand your protein. Do not rely only on a title, search snippet, or AI-generated summary.

The Brightspace assignment asks for brief relevance explanations for two sources and a more detailed annotation for the third. An AI-assisted literature-search tool may help you discover possible sources, but verify each citation in Google Scholar, PubMed, or the publisher's page before using it.

## Part VIII. Save with Git and GitHub

Save `README.md` in VS Code. Move to the Lab 3 repository:

```bash
cd /home/bio312-user/bio312-lab03-$MYGIT
```

Inspect changed and new files:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT status --short
```

Inspect your README answers before staging:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT diff -- README.md
```

Stage the README and the two result directories:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT add README.md myoglobin projects/$PROJECT_ID
```

Confirm what will be committed. `allprotein.fasta`, `allprotein.fasta.fai`, and `blast_db/` should not appear:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT status --short
```

Commit locally:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT commit -m "Complete Lab 3 protein BLAST analyses"
```

Push the commit to GitHub:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT push
```

If VS Code asks you to sign in to GitHub, complete the browser sign-in and run `git push` again. Do not enter passwords, tokens, private keys, IP addresses, or AWS credentials in any submitted file.

Verify the latest local commit:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT log -1 --oneline
```

Confirm a clean local working tree:

```bash
git -C /home/bio312-user/bio312-lab03-$MYGIT status
```

**Expected output:** `nothing to commit, working tree clean`

Refresh your repository on GitHub and confirm that the latest commit and result files are visible.

Finally, close the remote VS Code window and select **End Lab** in AWS Academy. The same saved VS Code connection can be reused next time.

## Final checklist

- [ ] I answered `L3-A01`–`L3-A10` in the 15-point Brightspace assignment.
- [ ] I completed the separate 10-point paper exit ticket during lab.
- [ ] I completed every [GitHub] block for the 10-point repository check and saved `README.md`.
- [ ] I combined 11 proteomes, examined the FASTA files, and built the BLAST database with `makeblastdb`.
- [ ] I retrieved both human queries with `samtools faidx`.
- [ ] I ran readable and tabular BLAST searches for myoglobin and my assigned protein.
- [ ] I counted hits before and after filtering and retrieved the filtered sequences with `seqkit grep`.
- [ ] My myoglobin and assigned-protein directories contain all files listed in Part VI.
- [ ] Both `.homologs.fas` files are ready for Lab 4.
- [ ] A TA checked that my assigned `.candidates.tsv` and `.homologs.fas` files are nonempty and contain the same number of records.
- [ ] I opened the separate literature-search assignment and saved citation information for at least three promising sources.
- [ ] The latest commit is visible on GitHub and `git status` reports a clean working tree.
- [ ] I ended the AWS Academy Learner Lab session.

## Current technical references

- [NCBI BLAST+ command-line manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
- [NCBI: building a local BLAST database](https://www.ncbi.nlm.nih.gov/books/NBK569841/)
- [Samtools faidx manual](https://www.htslib.org/doc/samtools-faidx.html)
- [SeqKit grep documentation](https://bioinf.shenwei.me/seqkit/usage/#grep)
- [UniProtKB help](https://www.uniprot.org/help/uniprotkb)
- [VS Code: working with GitHub](https://code.visualstudio.com/docs/sourcecontrol/github)
