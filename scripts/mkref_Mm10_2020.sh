#!/bin/bash
# Genome metadata
genome="mm10"
version="2020-A"

# Set up source and build directories
build="mm10-2020-A_build"
mkdir -p "$build"

# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"

fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.vM23.primary_assembly.annotation.gtf"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi

# Modify sequence headers in the Ensembl FASTA to match GENCODE format
fasta_modified="$build/$(basename "$fasta_in").modified"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from Ensembl gene, transcript, and exon IDs
gtf_modified="$build/$(basename "$gtf_in").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for filtering biotypes
BIOTYPE_PATTERN="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Construct the gene ID allowlist
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$gtf_modified" >> "$gtf_filtered"

# Extract intronic regions for RNA velocity analysis
gtf_introns="${build}/$(basename "$gtf_in").introns"
awk '$3 == "intron"' "$gtf_filtered" > "$gtf_introns"

# Create reference package for CellRanger
cellranger mkref --ref-version="$version" \
    --genome="$genome" --fasta="$fasta_modified" --genes="$gtf_filtered"

echo "Reference preparation complete. If you are running RNA velocity, use Velocyto after CellRanger."

