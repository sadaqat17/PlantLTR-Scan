#!/bin/bash
#SBATCH --job-name=Alphonso
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
#SBATCH -p ecobio,genouest

# -----------------------------
# Environment Setup
# -----------------------------
source /home/$USER/miniconda3/etc/profile.d/conda.sh
conda activate ltr_finder-1.07
# -----------------------------
# Help Function
# -----------------------------
function usage() {
    echo "Usage: $0 -i <input_directory> -s <species_name> [-t <threads>] [-b <bootstraps>]"
    exit 1
}

# -----------------------------
# Argument Parsing
# -----------------------------
THREADS=1  # Default
BOOTSTRAPS=""
SPECIES=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i) INPUT_DIR="$2"; shift 2 ;;
		-s) SPECIES="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -b) BOOTSTRAPS="$2"; shift 2 ;;
        *) usage ;;
    esac
done

# -----------------------------
# Input Validation
# -----------------------------
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory not provided."
    usage
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' does not exist."
    exit 1
fi

if [ -z "$SPECIES" ]; then
    echo "Error: Species name not provided."
    usage
fi

echo "Input directory: $INPUT_DIR"
echo "Species: $SPECIES"


# -----------------------------
# Input File Checks
# -----------------------------
FASTA_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) | head -n 1)
GFF_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f \( -iname "*.gff" -o -iname "*.gff3" -o -iname "*.gtf" \) | head -n 1)
PEP_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f -iname "*.pep" | head -n 1)

if [ -z "$FASTA_FILE" ] || [ -z "$GFF_FILE" ] || [ -z "$PEP_FILE" ]; then
    echo "Error: Required input files are missing."
    echo "FASTA: $FASTA_FILE"
    echo "GFF: $GFF_FILE"
    echo "PEP: $PEP_FILE"
    exit 1
fi

echo "All required files found:"
echo "- FASTA: $FASTA_FILE"
echo "- GFF: $GFF_FILE"
echo "- PEP: $PEP_FILE"

# -----------------------------
# LTR_FINDER
# -----------------------------
BASE_NAME=$(basename "$FASTA_FILE" | sed 's/\.[^.]*$//')
LTR_TMP="tmp/ltr_finder"
LTR_OUT="$LTR_TMP/${BASE_NAME}.finder.scn"
LTR_SUMMARY="Results/ltr_finder_summary.tsv"
LTR_BAR="Results/ltr_finder_barplot.svg"

mkdir -p "$LTR_TMP" Results/

echo "Running LTR_FINDER..."
python ./bin/run_ltrfinder_parallel.py -i "$FASTA_FILE" -o "$LTR_OUT" -t "$THREADS"
python ./bin/parse_ltr_finder.py -i "$LTR_OUT" -o "$LTR_SUMMARY"
python ./bin/ltr_finder_barplot.py -i "$LTR_SUMMARY" -o "$LTR_BAR"
conda deactivate
conda activate ltr_retriever-2.9.5
# -----------------------------
# LTR_RETRIEVER
# -----------------------------
LTR_RET_TMP="tmp/ltr_retriever"
LTR_RET_SUMMARY="Results/LTR_retriever_summary.tsv"
LTR_RET_BAR="Results/ltr_retriever_barplot.svg"
LTR_BED="$LTR_RET_TMP/ltrs.bed"
LTR_FASTA="$LTR_RET_TMP/final_ltrs.fasta"

mkdir -p "$LTR_RET_TMP"
cp "$LTR_OUT" "$LTR_RET_TMP/"

echo "Running LTR_retriever..."
cd "$LTR_RET_TMP" || exit 1
LTR_retriever -genome "../../$FASTA_FILE" -infinder "$(basename "$LTR_OUT")" -threads "$THREADS" -noanno || { echo "❌ LTR_retriever failed."; exit 1; }

rm *.nmtf.pass.list
PASS_LIST=$(find . -name "*.pass.list" | head -n 1)
[ -f "$PASS_LIST" ] || { echo ".pass.list not found!"; exit 1; }
cp "$PASS_LIST" "../../$LTR_RET_SUMMARY"
cd ../../

python ./bin/stacked_te_plot.py -i "$LTR_RET_SUMMARY" -o "$LTR_RET_BAR"
python ./bin/parse_ltr_passlist_to_bed.py --input "$LTR_RET_SUMMARY" --fasta "$FASTA_FILE" --bed "$LTR_BED" --output_fasta "$LTR_FASTA"

# -----------------------------
# TEsorter
# -----------------------------
TE_SORT_TMP="tmp/TEsorter"
TE_IN="$TE_SORT_TMP/final_ltrs.fasta"
TE_RAW="$TE_SORT_TMP/final_ltrs.fasta.rexdb-plant.cls.lib"
TE_CLEAN="$TE_SORT_TMP/final_ltrs.fasta.rexdb-plant.fasta"
TE_FINAL="Results/TEsorter_final_LTRs.fasta"
COPIA_FASTA="Results/copia.fasta"
GYPSY_FASTA="Results/gypsy.fasta"

mkdir -p "$TE_SORT_TMP"
cp "$LTR_FASTA" "$TE_IN"

echo "Running TEsorter in $TE_SORT_TMP..."
cd "$TE_SORT_TMP" || exit 1
TEsorter final_ltrs.fasta -p "$THREADS" -cov 10 -eval 1e-2 -db rexdb-plant
cd ../../

[ -f "$TE_RAW" ] || { echo "TEsorter output not found: $TE_RAW"; exit 1; }
awk '/^>/ {sub(/ .*/, "", $0)} {print}' "$TE_RAW" > "$TE_CLEAN"
cp "$TE_CLEAN" "$TE_FINAL"
awk '/^>/ {f=($0~/#LTR\/Copia/)} f' "$TE_CLEAN" > "$COPIA_FASTA"
awk '/^>/ {f=($0~/#LTR\/Gypsy/)} f' "$TE_CLEAN" > "$GYPSY_FASTA"
python ./bin/extract_chr_column.py -i ./tmp/TEsorter/final_ltrs.fasta.rexdb-plant.dom.gff3 -o ./tmp/TEsorter/final_ltrs.chrid.gff3

# -----------------------------
# Phylogenetics
# -----------------------------
ALIGN_DIR="tmp/alignment"
COPIA_ALN="$ALIGN_DIR/copia.aln.fasta"
GYPSY_ALN="$ALIGN_DIR/gypsy.aln.fasta"

mkdir -p "$ALIGN_DIR"

echo "Aligning sequences..."
mafft --thread "$THREADS" --auto --reorder "$COPIA_FASTA" > "$COPIA_ALN"
mafft --thread "$THREADS" --auto --reorder "$GYPSY_FASTA" > "$GYPSY_ALN"

# -----------------------------
# Phylogenetic Inference
# -----------------------------
PHYLO_DIR="tmp/phylogenetics"
COPIA_TREE_TMP="$PHYLO_DIR/copia.aln.fasta.treefile"
GYPSY_TREE_TMP="$PHYLO_DIR/gypsy.aln.fasta.treefile"
COPIA_TREE_FINAL="Results/copia.treefile"
GYPSY_TREE_FINAL="Results/gypsy.treefile"

mkdir -p "$PHYLO_DIR"

echo "Building trees with IQ-TREE..."
cp "$COPIA_ALN" "$PHYLO_DIR/"
cp "$GYPSY_ALN" "$PHYLO_DIR/"

cd "$PHYLO_DIR" || exit 1

if [ -n "$BOOTSTRAPS" ]; then
  iqtree -s "$(basename "$COPIA_ALN")" -m TEST -nt AUTO -keep-ident -B "$BOOTSTRAPS"
  iqtree -s "$(basename "$GYPSY_ALN")" -m TEST -nt AUTO -keep-ident -B "$BOOTSTRAPS"
else
  iqtree -s "$(basename "$COPIA_ALN")" -m TEST -nt AUTO -keep-ident
  iqtree -s "$(basename "$GYPSY_ALN")" -m TEST -nt AUTO -keep-ident
fi

# Move back to main directory
cd ../../

cp "$COPIA_TREE_TMP" "$COPIA_TREE_FINAL"
cp "$GYPSY_TREE_TMP" "$GYPSY_TREE_FINAL"

# -----------------------------
# GFF Comparison
# -----------------------------
GFF_COMPARE_DIR="tmp/gffcompare"
GFF_INPUT_REF="$GFF_FILE" 
GFF_PRED="tmp/TEsorter/final_ltrs.chrid.gff3"
GFF_OUTPUT_PREFIX="Results/gffcompare_out"

mkdir -p "$GFF_COMPARE_DIR"
cp "$GFF_PRED" "$GFF_COMPARE_DIR/"
cp "$GFF_INPUT_REF" "$GFF_COMPARE_DIR/"

cd "$GFF_COMPARE_DIR" || exit 1

echo "Running gffcompare..."
gffcompare -r "$(basename "$GFF_INPUT_REF")" -o "$(basename "$GFF_OUTPUT_PREFIX")" "$(basename "$GFF_PRED")" || { echo "❌ gffcompare failed."; exit 1; }

cd ../../
python ./bin/filter_tracking.py -i ./tmp/gffcompare/*.tracking -o ./Results/Gffcompare_summary.txt
python ./bin/fetch_proteins.py -i ./Results/Gffcompare_summary.txt -p "$PEP_FILE" -a "$GFF_FILE" -o ./Results/Gene_inserted_LTR_proteins.fasta
echo "Protein Seuences of Gene Inserted LTR-RTs retrived successfully!"
conda deactivate

# -----------------------------
# Pannzer2 Gene Ontology
# -----------------------------
conda activate pannzer2
mkdir ./tmp/Gene_ontology
cp -r ./Results/Gene_inserted_LTR_proteins.fasta ./tmp/Gene_ontology
cd ./tmp/Gene_ontology
echo "Running gene ontology analysis..."
python ../../bin/SANSPANZ.3/runsanspanz.py -R -o ",DE.out,GO.out,anno.out" -s "$SPECIES" < Gene_inserted_LTR_proteins.fasta
python ../../bin/ontology_summary.py -i GO.out -o ../../Results/GO
echo "Pipeline completed successfully!"
conda deactivate