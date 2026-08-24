#!/bin/bash

# Filepaths
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../../" && pwd)
FIGURES="${PROJECT_DIR}/dist/figures"
VENV="${PROJECT_DIR}/.venv"
BIN="${PROJECT_DIR}/bin"
DATA="${PROJECT_DIR}/data"
GENOMES="${DATA}/annotated_genomes"
MODELS="${SCRIPT_DIR}/models"
RECONSTRUCTION="${SCRIPT_DIR}/reconstruction"
LOGS="${SCRIPT_DIR}/logs"

# Make sure needed directories exist
mkdir -p "$FIGURES"
mkdir -p "$GENOMES"
mkdir -p "$MODELS"
mkdir -p "$RECONSTRUCTION"
mkdir -p "$LOGS"

# Ignore generated files in git
echo "*" > "${MODELS}/.gitignore"
echo "*" > "${RECONSTRUCTION}/.gitignore"
echo "*" > "${LOGS}/.gitignore"

# Load .env variables (if .env exists)
ENV_FILE="${PROJECT_DIR}/.env"
if [ -f "$env" ]; then
    set -a
    source "$ENV_FILE"
    set +a
fi

# Activate python venv
source "${VENV}/bin/activate"


# Download annotated genomes
ACCESSION_COL=4 # the column with accession numbers is the 4th column
START_ROW=2 # skip the header row
accessions=$(tail -n "+${START_ROW}" "${DATA}/model_info.csv" | cut -d ',' -f $ACCESSION_COL)
for accession in $accessions; do
    genome_path="${GENOMES}/${accession}.faa"
    if [ ! -f "$genome_path" ]; then
        echo "Downloading annotated genome - ${accession}..."
        python "${BIN}/download_genome.py" "$accession" --out "$genome_path"
    fi
done


# Get ModelSEED universal models from Reconstructor
seed_gramneg_universe="${RECONSTRUCTION}/seed_gramneg_universe.sbml.gz"
seed_grampos_universe="${RECONSTRUCTION}/seed_grampos_universe.sbml.gz"
if [[ ! -f "$seed_gramneg_universe" || ! -f "$seed_grampos_universe" ]]; then
    echo "Getting ModelSEED universal models from Reconstructor..."
    python "${BIN}/get_reconstructor_universe.py" \
        --gramneg "$seed_gramneg_universe" \
        --grampos "$seed_grampos_universe"
fi


# Create template models for ModelSEED namespace
seed_gramneg_template="${RECONSTRUCTION}/seed_gramneg_template.sbml"
seed_grampos_template="${RECONSTRUCTION}/seed_grampos_template.sbml"
if [[ ! -f "$seed_gramneg_template" || ! -f "$seed_grampos_template" ]]; then
    echo "Preparing template models for ModelSEED namespace..."
    python "${BIN}/get_reconstructor_templates.py" \
        "$seed_universe" \
        --gramneg "$seed_gramneg_template" \
        --grampos "$seed_grampos_template"
fi


# Get KEGG protein sequences
kegg_seqs="${RECONSTRUCTION}/kegg_proteins.dmnd"
if [ ! -f "$kegg_seqs" ]; then
    echo "Getting KEGG DIAMOND database and fixing gene IDs for CarveMe..."
    module load diamond/2.1.6

    # Download DB
    db_tmp_path="${RECONSTRUCTION}/screened_kegg_prokaryotes_pep_db.dmnd"
    db_url="https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
    curl -L "$db_url" -o "$db_tmp_path"

    # Process gene IDs in DB
    diamond getseq \
        -d "$db_tmp_path" \
        | sed -E 's/^>([a-z]+):(\w+)/>\1.\2/' \
        | diamond makedb -d "$kegg_seqs" --quiet

    # Remove temporary DB
    rm "$db_tmp_path"

    module unload diamond
fi


# Get BiGG universal models from CarveMe
echo "Getting BiGG universal models from CarveMe..."
get_carveme_universe() {
    # Arg 1: url for universal model
    # Arg 2: output path
    python "${BIN}/get_carveme_universe.py" "$1" --out "$2"
}

bigg_generic_universe="${RECONSTRUCTION}/bigg_universe_bacteria.xml.gz"
generic_url="https://github.com/cdanielmachado/carveme/raw/refs/heads/master/carveme/data/generated/universe_bacteria.xml.gz"
[ ! -f "$bigg_generic_universe" ] && get_carveme_universe "$generic_url" "$bigg_generic_universe"

bigg_gramneg_universe="${RECONSTRUCTION}/bigg_universe_gramneg.xml.gz"
gramneg_url="https://github.com/cdanielmachado/carveme/raw/refs/heads/master/carveme/data/generated/universe_gramneg.xml.gz"
[ ! -f "$bigg_gramneg_universe" ] && get_carveme_universe "$gramneg_url" "$bigg_gramneg_universe"

bigg_grampos_universe="${RECONSTRUCTION}/bigg_universe_grampos.xml.gz"
grampos_url="https://github.com/cdanielmachado/carveme/raw/refs/heads/master/carveme/data/generated/universe_grampos.xml.gz"
[ ! -f "$bigg_grampos_universe" ] && get_carveme_universe "$grampos_url" "$bigg_grampos_universe"


# Create template models for BiGG namespace
echo "Preparing template models for BiGG namespace..."
make_bigg_template () {
    python "${BIN}/make_bigg_template.py" "${DATA}/biomass_db.tsv" --out "$1" --biomass $2 --template-id $3 --threads 5
}

bigg_gramneg_template="${RECONSTRUCTION}/bigg_gramneg_template.sbml"
[ ! -f "$bigg_gramneg_template" ] && make_bigg_template "$bigg_gramneg_template" @gramneg gram_neg_template

bigg_grampos_template="${RECONSTRUCTION}/bigg_grampos_template.sbml"
[ ! -f "$bigg_grampos_template" ] && make_bigg_template "$bigg_grampos_template" @grampos gram_pos_template

bigg_generic_template="${RECONSTRUCTION}/bigg_generic_template.sbml"
[ ! -f "$bigg_generic_template" ] && make_bigg_template "$bigg_generic_template" @bacteria generic_bacteria_template


# Get BiGG protein sequences from CarveMe
bigg_seqs="${RECONSTRUCTION}/bigg_proteins.dmnd"
if [ ! -f "$bigg_seqs" ]; then
    echo "Downloading BiGG DIAMOND database..."
    bigg_seqs_url="https://github.com/cdanielmachado/carveme/raw/refs/heads/master/carveme/data/generated/bigg_proteins.dmnd"
    curl -L "$bigg_seqs_url" -o "$bigg_seqs"
fi


# Prepare files mapping genes to reactions
bigg_gprs="${RECONSTRUCTION}/bigg_gprs.csv.gz"
bigg_gene_rxn_map="${RECONSTRUCTION}/bigg_gene_rxn_map.json.gz"
kegg_gprs="${RECONSTRUCTION}/kegg_seed_gprs.csv"
kegg_gene_rxn_map="${RECONSTRUCTION}/kegg_gene_rxn_map.json.gz"
if [[ ! ( -f "$bigg_gprs" && -f "$bigg_gene_rxn_map" && -f "$kegg_gprs" && -f "$kegg_gene_rxn_map" ) ]]; then
    echo "Preparing GPR files..."
    python "${BIN}/prep_GPR_files.py" \
        --bigg-gprs "$bigg_gprs" \
        --bigg-map "$bigg_gene_rxn_map" \
        --kegg-gprs "$kegg_gprs" \
        --kegg-map "$kegg_gene_rxn_map"
fi


# Create an empty annotation file (because we don't need gene annotations for this analysis)
annotation_file="${RECONSTRUCTION}/gene_annotations.tsv"
if [ ! -f "$annotation_file" ]; then
    echo "Creating empty annotation file..."
    printf $'model\tgene\tannotation\n' > "$annotation_file"
    echo
fi


# Job files
queue_reconstructor="${BIN}/queue_reconstructor.py"
queue_carveme="${BIN}/queue_carveme.py"
array_job="${BIN}/slurm/job_array.sh"
reconstructor_job="${BIN}/slurm/reconstructor.sh"
carveme_job="${BIN}/slurm/carveme.sh"


# Queue jobs to build models with Reconstructor + KEGG/ModelSEED data
rks_flag="${RECONSTRUCTION}/reconstructor_kegg_seed/submitted.jobid"
if [ ! -f "$rks_flag" ]; then
    echo "Queuing jobs to build models with Reconstructor + KEGG/ModelSEED data..."
    python "$queue_reconstructor" \
        "${DATA}/model_info.csv" \
        --genomes "$GENOMES" \
        --reconstruction "${RECONSTRUCTION}/reconstructor_kegg_seed" \
        --models "${MODELS}/reconstructor_kegg_seed" \
        --logs "${LOGS}/reconstructor_kegg_seed" \
        --venv "${VENV}/bin/activate" \
        --gene-rxn-map "$kegg_gene_rxn_map" \
        --diamond-db "$kegg_seqs" \
        --media-db "${DATA}/mediadb.tsv" \
        --gapfill minimal[seed] \
        --base-inputs complete[seed] \
        --generic-universe "$seed_gramneg_universe" \
        --gramneg-universe "$seed_gramneg_universe" \
        --grampos-universe "$seed_grampos_universe" \
        --generic-template "$seed_gramneg_template" \
        --gramneg-template "$seed_gramneg_template" \
        --grampos-template "$seed_grampos_template" \
        --reconstructor-job "$reconstructor_job" \
        --array-job "$array_job" \
        --save-jobid "$rks_flag"
fi

rks_jobid=$(cat "$rks_flag")


# Queue jobs to build models with CarveMe + KEGG/ModelSEED data
cks_flag="${RECONSTRUCTION}/carveme_kegg_seed/submitted.jobid"
if [ ! -f "$cks_flag" ]; then
    echo "Queuing jobs to build models with CarveMe + KEGG/ModelSEED data..."
    python "$queue_carveme" \
        "${DATA}/model_info.csv" \
        --genomes "$GENOMES" \
        --reconstruction "${RECONSTRUCTION}/carveme_kegg_seed" \
        --models "${MODELS}/carveme_kegg_seed" \
        --logs "${LOGS}/carveme_kegg_seed" \
        --venv "${VENV}/bin/activate" \
        --gprs "$kegg_gprs" \
        --annotations "$annotation_file" \
        --diamond-db "$kegg_seqs" \
        --media-db "${DATA}/mediadb.tsv" \
        --gapfill minimal[seed] \
        --generic-universe "$seed_gramneg_universe" \
        --gramneg-universe "$seed_gramneg_universe" \
        --grampos-universe "$seed_grampos_universe" \
        --carveme-job "$carveme_job" \
        --array-job "$array_job" \
        --save-jobid "$cks_flag"
fi

cks_jobid=$(cat "$cks_flag")


# Queue jobs to build models with Reconstructor + BiGG data
rb_flag="${RECONSTRUCTION}/reconstructor_bigg/submitted.jobid"
if [ ! -f "$rb_flag" ]; then
    echo "Queuing jobs to build models with Reconstructor + BiGG data..."
    python "$queue_reconstructor" \
        "${DATA}/model_info.csv" \
        --genomes "$GENOMES" \
        --reconstruction "${RECONSTRUCTION}/reconstructor_bigg" \
        --models "${MODELS}/reconstructor_bigg" \
        --logs "${LOGS}/reconstructor_bigg" \
        --venv "${VENV}/bin/activate" \
        --gene-rxn-map "$bigg_gene_rxn_map" \
        --diamond-db "$bigg_seqs" \
        --media-db "${DATA}/mediadb.tsv" \
        --gapfill minimal[bigg] \
        --base-inputs complete[bigg] \
        --generic-universe "$bigg_gramneg_universe" \
        --gramneg-universe "$bigg_gramneg_universe" \
        --grampos-universe "$bigg_grampos_universe" \
        --generic-template "$bigg_gramneg_template" \
        --gramneg-template "$bigg_gramneg_template" \
        --grampos-template "$bigg_grampos_template" \
        --reconstructor-job "$reconstructor_job" \
        --array-job "$array_job" \
        --save-jobid "$rb_flag"
fi

rb_jobid=$(cat "$rb_flag")


# Queue jobs to build models with CarveMe + BiGG data
cb_flag="${RECONSTRUCTION}/carveme_bigg/submitted.jobid"
if [ ! -f "$cb_flag" ]; then
    echo "Queuing jobs to build models with CarveMe + BiGG data..."
    python "$queue_carveme" \
        "${DATA}/model_info.csv" \
        --genomes "$GENOMES" \
        --reconstruction "${RECONSTRUCTION}/carveme_bigg" \
        --models "${MODELS}/carveme_bigg" \
        --logs "${LOGS}/carveme_bigg" \
        --venv "${VENV}/bin/activate" \
        --gprs "$bigg_gprs" \
        --annotations "$annotation_file" \
        --diamond-db "$bigg_seqs" \
        --media-db "${DATA}/mediadb.tsv" \
        --gapfill minimal[bigg] \
        --generic-universe "$bigg_gramneg_universe" \
        --gramneg-universe "$bigg_gramneg_universe" \
        --grampos-universe "$bigg_grampos_universe" \
        --carveme-job "$carveme_job" \
        --array-job "$array_job" \
        --save-jobid "$cb_flag"
fi

cb_jobid=$(cat "$cb_flag")
