#!/bin/bash

# Filepaths
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../../" && pwd)
FIGURES="${PROJECT_DIR}/dist/figures"
BIN="${PROJECT_DIR}/bin"
DATA="${PROJECT_DIR}/data"
GENOMES="${DATA}/annotated_genomes"
BIGG_SEQS="${DATA}/bigg_seqs"
MODELS="${SCRIPT_DIR}/models"
RECONSTRUCTION="${SCRIPT_DIR}/reconstruction"
LOGS="${SCRIPT_DIR}/logs"
RESULTS="${SCRIPT_DIR}/results"

# Make sure needed directories exist
mkdir -p "$FIGURES"
mkdir -p "$GENOMES"
mkdir -p "$BIGG_SEQS"
mkdir -p "$MODELS"
mkdir -p "$RECONSTRUCTION"
mkdir -p "$LOGS"
mkdir -p "$RESULTS"

# Ignore generated files in git
echo "*" > "${BIGG_SEQS}/.gitignore"
echo "*" > "${MODELS}/.gitignore"
echo "*" > "${RECONSTRUCTION}/.gitignore"
echo "*" > "${LOGS}/.gitignore"
echo "*" > "${RESULTS}/.gitignore"

# Load .env variables (if .env exists)
ENV_FILE="${PROJECT_DIR}/.env"
if [ -f "$ENV_FILE" ]; then
    set -a
    source "$ENV_FILE"
    set +a
fi

# Set default python path if not set in .env
if [ -z "$PYTHON" ]; then
    PYTHON="${PROJECT_DIR}/.venv/bin/python"
fi

# Assume diamond on path if not set in .env
if [ -z "$DIAMOND" ]; then
    DIAMOND="diamond"
fi


# Download annotated genomes
ACCESSION_COL=4 # the column with accession numbers is the 4th column
START_ROW=2 # skip the header row
accessions=$(tail -n "+${START_ROW}" "${DATA}/model_info.csv" | cut -d ',' -f $ACCESSION_COL)
for accession in $accessions; do
    genome_path="${GENOMES}/${accession}.faa"
    if [ ! -f "$genome_path" ]; then
        echo "Downloading annotated genome - ${accession}..."
        "$PYTHON" "${BIN}/download_genome.py" "$accession" --out "$genome_path"
    fi
done


# Get BiGG universal models from CarveMe
echo "Getting BiGG universal models from CarveMe..."
get_carveme_universe() {
    # Arg 1: url for universal model
    # Arg 2: output path
    "$PYTHON" "${BIN}/get_carveme_universe.py" "$1" --out "$2"
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
    "$PYTHON" "${BIN}/make_bigg_template.py" "${DATA}/biomass_db.tsv" --out "$1" --biomass $2 --template-id $3 --threads 5
}

bigg_gramneg_template="${RECONSTRUCTION}/bigg_gramneg_template.sbml"
[ ! -f "$bigg_gramneg_template" ] && make_bigg_template "$bigg_gramneg_template" @gramneg gram_neg_template

bigg_grampos_template="${RECONSTRUCTION}/bigg_grampos_template.sbml"
[ ! -f "$bigg_grampos_template" ] && make_bigg_template "$bigg_grampos_template" @grampos gram_pos_template

bigg_generic_template="${RECONSTRUCTION}/bigg_generic_template.sbml"
[ ! -f "$bigg_generic_template" ] && make_bigg_template "$bigg_generic_template" @bacteria generic_bacteria_template


# Prepare the gene -> reaction map (can be the same for all datsets because DIAMOND db changes)
gene_rxn_map="${RECONSTRUCTION}/gene_rxn_map.json.gz"
if [ ! -f "$gene_rxn_map" ]; then
    echo "Preparing gene -> reaction map..."

    # Get list of BiGG model IDs that map will be based on
    bigg_ids=()
    {
        read -r header # skip header
        while IFS="," read -r _ _ bigg_id _
        do
            bigg_ids+=("$bigg_id")
        done
    } < "${DATA}/bigg_bacteria.csv"

    # Build map
    "$PYTHON" "${BIN}/make_gene_rxn_map.py" "${bigg_ids[@]}" --out "$gene_rxn_map"
fi


# Get all the protein sequences for species in BiGG
{
    read -r header # skip header
    while IFS="," read -r _ _ bigg_id refseq_accession
    do
        fasta="${BIGG_SEQS}/${bigg_id}.faa"
        if [ ! -f "$fasta" ]; then
            echo "Getting protein sequences for ${bigg_id} model..."
            "$PYTHON" "${BIN}/download_bigg_seqs.py" \
                "$bigg_id" \
                --refseq-accession "$refseq_accession" \
                --out "$fasta"
        fi
    done
} < "${DATA}/bigg_bacteria.csv"


# Job files
queue_reconstructor="${BIN}/queue_reconstructor.py"
array_job="${BIN}/slurm/job_array.sh"
reconstructor_job="${BIN}/slurm/reconstructor.sh"

job_ids=()


# For each annotation set, prep resource files and build models
{
    read -r header # read the header of the annotation set csv
    while IFS="," read -r set_id fraction seed namespace
    do
        # Directories for annotation set
        work_dir="${RECONSTRUCTION}/${set_id}"
        model_dir="${MODELS}/${set_id}"
        log_dir="${LOGS}/${set_id}"

        submitted_flag="${work_dir}/submitted.jobid"
        if [ ! -f "$submitted_flag" ]; then

            # Create directories
            [ ! -d "$work_dir" ] && mkdir "$work_dir"
            [ ! -d "$model_dir" ] && mkdir "$model_dir"
            [ ! -d "$log_dir" ] && mkdir "$log_dir"

            # Select species for this resource set
            selected_species="${work_dir}/selected_species.txt"
            if [ ! -f "$selected_species" ]; then
                echo "Selecting bacteria species for annotation set ${set_id}..."
                echo "Species fraction = ${fraction}"
                "$PYTHON" "${SCRIPT_DIR}/select_species.py" \
                    "${DATA}/bigg_bacteria.csv" \
                    --fraction "$fraction" \
                    --seed "$seed" \
                    --out "$selected_species"
            fi

            # Build DIAMOND database
            diamond_db="${work_dir}/protein_seqs.dmnd"
            if [ ! -f "$diamond_db" ]; then
                echo "Creating DIAMOND database for annotation set ${set_id}..."
                diamond_fasta="${work_dir}/protein_seqs.fasta"
                :> "$diamond_fasta"

                # Add sequences from selected species to fasta file
                while IFS= read -r bigg_id || [ -n "$bigg_id" ]; do
                    cat "${BIGG_SEQS}/${bigg_id}.faa" >> "$diamond_fasta"
                done < "$selected_species"

                # Convert fasta file to DIAMOND db
                "$DIAMOND" makedb --in "$diamond_fasta" --db "$diamond_db"
            fi

            # Submit jobs to build models
            echo "Submitting job to build models for annotation set ${set_id}..."
            "$PYTHON" "$queue_reconstructor" \
                "${DATA}/model_info.csv" \
                --genomes "$GENOMES" \
                --reconstruction "$work_dir" \
                --models "$model_dir" \
                --logs "$log_dir" \
                --gene-rxn-map "$gene_rxn_map" \
                --diamond-db "$diamond_db" \
                --diamond "$DIAMOND" \
                --media-db "${DATA}/mediadb.tsv" \
                --gapfill "minimal[bigg]" \
                --base-inputs "complete[bigg]" \
                --generic-universe "$bigg_generic_universe" \
                --gramneg-universe "$bigg_gramneg_universe" \
                --grampos-universe "$bigg_grampos_universe" \
                --generic-template "$bigg_generic_template" \
                --gramneg-template "$bigg_gramneg_template" \
                --grampos-template "$bigg_grampos_template" \
                --reconstructor-job "$reconstructor_job" \
                --array-job "$array_job" \
                --save-jobid "$submitted_flag"

            echo
        fi

        # Record the slurm job id
        job_ids+=($(cat "$submitted_flag"))

    done
} < "${SCRIPT_DIR}/datasets.csv"


# Prepare dependencies
deps=""
for job_id in "${job_ids[@]}"; do
    if squeue -j "$job_id" -h &>/dev/null; then
        deps+=":$job_id"
    fi
done

if [ -n "$deps" ]; then
    deps="afterok${deps}"
fi


# Queue job to process generated models
rxn_matrix="${RESULTS}/reaction_matrix.csv"
if [ ! -f "$rxn_matrix" ]; then
    echo "Queuing job to process generated models and create reaction matrix..."

    # Command arguments
    CMD=(
        "$PYTHON"
        "${BIN}/reaction_matrix.py"
        --model-info "${DATA}/model_info.csv"
        --set-info "${SCRIPT_DIR}/datasets.csv"
        --models "$MODELS"
        --out "$rxn_matrix"
        --include-bigg
    )

    # Submit job
    sbatch \
        --output="${LOGS}/reaction-matrix_%j.out" \
        --dependency="$deps" \
        --time=0:25:00 \
        --wrap="${CMD[*]}"
fi


# Queue job to run substrate utilization predictions
growth_data="${RESULTS}/growth_data.csv"
if [ ! -f "$growth_data" ]; then
    echo "Queuing job to run substrate utilization predictions..."

    # Command arguments
    CMD=(
        "$PYTHON"
        "${BIN}/substrate_utilization.py"
        --model-info "${DATA}/model_info.csv"
        --set-info "${SCRIPT_DIR}/datasets.csv"
        --models "$MODELS"
        --data "${DATA}/substrate_utilization_data.csv"
        --media-db "${DATA}/mediadb.tsv"
        --out "$growth_data"
    )

    # Submit job
    sbatch \
        --output="${LOGS}/substrate-utilization_%j.out" \
        --dependency="$deps" \
        --time=0:30:00 \
        --wrap="${CMD[*]}"
fi


# Compile metrics for resource sets
dataset_metrics="${RESULTS}/dataset_metrics.csv"
if [ ! -f "$dataset_metrics" ]; then
    echo "Calculating resource set coverage and diversity..."
    "$PYTHON" "${SCRIPT_DIR}/calculate_set_metrics.py" \
        "${SCRIPT_DIR}/datasets.csv" \
        --out  "$dataset_metrics" \
        --reconstruction-dir "$RECONSTRUCTION" \
        --bigg-info "${DATA}/bigg_bacteria.csv"
fi
