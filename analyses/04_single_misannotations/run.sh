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


# Prepare the DIAMOND database
diamond_db="${RECONSTRUCTION}/protein_seqs.dmnd"
if [ ! -f "$diamond_db" ]; then
    echo "Preparing DIAMOND database..."
    diamond_fasta="${RECONSTRUCTION}/protein_seqs.fasta"
    :> "$diamond_fasta"

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
            cat "$fasta" >> "$diamond_fasta"
        done
    } < "${DATA}/bigg_bacteria.csv"

    # Convert fasta file to DIAMOND db
    "$DIAMOND" makedb --in "$diamond_fasta" --db "$diamond_db"
fi


# Prepare the base gene -> reaction map
base_gene_rxn_map="${RECONSTRUCTION}/base_gene_rxn_map.json.gz"
if [ ! -f "$base_gene_rxn_map" ]; then
    echo "Preparing base gene -> reaction map..."

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
    "$PYTHON" "${BIN}/make_gene_rxn_map.py" "${bigg_ids[@]}" --out "$base_gene_rxn_map"
fi


# Make an annotation pool that has all unique annotations from the base gene -> rxn map
annotation_pool="${RECONSTRUCTION}/unique_annotations.json"
if [ ! -f "$annotation_pool" ]; then
    "$PYTHON" "${SCRIPT_DIR}/get_unique_annotations.py" "$base_gene_rxn_map" --out "$annotation_pool"
fi


# Determine gram stain of each organism (for picking resource files) and add random seeds
organism_info="${RECONSTRUCTION}/organism_info.csv"
if [ ! -f "$organism_info" ]; then
    echo "Determining gram stain of organisms for picking resource files..."
    "$PYTHON" "${SCRIPT_DIR}/prep_organism_info.py" \
        "${DATA}/model_info.csv" \
        --extra "${SCRIPT_DIR}/random_seeds.csv" \
        --out "$organism_info"
fi


# For each organism, build models and run substrate utilization predictions
{
    read -r header # skip the header row
    while IFS="," read -r organism_id accession bigg_model gram_stain seed
    do
        # Directories for annotation set
        work_dir="${RECONSTRUCTION}/${organism_id}"
        model_dir="${MODELS}/${organism_id}"
        log_dir="${LOGS}/${organism_id}"

        # Create directories
        [ ! -d "$work_dir" ] && mkdir "$work_dir"
        [ ! -d "$model_dir" ] && mkdir "$model_dir"
        [ ! -d "$log_dir" ] && mkdir "$log_dir"


        RUN_REMAINING=0     # Flag to signal at what point all remaining steps need to be run (because an earlier step was run)


        # Run blast on input genome for organism
        blast_flag="${work_dir}/blast.submitted.jobid"
        blast_output="${work_dir}/${accession}.blast.out"
        blast_jobid=""
        if [ ! -f "$blast_flag" ]; then
            echo "Submitting BLAST job for ${organism_id}..."
            blast_cmd=(
                "$DIAMOND"
                blastp
                --db "$diamond_db"
                --query "${GENOMES}/${accession}.faa"
                --out "$blast_output"
                --more-sensitive
                --max-target-seqs 1
            )
            blast_jobid=$(
                sbatch \
                    --parsable \
                    --job-name=blast \
                    --output="${log_dir}/%x-%j.out" \
                    --time=5:00 \
                    --mem=512M \
                    --ntasks=1 \
                    --cpus-per-task=4 \
                    --wrap="${blast_cmd[*]}"
            )

            echo $blast_jobid > "$blast_flag"
            RUN_REMAINING=1
        elif [ ! -f "$blast_output" ]; then     # job submitted but not done yet
            blast_jobid=$(cat "$blast_flag")
        fi

        # Dependency on blast job for next job
        dep=""
        if [ -n "$blast_jobid" ]; then
            dep=afterok:${blast_jobid}
        fi

        # Define the annotation sets for the organism
        sets_flag="${work_dir}/define_sets.submitted.jobid"
        annotation_sets="${work_dir}/annotation_sets.csv"
        sets_jobid=""
        if [[ ! -f "$sets_flag" || $RUN_REMAINING -eq 1 ]]; then
            echo "Submitting job to define annotation sets for ${organism_id}..."
            cmd=(
                "$PYTHON"
                "${SCRIPT_DIR}/define_annotation_sets.py"
                "$blast_output"
                --out "$annotation_sets"
                --annotation-pool "$annotation_pool"
                -n 3
                --seed $seed
            )
            sets_jobid=$(
                sbatch \
                    --parsable \
                    --job-name=define_sets \
                    --output="${log_dir}/%x-%j.out" \
                    --dependency=${dep} \
                    --time=5:00 \
                    --mem=512M \
                    --wrap="${cmd[*]}"
            )

            echo $sets_jobid > "$sets_flag"
            RUN_REMAINING=1
        elif [ ! -f "$annotation_sets" ]; then      # job submitted but not done yet
            sets_jobid=$(cat "$sets_flag")
        fi

        # Dependency on define_sets job for next job
        dep=""
        if [ -n "$sets_jobid" ]; then
            dep=afterok:${sets_jobid}
        fi

        # Build the models
        build_flag="${work_dir}/build.submitted.jobid"
        finished_flag="${work_dir}/build.completed.jobid"
        build_jobid=""
        if [[ ! -f "$build_flag" || $RUN_REMAINING -eq 1 ]]; then
            echo "Submitting job to build models for ${organism_id}..."

            # Pick resource files based on gram stain
            if [ "$gram_stain" = "negative" ]; then
                universe="$bigg_gramneg_universe"
                template="$bigg_gramneg_template"
            elif [ "$gram_stain" = "positive" ]; then
                universe="$bigg_grampos_universe"
                template="$bigg_grampos_template"
            else
                universe="$bigg_generic_universe"
                template="$bigg_generic_template"
            fi

            cmd=(
                "$PYTHON"
                "${SCRIPT_DIR}/build_model_set.py"
                "$annotation_sets"
                --out "$model_dir"
                --blast "$blast_output"
                --universe "$universe"
                --gene-rxn-map "$base_gene_rxn_map"
                --template "$template"
                --mediadb "${DATA}/mediadb.tsv"
                --processes 40
            )
            if [ $RUN_REMAINING -eq 1 ]; then
                cmd+=(--rebuild)
            fi

            build_jobid=$(
                sbatch \
                    --parsable \
                    --job-name=build_${organism_id} \
                    --output="${log_dir}/%x-%j.out" \
                    --dependency=${dep} \
                    --time=4:00:00 \
                    --mem-per-cpu=4G \
                    --ntasks=1 \
                    --cpus-per-task=40 \
                    --wrap="${cmd[*]}"
            )

            # Job to create flag file after build job finishes successfully
            if [ -f "$finished_flag" ]; then
                rm "$finished_flag"
            fi
            sbatch \
                --quiet \
                --job-name=post-build_${organism_id} \
                --output="${log_dir}/%x-%j.out" \
                --dependency=afterok:${build_jobid} \
                --time=1:00 \
                --wrap="echo ${build_jobid} > ${finished_flag}"

            echo $build_jobid > "$build_flag"
            RUN_REMAINING=1
        elif [ ! -f "$finished_flag" ]; then
            build_jobid=$(cat "$build_flag")
        fi

        # Dependency on build job for next jobs
        dep=""
        if [ -n "$build_jobid" ]; then
            dep=afterok:${build_jobid}
        fi

        # Process reaction content
        reactions_flag="${work_dir}/reactions.submitted.jobid"
        reactions_out="${MODELS}/reactions_${organism_id}.csv"
        if [[ ! -f "$reactions_flag" || $RUN_REMAINING -eq 1 ]]; then
            echo "Submitting job to process reactions for ${organism_id}..."
            cmd=(
                "$PYTHON"
                "${SCRIPT_DIR}/process_reactions.py"
                $organism_id
                --out "$reactions_out"
                --set-info "$annotation_sets"
                --models "$model_dir"
                --include-bigg $bigg_model
            )
            reactions_jobid=$(
                sbatch \
                    --parsable \
                    --job-name=reactions \
                    --output="${log_dir}/%x-%j.out" \
                    --dependency=${dep} \
                    --time=10:00 \
                    --mem=4G \
                    --wrap="${cmd[*]}"
            )
        fi

        # Substrate utilization predictions
        predictions_flag="${work_dir}/predictions.submitted.jobid"
        predictions_out="${MODELS}/growth_data_${organism_id}.csv"
        if [[ ! -f "$predictions_flag" || $RUN_REMAINING -eq 1 ]]; then
            echo "Submitting job to run substrate utilization predictions for ${organism_id}..."
            cmd=(
                "$PYTHON"
                "${SCRIPT_DIR}/substrate_utilization.py"
                $organism_id
                --out "$predictions_out"
                --set-info "$annotation_sets"
                --models "$model_dir"
                --exp-data "${DATA}/substrate_utilization_data.csv"
                --mediadb "${DATA}/mediadb.tsv"
                --processes 40
            )
            predictions_jobid=$(
                sbatch \
                    --parsable \
                    --job-name=sub_util \
                    --output="${log_dir}/%x-%j.out" \
                    --dependency=${dep} \
                    --time=30:00 \
                    --mem-per-cpu=1G \
                    --ntasks=1 \
                    --cpus-per-task=40 \
                    --wrap="${cmd[*]}"
            )
            echo $predictions_jobid > "$predictions_flag"
        fi
    done
} < "$organism_info"
