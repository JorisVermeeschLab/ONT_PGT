#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --get-user-env="60L"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="15"
#SBATCH --ntasks="15"
#SBATCH --time="2-02:00:00"

# submit job
# sbatch --cluster=genius Clair3_refcall_hg38.slurm

SAMPLE="xx"
PROJECT_DIR="xx";
SAMPLE_DIR="$PROJECT_DIR/xx"
OUTPUT_DIR="$PROJECT_DIR/xx";
mkdir -p $OUTPUT_DIR

# Clair3 is a small variant caller for long-reads
conda activate clair3_single

BAM="$SAMPLE_DIR/$SAMPLE.trimmed.sorted.bam"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"
MODEL_PATH="/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/clair3_models"

# WhatsHap is a read-based phasing tool.  (default: do not phase indels)
# --print_ref_calls # Show reference calls (0/0) in vcf file, default: disable. 
#--qual=2  #mark variants with GQ>2 as PASS, GQ<=2 as LowQual. 2 is the tested default cutoff.

run_clair3.sh \
--bam_fn=$BAM \
--ref_fn=$REF \
--output=$OUTPUT_DIR \
--threads=15 \
--platform="ont" \
--model_path=$MODEL_PATH/xx \
--sample_name=$SAMPLE \
--remove_intermediate_dir \
--print_ref_calls \
--qual=2 