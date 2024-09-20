#!/bin/bash -l
#SBATCH --account="lp_joris_vermeesch"
#SBATCH --chdir="/user/leuven/339/vsc33900"
#SBATCH --error="xx/%x.e%A"
#SBATCH --output="xx/%x.o%A"
#SBATCH --export="NONE"
#SBATCH --get-user-env="60L"
#SBATCH --mail-type="BEGIN,END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="yan.zhao1@student.kuleuven.be"
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="12"
#SBATCH --ntasks="12"
#SBATCH --time="2-18:00:00"

# submit job
# sbatch --cluster=genius 7_Clair3_trio_newmodel.slurm


SAMPLE="" 
F=""
M=""

tmpdir="$VSC_SCRATCH/$SAMPLE"
mkdir -p $tmpdir
export TMPDIR=$tmpdir

PROJECT_DIR="";
OUTPUT_DIR="$PROJECT_DIR/7_Clair3_trio";
mkdir -p $OUTPUT_DIR

BAM="$PROJECT_DIR/5_mapped/$SAMPLE.trimmed.sorted.bam"

BAM_F="xx/$F.trimmed.sorted.bam"
BAM_M="xx/$M.trimmed.sorted.bam"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"
MODEL_PATH_clair3="/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/clair3_models"
MODEL_PATH_clair3_trio="/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/clair3_trio_models"

conda activate clair3
# use clair3 model r1041_e82_400bps_sup_v400 (Dorado)

/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/Clair3-Trio/run_clair3_trio.sh \
  --bam_fn_c=$BAM \
  --bam_fn_p1=$BAM_F \
  --bam_fn_p2=$BAM_M \
  --output=$OUTPUT_DIR \
  --ref_fn=$REF \
  --threads=30 \
  --model_path_clair3="$MODEL_PATH_clair3/r1041_e82_400bps_sup_v400" \
  --model_path_clair3_trio="$MODEL_PATH_clair3_trio/c3t_hg002_dna_r1041_e82_400bps_sup" \
  --sample_name_c=$SAMPLE \
  --sample_name_p1=$F \
  --sample_name_p2=$M