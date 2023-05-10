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
#SBATCH --nodes="2"
#SBATCH --ntasks-per-node="15"
#SBATCH --ntasks="30"
#SBATCH --time="3-00:00:00"

# submit job
# sbatch --cluster=genius Clair3_trio_hg38.slurm


tmpdir="$VSC_SCRATCH/xx"
mkdir -p $tmpdir
export TMPDIR=$tmpdir

PROJECT_DIR="xx";
OUTPUT_DIR="$PROJECT_DIR/xx";
mkdir -p $OUTPUT_DIR

# Clair3 is a small variant caller for long-reads
conda activate clair3_single

BAM="$PROJECT_DIR/xx"

BAM_GM24149_HG003="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24149/6_GRCh38/6_mapped/Father_GM24149_HG003.rg.changed.bam"
BAM_GM24143_HG004="/lustre1/project/stg_00019/research/yan/nanopore_data/00_04_resource_parental_nanopore_data/GM24143/6_GRCh38/1_mapped/Mother_GM24143_HG004.rg.changed.bam"
REF="/lustre1/project/stg_00019/genome/homo_sapiens/hg38/genome.fa"
MODEL_PATH_clair3="/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/clair3_models"
MODEL_PATH_clair3_trio="/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/clair3_trio_models"

#--qual=2  #mark variants with GQ>2 as PASS, GQ<=2 as LowQual. 2 is the tested default cutoff.

conda activate clair3

/vsc-hard-mounts/leuven-data/339/vsc33900/A00_software/Clair3-Trio/run_clair3_trio.sh \
  --bam_fn_c=$BAM \
  --bam_fn_p1=$BAM_GM24149_HG003 \
  --bam_fn_p2=$BAM_GM24143_HG004 \
  --output=$OUTPUT_DIR \
  --ref_fn=$REF \
  --threads=30 \
  --model_path_clair3="$MODEL_PATH_clair3/xx" \
  --model_path_clair3_trio="$MODEL_PATH_clair3_trio/c3t_hg002_r941_prom_sup_g5014" \
  --sample_name_c=HG002 \
  --sample_name_p1=HG003 \
  --sample_name_p2=HG004 \
  --remove_intermediate_dir \
  --qual=2  