#!/bin/bash
#SBATCH --job-name=03_coverage_albertine_rift	# Job name
#SBATCH --partition=bi          # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jccooper@ku.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12              # Threads
#SBATCH --mem-per-cpu=4g            # Job memory request
#SBATCH --time=15-00:00:00	# time in days-hours:min:sec
#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/    # Directory
#SBATCH --output=03_coverage_albertine_rift.log     # Standard output and error log
#SBATCH --error=03_coverage_albertine_rift.err # error file

module load java bwa samtools

cd ~/scratch/albertine_rift/01_bam_files/

samtools depth -a Chamaetylas_poliophrys_FMNH346401_final.bam \
Chamaetylas_poliophrys_FMNH355643_final.bam \
Chamaetylas_poliophrys_FMNH385038_final.bam \
Chamaetylas_poliophrys_FMNH443853_final.bam \
Chamaetylas_poliophrys_FMNH443854_final.bam \
Chamaetylas_poliophrys_FMNH443856_final.bam \
Chamaetylas_poliophrys_FMNH450510_final.bam \
Chamaetylas_poliophrys_FMNH450511_final.bam \
Chamaetylas_poliophrys_FMNH450512_final.bam \
Cossypha_archeri_FMNH355604_final.bam \
Cossypha_archeri_FMNH358029_final.bam \
Cossypha_archeri_FMNH385023_final.bam \
Cossypha_archeri_FMNH438819_final.bam \
Cossypha_archeri_FMNH438820_final.bam \
Cossypha_archeri_FMNH443846_final.bam \
Cossypha_archeri_FMNH450492_final.bam \
Cossypha_semirufa_EB062_final.bam \
Melaenornis_chocolatinus_EB067_final.bam \
Phylloscopus_laetus_FMNH346427_final.bam \
Phylloscopus_laetus_FMNH355942_final.bam \
Phylloscopus_laetus_FMNH385174_final.bam \
Phylloscopus_laetus_FMNH443920_final.bam \
Phylloscopus_laetus_FMNH443922_final.bam \
Phylloscopus_laetus_FMNH450530_final.bam \
Phylloscopus_laetus_FMNH450531_final.bam \
Phylloscopus_umbrovirens_EB019_final.bam > ../lacustrine_coverage.txt

# break up the depth files into single column files for each individual (locations dropped)

cd ..

while read -r name1 number1; do
	number2=$((number1 + 2));
  cut lacustrine_coverage.txt -f $number2 > ${name1}_depth.txt;
done < dryo_popmap.txt
