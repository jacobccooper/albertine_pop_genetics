#!/bin/bash
#SBATCH --job-name=04_flagstat_albertine_rift	# Job name
#SBATCH --partition=bi          # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jccooper@ku.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4              # Threads
#SBATCH --mem-per-cpu=8g            # Job memory request
#SBATCH --time=15-00:00:00	# time in days-hours:min:sec
#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/    # Directory
#SBATCH --output=04_flagstat_albertine_rift.log     # Standard output and error log
#SBATCH --error=04_flagstat_albertine_rift.err # error file

module load samtools

workdir=~/scratch/albertine_rift

cd ${workdir}/01_bam_files

for i in $( ls *final.bam); do
	samtools flagstat ${workdir}/01_bam_files/${i} > ${workdir}/01_bam_files/${i%_final.bam}_flagstat.txt
done
