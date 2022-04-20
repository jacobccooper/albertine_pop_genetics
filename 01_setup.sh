cd /media/j141c380/My\ Passport/Sequences/albertine_rift/

# make all directories 
mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 10_align_script
mkdir 12_filter
mkdir 04_relernn
mkdir 05_trees100kbp
mkdir 05_trees100kbp/windows
mkdir 06_admixture
mkdir 07_admixture_windows
mkdir 07_admixture_windows/windows

# move all the raw fastq data to the 00_fastq directory with a file transfer program

# change directory to the raw data directory
cd /media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq

# skipping cat steps from Manthey
# I only have two reads per sample