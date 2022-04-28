#!/bin/bash
#SBATCH --job-name=02_align_albertine_rift      # Job name
#SBATCH --partition=bi          # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jccooper@ku.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12		# Threads
#SBATCH --mem-per-cpu=4g            # Job memory request
#SBATCH --time=15-00:00:00 	# time in days-hours:min:sec
#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/ 	# Directory
#SBATCH --output=02_align_albertine_rift.log     # Standard output and error log
#SBATCH --error=02_align_albertine_rift.err # error file

# aligning script

# intel ??

module load java bwa samtools

#export SINGULARITY_CACHEDIR="${workdir}/singularity-cache"

# define main working directory
#workdir=${workdir}/albertine_rift

# word DIR doesn't work - space in name of hard drive!

cd ~/work/albertine_rift/

# save filelist
# ls 00_fastq/*_R1.fastq.gz > basenames.txt
# edited this slightly by hand! removed the tail end of the phrase

export SINGULARITY_CACHEDIR="${workdir}/singularity_cachedir"
workdir=/home/j141c380/work

# get list of files to be analyzed

# define the reference genome
refplatysteira=${workdir}/albertine_rift/00_fastq/references/GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna
refphylloscopus=${workdir}/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna
refficedula=${workdir}/albertine_rift/00_fastq/references/GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna

# writing as loop per file
# each has a different reference file

# Phylloscopus

LINES=$(cat basenames_phyllo.txt)
# for debugging
# LINE=Phylloscopus_laetus_FMNH346427
refgenome=${workdir}/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna

for LINE in $LINES
do
	# run bbduk
	# AUTOMATICALLY USES 7 THREADS
	${workdir}/bbmap/bbduk.sh \
		in1=${workdir}/albertine_rift/00_fastq/${LINE}_R1.fastq.gz \
		in2=${workdir}/albertine_rift/00_fastq/${LINE}_R2.fastq.gz \
		out1=${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		out2=${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
		ref=${workdir}/bbmap/resources/adapters.fa hdist=1 tbo tpe
	
	# run bwa mem
	# specifies 12 cores
	bwa mem -t 12 ${workdir}/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		> ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
		
	# convert sam to bam
	samtools view -b -S -o ${workdir}/albertine_rift/01_bam_files/${LINE}.bam ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# remove sam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# clean up the bam file
	../gatk-4.2.3.0/gatk CleanSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam

	# remove the raw bam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.bam
	
	# sort the cleaned file
	../gatk-4.2.3.0/gatk SortSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam --SORT_ORDER coordinate
	
	# remove cleaned bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam
	
	# add read groups to sorted and cleaned bam file
	../gatk-4.2.3.0/gatk AddOrReplaceReadGroups -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${LINE}
	
	# remove cleaned and sorted bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam
	
	# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
	../gatk-4.2.3.0/gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 \
		-M ${workdir}/albertine_rift/01_bam_files/${LINE}_markdups_metric_file.txt -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
	
	# remove sorted, cleaned, and read grouped bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam
	
	# index the final bam file
	samtools index ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
done

# Batis

LINES=$(cat basenames_batis.txt)
# for debugging
# LINE=Phylloscopus_laetus_FMNH346427
refgenome=${workdir}/albertine_rift/00_fastq/references/GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna

for LINE in $LINES
do
	# run bbduk
	# AUTOMATICALLY USES 7 THREADS
	${workdir}/bbmap/bbduk.sh \
		in1=${workdir}/albertine_rift/00_fastq/${LINE}_R1.fastq.gz \
		in2=${workdir}/albertine_rift/00_fastq/${LINE}_R2.fastq.gz \
		out1=${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		out2=${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
		ref=${workdir}/bbmap/resources/adapters.fa hdist=1 tbo tpe
	
	# run bwa mem
	# specifies 12 cores
	bwa mem -t 12 ${workdir}/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		> ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
		
	# convert sam to bam
	samtools view -b -S -o ${workdir}/albertine_rift/01_bam_files/${LINE}.bam ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# remove sam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# clean up the bam file
	../gatk-4.2.3.0/gatk CleanSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam

	# remove the raw bam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.bam
	
	# sort the cleaned file
	../gatk-4.2.3.0/gatk SortSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam --SORT_ORDER coordinate
	
	# remove cleaned bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam
	
	# add read groups to sorted and cleaned bam file
	../gatk-4.2.3.0/gatk AddOrReplaceReadGroups -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${LINE}
	
	# remove cleaned and sorted bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam
	
	# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
	../gatk-4.2.3.0/gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 \
		-M ${workdir}/albertine_rift/01_bam_files/${LINE}_markdups_metric_file.txt -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
	
	# remove sorted, cleaned, and read grouped bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam
	
	# index the final bam file
	samtools index ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
done

# Muscicapidae

LINES=$(cat basenames_musci.txt)
# for debugging
# LINE=Phylloscopus_laetus_FMNH346427
refgenome=${workdir}/albertine_rift/00_fastq/references/GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna

for LINE in $LINES
do
	# run bbduk
	# AUTOMATICALLY USES 7 THREADS
	${workdir}/bbmap/bbduk.sh \
		in1=${workdir}/albertine_rift/00_fastq/${LINE}_R1.fastq.gz \
		in2=${workdir}/albertine_rift/00_fastq/${LINE}_R2.fastq.gz \
		out1=${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		out2=${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
		ref=${workdir}/bbmap/resources/adapters.fa hdist=1 tbo tpe
	
	# run bwa mem
	# specifies 12 cores
	bwa mem -t 12 ${workdir}/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R1.fastq.gz \
		${workdir}/albertine_rift/01_cleaned/${LINE}_R2.fastq.gz \
		> ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
		
	# convert sam to bam
	samtools view -b -S -o ${workdir}/albertine_rift/01_bam_files/${LINE}.bam ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# remove sam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.sam
	
	# clean up the bam file
	../gatk-4.2.3.0/gatk CleanSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam

	# remove the raw bam
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}.bam
	
	# sort the cleaned file
	../gatk-4.2.3.0/gatk SortSam -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam -O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam --SORT_ORDER coordinate
	
	# remove cleaned bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned.bam
	
	# add read groups to sorted and cleaned bam file
	../gatk-4.2.3.0/gatk AddOrReplaceReadGroups -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${LINE}
	
	# remove cleaned and sorted bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted.bam
	
	# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
	../gatk-4.2.3.0/gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 \
		-M ${workdir}/albertine_rift/01_bam_files/${LINE}_markdups_metric_file.txt -I ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam \
		-O ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
	
	# remove sorted, cleaned, and read grouped bam file
	rm ${workdir}/albertine_rift/01_bam_files/${LINE}_cleaned_sorted_rg.bam
	
	# index the final bam file
	samtools index ${workdir}/albertine_rift/01_bam_files/${LINE}_final.bam
done
