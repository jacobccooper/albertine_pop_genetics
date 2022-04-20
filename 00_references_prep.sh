# this will help get one member per genus to a reference genome
# relies on https://speciationgenomics.github.io/mapping_reference/
# based on code by JD Manthey

# set directory
cd /media/j141c380/My Passport/Sequences/albertine_rift/00_fastq/references/

# reference is GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna.gz

# https://github.com/broadinstitute/picard
# java -jar /home/j141c380/Documents/picard/build/libs/picard.jar 

gunzip GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna.gz

samtools faidx GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna

bwa index GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna

# note that a space in the path name for "My Passport" throws an error

java -jar /home/j141c380/Documents/picard/build/libs/picard.jar \
	CreateSequenceDictionary \
	R=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna \
	O=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.dict

##### Platysteira

gunzip GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna.gz

samtools faidx GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna

bwa index GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna

# note that a space in the path name for "My Passport" throws an error

java -jar /home/j141c380/Documents/picard/build/libs/picard.jar \
	CreateSequenceDictionary \
	R=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna \
	O=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.dict

##### Phylloscopus

gunzip GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna.gz

samtools faidx GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna

bwa index GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna

# note that a space in the path name for "My Passport" throws an error

java -jar /home/j141c380/Documents/picard/build/libs/picard.jar \
	CreateSequenceDictionary \
	R=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna \
	O=/media/j141c380/My\ Passport/Sequences/albertine_rift/00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.dict

