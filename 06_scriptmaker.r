# File for creating sh files and analyzing data
# adapted from JD Manthey code by JC Cooper on 19 May 2022
##########################################################

# popmap = individual base names of fastq files, one line per individual
# make sure reference is indexed with bwa and samtools before use, and use CreateSequenceDictionary in GATK
# the following was all run manually in R
# not run as a script

options(scipen=999)

library(data.table) # for like function

### Genera
# Chamaetylas
# Cossypha
# Melaenornis
# Phylloscopus

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
##############    START OPTIONS     #####################################

# project specific options
project_directory <- "~/scratch/albertine_rift/"
setwd(project_directory) # just to know where we are and read popmap
directory_name <- "11_genotype_lacustrine"
phyllo_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
								    'GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna')
platy_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
										  'GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna')
musci_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
										  'GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna')

cluster <- "bi"
output_name <- "lacustrine_genotype"
popmap <- "popmap.txt"
individuals <- read.table(popmap, sep="\t")
# faidx <- read.table("dryobates_chromosomes.fai", stringsAsFactors=F)
faidx_phyllo <- read.table(paste0(project_directory,
								  '00_fastq/references/GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna.fai'),
						   stringsAsFactors=F)
faidx_platy <- read.table(paste0(project_directory,
								  '00_fastq/references/GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna.fai'),
						   stringsAsFactors=F)
faidx_musci <- read.table(paste0(project_directory,
								  '00_fastq/references/GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna.fai'),
						   stringsAsFactors=F)

# we will load these with module and run using java
# singularity_cache <- "/lustre/work/jmanthey/singularity-cachedir"
# name_of_gatk_singularity_image <- "gatk_4.2.3.0.sif"

# define minimum and maximum genotyping job sizes
min_scaffold_size <- 950000
max_genotype_job_size <- 5000000
max_individual_genotype_job_size <- 100000000

# define number of cores for each step
ncores_step1 <- 8
ncores_step2 <- 8
ncores_step3 <- 12
# define memory per core for each step
mem_step1 <- "8G"
mem_step2 <- "8G"
mem_step3 <- "8G"

##############    END OPTIONS     #######################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

# calculate total java memory per job (= memory available * 0.9)
total_mem1 <- ceiling(ncores_step1 * as.numeric(strsplit(mem_step1, "G")[[1]]) * 0.9)
total_mem2 <- ceiling(ncores_step2 * as.numeric(strsplit(mem_step2, "G")[[1]]) * 0.9)
total_mem3 <- ceiling(ncores_step3 * as.numeric(strsplit(mem_step3, "G")[[1]]) * 0.9)

# make directories
dir.create(directory_name)
dir.create(paste0(directory_name, "/01_gatk_split"))
dir.create(paste0(directory_name, "/02b_gatk_database"))
dir.create(paste0(directory_name, "/03b_group_genotype_database"))

# subset the index
# FAIDX steps must be repeated for each reference genome
# faidx <- faidx[faidx[,2] >= min_scaffold_size, ]
faidx_phyllo=faidx_phyllo[faidx_phyllo[,2] >= min_scaffold_size, ]
faidx_platy=faidx_platy[faidx_platy[,2] >= min_scaffold_size, ]
faidx_musci=faidx_musci[faidx_musci[,2] >= min_scaffold_size, ]

# finds scaffolds too big to genotype at once

faidx_phyllo_keep <- faidx_phyllo[faidx_phyllo[,2] < max_individual_genotype_job_size,1:2]
faidx_phyllo_change <- faidx_phyllo[faidx_phyllo[,2] >= max_individual_genotype_job_size,1:2]

faidx_platy_keep <- faidx_platy[faidx_platy[,2] < max_individual_genotype_job_size,1:2]
faidx_platy_change <- faidx_platy[faidx_platy[,2] >= max_individual_genotype_job_size,1:2]

faidx_musci_keep <- faidx_musci[faidx_musci[,2] < max_individual_genotype_job_size,1:2]
faidx_musci_change <- faidx_musci[faidx_musci[,2] >= max_individual_genotype_job_size,1:2]

# paste the interval to use to each of the faidx objects
faidx_phyllo_keep=cbind(faidx_phyllo_keep,faidx_phyllo_keep[,1],rep(1,nrow(faidx_phyllo_keep)),faidx_phyllo_keep[,2])
faidx_platy_keep=cbind(faidx_platy_keep,faidx_platy_keep[,1],rep(1,nrow(faidx_platy_keep)),faidx_platy_keep[,2])
faidx_musci_keep=cbind(faidx_musci_keep,faidx_musci_keep[,1],rep(1,nrow(faidx_musci_keep)),faidx_musci_keep[,2])

faidx_phyllo_change=cbind(faidx_phyllo_change,rep("x",nrow(faidx_phyllo_change)))
faidx_platy_change=cbind(faidx_platy_change,rep("x",nrow(faidx_platy_change)))
faidx_musci_change=cbind(faidx_musci_change,rep("x",nrow(faidx_musci_change)))

faidx_manip=function(faidx_change){
	new_faidx_change=c()
	for(a in 1:nrow(faidx_change)) {
		a_rep <- faidx_change[a,]
		a_breaks <- floor(as.numeric(a_rep[1,2]) / 2)
		a_break1 <- c(a_rep[1,1], a_breaks, paste0(a_rep[1,1], ":1-", a_breaks), 1, a_breaks)
		a_break2 <- c(paste0(a_rep[1,1], "b"), as.numeric(a_rep[1,2]) - a_breaks, paste0(a_rep[1,1], ":", a_breaks + 1, "-", as.numeric(a_rep[1,2])), a_breaks + 1, as.numeric(a_rep[1,2]))
		new_faidx_change <- rbind(new_faidx_change, a_break1, a_break2)
	}
	return(new_faidx_change)
}

# note two don't have anything in the change category
# faidx_phyllo_change_new=faidx_manip(faidx_phyllo_change)
# faidx_platy_change_new=faidx_manip(faidx_platy_change)
faidx_musci_change_new=faidx_manip(faidx_musci_change)

col_rename=function(faidx){
	colnames(faidx)=c("id", "length", "interval", "start", "end")
	return(faidx)
}

faidx_phyllo_keep2=col_rename(faidx_phyllo_keep)
# faidx_phyllo_change=col_rename(faidx_phyllo_change)
faidx_platy_keep2=col_rename(faidx_platy_keep)
# faidx_platy_change=col_rename(faidx_platy_change)
faidx_musci_keep2=col_rename(faidx_musci_keep)
faidx_musci_change2=col_rename(faidx_musci_change_new)

# only muscicapa gets rbound

faidx_musci2=rbind(faidx_musci_keep2,faidx_musci_change2)
faidx_phyllo2=faidx_phyllo_keep2
faidx_platy2=faidx_platy_keep2

# characterize, literally
characterize=function(faidx){
	faidx[,3] <- as.character(faidx[,3])
	faidx[,1] <- as.character(faidx[,1])
	faidx[,2] <- as.numeric(faidx[,2])
	faidx[,4] <- as.numeric(faidx[,4])
	faidx[,5] <- as.numeric(faidx[,5])
	faidx <- na.omit(faidx)
	return(faidx)
}

faidx_phyllo2=characterize(faidx_phyllo2)
faidx_platy2=characterize(faidx_platy2)
faidx_musci2=characterize(faidx_musci2)

phyllo_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
                                           'GCA_002305835.1_ASM230583v1_genomic_Phylloscopus-trochilus.fna')
platy_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
                                          'GCA_013397175.1_ASM1339717v1_genomic_Platysteira-castanea.fna')
musci_reference_genome_location <- paste0(project_directory,'00_fastq/references/',
                                          'GCF_000247815.1_FicAlb1.5_genomic_Ficedula-albicollis.fna')

# create a helper file for each taxonomic grouping

for(i in 1:3){
  if(i==1){
    individuals.x=individuals[which(individuals$V1%like%'Chamaetylas'|
                                      individuals$V1%like%'Cossypha'|
                                      individuals$V1%like%'Melaenornis'),]
    faidx=faidx_musci2
    fname='musci'
  }
  if(i==2){
    individuals.x=individuals[which(individuals$V1%like%'Phylloscopus'),]
    faidx=faidx_phyllo2
    fname='phyllo'
  }
  if(i==3){
    individuals.x=individuals[which(individuals$V1%like%'Batis'|
                                      individuals$V1%like%'Platysteira'),]
    faidx=faidx_platy2
    fname='platy'
  }
  individuals.x=as.data.frame(individuals.x)
  for(a in 1:nrow(individuals.x)){
    if(a == 1) {
      helper1 <- cbind(faidx[,1], as.character(faidx[,3]), rep(individuals.x[a,1], nrow(faidx)))
    } else {
      helper1 <- rbind(helper1, cbind(faidx[,1], as.character(faidx[,3]), rep(individuals.x[a,1], nrow(faidx))))
    }
  }
  # write the chromosome/scaffold to genotype for each job
  write(helper1[,2], file=paste0(directory_name, "/01_gatk_split/",fname,"_helper1.txt"), ncolumns=1)
  # write the individual to genotype for each job
  write(helper1[,3], file=paste0(directory_name, "/01_gatk_split/",fname,"_helper2.txt"), ncolumns=1)
  # write the chromosome name for each job
  write(helper1[,1], file=paste0(directory_name, "/01_gatk_split/",fname,"_helper1b.txt"), ncolumns=1)
  
  # write the helper files for the 2nd genotyping step
  # write the chromosome/scaffold to database for each job
  write(faidx[,1], file=paste0(directory_name, "/02b_gatk_database/",fname,"_helper3.txt"), ncolumns=1)
  write(faidx[,3], file=paste0(directory_name, "/02b_gatk_database/",fname,"_helper3b.txt"), ncolumns=1)
  
  # write the helper files for the 3rd genotyping step
  helper <- c()
  job_suffixes <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
  for(a in 1:nrow(faidx)) {
    a_rep <- faidx[a,]
    job_number <- ceiling(a_rep[1,2] / max_genotype_job_size)
    job_number <- job_suffixes[1:job_number]
    
    if(length(job_number) > 1) {
      # define interval to group genotype if more than one job
      interval_start <- a_rep[1,4]
      interval_end <- a_rep[1,4] + max_genotype_job_size - 1
      for(b in 1:length(job_number)) {
        helper <- rbind(helper, c(faidx[a,1], 
                                  paste0(strsplit(faidx[a,3], ":")[[1]][1], ":", interval_start, "-", interval_end),
                                  paste0(faidx[a,1], "__", job_number[b])
        ))
        interval_start <- interval_start + max_genotype_job_size
        if((b+1) != length(job_number)) {
          interval_end <- interval_end + max_genotype_job_size
        } else {
          interval_end <- faidx[a,5]
        }
      }
    } else {
      helper <- rbind(helper, c(faidx[a,1], 
                                paste0(faidx[a,1], ":", 1, "-", faidx[a,2]),
                                paste0(faidx[a,1])
      ))
    }
  }
  # write the chromosome/scaffold to genotype for each job
  write(helper[,1], file=paste0(directory_name, "/03b_group_genotype_database/",fname,"_helper4.txt"), ncolumns=1)
  # write the interval range to genotype for each job
  write(helper[,2], file=paste0(directory_name, "/03b_group_genotype_database/",fname,"_helper5.txt"), ncolumns=1)
  # write the output base name for each job
  write(helper[,3], file=paste0(directory_name, "/03b_group_genotype_database/",fname,"_helper6.txt"), ncolumns=1)
}

########################################

# create single script code to run all individuals of a specific group
### NOTE ###

# I moved all BAM files into three subfolders, one for each taxonomic grouping
# this folder is then referenced below

for(i in 1:3){
  if(i==1){
    fname='musci'
    reference_genome_location=musci_reference_genome_location
  }
  if(i==2){
    fname='phyllo'
    reference_genome_location=phyllo_reference_genome_location
  }
  if(i==3){
    fname='platy'
    reference_genome_location=platy_reference_genome_location
  }
  
  ########################################
  
  # step 1
  # genotype all individuals using GATK, one array job using the above two helper files
  
  a.script <- paste0(directory_name, "/01_gatk_split/",fname,"_step1_array.sh")
  write("#!/bin/sh", file=a.script)
  write("#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/", file=a.script, append=T)
  write(paste0("#SBATCH --job-name=", "step1_",fname), file=a.script, append=T)
  write("#SBATCH --nodes=1", file=a.script, append=T)
  write(paste0("#SBATCH --ntasks=", ncores_step1), file=a.script, append=T)
  write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
  write("#SBATCH --time=15-00:00:00", file=a.script, append=T)
  write(paste0("#SBATCH --mem-per-cpu=", mem_step1), file=a.script, append=T)
  write(paste0("#SBATCH --array=1-", nrow(helper1)), file=a.script, append=T)
  write("", file=a.script, append=T)
  write("module load java", file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper1.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("ind_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper2.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("name_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper1b.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  
  #gatk 4.0
  a_name <- paste0(project_directory,"/01_bam_files/",fname,"/${ind_array}","_final.bam")
  # alter to reflect different reference genomes
  gatk_command <- paste0('~/work/gatk-4.2.3.0/gatk gatk --java-options "-Xmx', 
                         total_mem1, 'g" HaplotypeCaller -R ', reference_genome_location, 
                         " -I ", a_name, " -ERC GVCF -O ", project_directory, 
                         "/02_vcf/",fname, "/${name_array}", "._${ind_array}_.g.vcf", 
                         " --QUIET --intervals ", "${chr_array}")
  write(gatk_command, file=a.script, append=T)

  #######################################
  
  # step 2
  # create genotyping database for each of the chromosomes
  a.script <- paste0(directory_name, "/02b_gatk_database/",fname,"_step2_array.sh")
  
  write("#!/bin/sh", file=a.script)
  write("#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/", file=a.script, append=T)
  write(paste0("#SBATCH --job-name=", "step2_",fname), file=a.script, append=T)
  write("#SBATCH --nodes=1", file=a.script, append=T)
  write(paste0("#SBATCH --ntasks=", ncores_step2), file=a.script, append=T)
  write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
  write("#SBATCH --time=48:00:00", file=a.script, append=T)
  write(paste0("#SBATCH --mem-per-cpu=", mem_step2), file=a.script, append=T)
  write(paste0("#SBATCH --array=1-", nrow(faidx)), file=a.script, append=T)
  write("", file=a.script, append=T)
  write("module load java", file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("name_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper3.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper3b.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  
  #make list of all vcfs to database
  for(b in 1:nrow(individuals)) {
    if(b == 1) {
      vcf_total <- paste0("-V ",project_directory,"/02_vcf/",fname,"/${name_array}","._",individuals[b,1], "_.g.vcf")
    } else {
      vcf_total <- paste0(vcf_total, " -V ", project_directory, "/02_vcf/",
                          fname,"/${name_array}", "._", individuals[b,1], "_.g.vcf")
    }
  }
  
  #gatk 4.0
  gatk_command=paste0('~/work/gatk-4.2.3.0/gatk gatk --java-options "-Xmx',
                      total_mem2,'g" GenomicsDBImport --genomicsdb-shared-posixfs-optimizations ', 
                      vcf_total, " --genomicsdb-workspace-path ", 
                      project_directory, "/02_vcf/",fname,"/${name_array}", " -L ", "${interval_array}")
  
  write(gatk_command, file=a.script, append=T)
  
  #################################
  
  # step 3
  # group genotyping for each interval
  a.script <- paste0(directory_name, "/03b_group_genotype_database/",fname,"_step3_array.sh")
  write("#!/bin/sh", file=a.script)
  write("#SBATCH --chdir=/panfs/pfs.local/home/j141c380/work/", file=a.script, append=T)
  write(paste0("#SBATCH --job-name=", "step3_",fname), file=a.script, append=T)
  write("#SBATCH --nodes=1", file=a.script, append=T)
  write(paste0("#SBATCH --ntasks=", ncores_step3), file=a.script, append=T)
  write(paste0("#SBATCH --partition=", cluster), file=a.script, append=T)
  write("#SBATCH --time=48:00:00", file=a.script, append=T)
  write(paste0("#SBATCH --mem-per-cpu=", mem_step3), file=a.script, append=T)
  write(paste0("#SBATCH --array=1-", nrow(helper)), file=a.script, append=T)
  write("", file=a.script, append=T)
  write("module load java", file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper4.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper5.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  write(paste0("name_array=$( head -n${SLURM_ARRAY_TASK_ID} ",fname,"_helper6.txt | tail -n1 )"), file=a.script, append=T)
  write("", file=a.script, append=T)
  
  #gatk 4.0
  gatk_command <- paste0('~/work/gatk-4.2.3.0/gatk gatk --java-options "-Xmx', 
                         total_mem3, 'g" GenotypeGVCFs --genomicsdb-shared-posixfs-optimizations -R ', 
                         reference_genome_location, " -V gendb://", project_directory, 
                         "/02_vcf/",fname,"/${chr_array}", " --include-non-variant-sites -O ", 
                         project_directory, "/03_vcf/",fname, "/${name_array}", ".g.vcf", " -L ", "${interval_array}")
  write(gatk_command, file=a.script, append=T)
}
