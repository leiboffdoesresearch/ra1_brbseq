#load info from config file
#config file contains sample name info, genome file info, and hisat parameters
#consider adding featureCounts summary factors as config parameters
configfile: "config.yaml"

#use minimal debian containerized environment with conda
#useful for OS standardization
container: "docker://continuumio/miniconda3:4.5.11"

#this rule looks for all the final files
#drives the back propagation of all intermediate rules
#elimintates the need for specifying inputs or outputs on command line
rule all:
    input:
        expand('{sample}_solo/features.tsv', sample=config["samples"])

############ FOR BRBseq ############
#sample barcode and UMI info is in R1
#pattern: (?P<cell_1>.{6})(?P<umi_1>.{10}[ACG]{5}).*
#real cDNA info comes from R2
#run trimming, alignment, etc on R2 only
####################################

############## Outline #############
#process barcodes, align, and quanitfy with STARsolo
####################################

#try STARsolo
# CB UMI Simple
rule STARsolo:
    input:
        R1 = 'compressed_reads/{sample}_R1.fastq.gz'
        R2 = 'compressed_reads/{sample}_R2.fastq.gz'
    output:
        star_out = '{sample}_STAR'
        solo_out = '{sample}_solo'
    conda:
        "envs/hisat2.yaml"
    threads: 40 # how many?
    params: 
        genome_index = config["genome_info"]["star_index"], #genome index
        whitelist = config["library_info"]["whitelist"] #barcode whitelist
    shell:
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {params.genome_index} "
        "--readFilesIn {input.R2} {input.R1} " #order required by solo
        "--readFilesCommand zcat "
        "--outFilterType BySJout "
        "--outFilterMultimapNmax 20 "
        "--outFilterMultimapNmax 20 "
        "--alignSJoverhangMin 8 "
        "--alignSJDBoverhangMin 1 "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFileNamePrefix {output.star_out} "
        "--outSAMtype BAM SortedByCoordinate "
        "--alignIntronMin 60 "
        "--alignIntronMax 50000 "
        "--soloBarcodeMate 0 "
        "--soloType CB_UMI_Simple "
        "--soloBarcodeReadLength 0 "
        "--soloCBstart 1   --soloCBlen 6 "
        "--soloUMIstart 7   --soloUMIlen 15 "
        "--soloCBwhitelist {params.whitelist} " 
        "--soloOutFileNames {output.solo_out} "
        "--soloStrand Reverse"