#First step is to generate the Cellranger output from raw fastq data. Follow the tutorials on sratoolkit to download and split the fastqs
#Then proceed to 10X tutorials on Cellranger - you need to install it and generate the output with cellranger count 
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

#Install samtools and velocyto

samtools sort -t CB -O BAM -o/mnt/с/Users/Emil/10X/Baoyi_gp_test/10X_20_055/outs/cellsorted_possorted_genome_bam.bam 
      /mnt/C/Users/Emil/10X/Baoyi_gp_test/10X_20_055/outs/possorted_genome_bam.bam #generate the cellsorted.bam file. If you have enough RAM, you may skip this step
      
#Follow velocyto tutorial how to obtain the mask file and reference genome annotation
#The appropriate genome annotation version could be found in web_summary file from Cellranger output folder
velocyto run10x -m /mnt/c/Users/Emil/10X/Baoyi_gp_test/mm10_rmsk.gtf 
/mnt/c/Users/Emil/10X/Baoyi_gp_test/10X_21_055 /mnt/c/Users/Emil/10X/Baoyi_gp_test/gencode.vM27.annotation.gtf

#Now the .loom file should be generated in a new folder 'velocyto' in Cellranger output folder. This file will be used for the RNA Velocity analysis
