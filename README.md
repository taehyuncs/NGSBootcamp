# NGSBootcamp
NGS Bootcamp for whole exome and genome sequencing analysis

Single Neucleotide Variants (SNVs) and small Insertion and Deletion (Indel) analysis

# Very general workflow for variant analysis
# 1. fastq files from your sample
# 2. Run FastQC for quality control
#  2a) if qualtify of bases (in particular at the end of reads) are bad, perform trimming.
#  2b) if there are adapter sequences present in your fastq files, remove adapter sequence.
# => OK. fastq quality is good. Lets go to step 3
# 3. Perform Alignemnt using your favorite aligner (e.g., BWA, etc.) and generate SAM or BAM files.
# 4. Sort bam file and generate an index file for your bam file.
# 5. Variant calling using VarScan (http://dkoboldt.github.io/varscan/) Please read manual (http://dkoboldt.github.io/varscan/using-varscan.html) and their original publication (http://www.ncbi.nlm.nih.gov/pubmed/22300766)
#  5a) Generate mpileup data.
#    => What is (m)pileup data? Go to the following link: http://samtools.sourceforge.net/pileup.shtml
#       mpileup data will be used to call SNVs and Indels using VarScan
#  5b) Call variants using VarScan with mpileup data
# 6. Now we have variant calls (e.g., a list of SNVs and Indels). We want to annotate them to see whether those variants are reported by someone eles or associated with known disease phenotype and/or see what would be functional impact.
#  6a) We will use snpEff tool (http://snpeff.sourceforge.net) and please read snpEff readme file (http://snpeff.sourceforge.net/SnpEff.html)
# 7. Interpretation. Now we have annotated variant calls and need to make sense out of it.
# 8. It is important to visualize your bam file to double check whether variants called by computational algorithms are real(?).
# Utilmately, you want to validate your variants using Sanger or etc.



########### Part 1: Variant calling using GATK #############

##########################################################################################
#         GATK Variant Call                                                              #
# https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS #
##########################################################################################

gatk=/home/thwang/software/GenomeAnalysisTK.jar #GATK is a computational tool to call variants
picard=/home/thwang/software/picard.jar #a tool for post processing for alignment files.

samtools=/home/thwang/software/samtools-1.2/bin/samtools #Samtools
reference=/home/thwang/reference/hg19.fa #Reference fastq file
java=/data/bootcamp/seqprg/jre1.8.0_101/bin/java #Java
varscan=/home/thwang/software/VarScan.v2.3.9.jar #VarScan is another computational tool to call variants
snpEFF=/home/thwang/software/snpEff/snpEff.jar #snpEff is a tool for variant annotation
Nbam=/home/thwang/test/normal.bam # Alignment file from normal tissue 
Tbam=/home/thwang/test/tumor.bam # Alignment file from tumor tissue 
outdir=/home/thwang/test/new # Output dir

sample=tumor
outdir=/home/thwang/test/new


$java -jar ${picard} SortSam INPUT=${outdir}/${sample}.sam OUTPUT=${outdir}/${sample}_sorted_reads.bam SORT_ORDER=coordinate
$java -jar ${picard} MarkDuplicates INPUT=${outdir}/${sample}_sorted_reads.bam OUTPUT=${outdir}/${sample}_dedup_reads.bam METRICS_FILE=metrics.txt
$samtools index ${outdir}/${sample}_dedup_reads.bam
$java -jar ${gatk} -T RealignerTargetCreator -R /home/thwang/reference/hg19.fa -I ${outdir}/${sample}_dedup_reads.bam -known /home/thwang/reference/1000G_phase1.indels.b37.vcf -known /home/thwang/reference/Mills_and_1000G_gold_standard.indels.b37.vcf -o ${outdir}/${sample}_target_intervals.list -L 17:7552323-7616951
$java -jar ${gatk} -I ${outdir}/${sample}_dedup_reads.bam -R /home/thwang/reference/hg19.fa -T IndelRealigner -targetIntervals ${outdir}/${sample}_target_intervals.list -known /home/thwang/reference/1000G_phase1.indels.b37.vcf -known /home/thwang/reference/Mills_and_1000G_gold_standard.indels.b37.vcf -o ${outdir}/${sample}_realigned.bam -L 17:7552323-7616951
$java -jar ${gatk} -T BaseRecalibrator -R /home/thwang/reference/hg19.fa -I ${outdir}/${sample}_realigned.bam -knownSites /home/thwang/reference/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites /home/thwang/reference/1000G_phase1.indels.b37.vcf -knownSites /home/thwang/reference/dbsnp_138.b37.vcf -o ${outdir}/${sample}_recal.table -L 17:7552323-7616951
$java -jar ${gatk} -T PrintReads -R /home/thwang/reference/hg19.fa -I ${outdir}/${sample}_realigned.bam -BQSR ${outdir}/${sample}_recal.table -o ${outdir}/${sample}_recal.bam -L 17:7552323-7616951

java -Djava.io.tmpdir=/qbrc/home/ygao/tmp -jar /qbrc/home/ygao/reference_known/GenomeAnalysisTK.jar -T HaplotypeCaller -R /qbrc/home/ygao/reference_known/hg19.fa -I ${outdir}/${sample}_recal.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ${outdir}/${sample}.vcf

$java -jar ${SnpSift} annotate ${snpEFFanno}/CosmicCodingMuts.vcf ${outdir}/${sample} > ${outdir}/${sample}.tmp.COSMIC.coding.snpEFF
$java -jar ${SnpSift} annotate ${snpEFFanno}/CosmicNonCodingVariants.vcf ${outdir}/${sample}.tmp.COSMIC.coding.snpEFF > ${outdir}/${sample}.tmp.COSMIC.coding.noncoding.snpEFF
$java -jar ${SnpSift} annotate ${snpEFFanno}/dbsnp/All_20160408.vcf ${outdir}/${sample}.tmp.COSMIC.coding.snpEFF > ${outdir}/${sample}.tmp.COSMIC.coding.noncoding.dbsnp.snpEFF
$java -jar ${SnpSift} dbnsfp -v ${snpEFFanno}/dbNSFP.txt.gz ${outdir}/${sample}.tmp.COSMIC.coding.noncoding.dbsnp.snpEFF > ${outdir}/${sample}.tmp.COSMIC.coding.noncoding.dbsnp.dnNSFP.snpEFF
$java -jar ${SnpSift} gwasCat ${snpEFFanno}/gwascatalog.txt ${outdir}/${sample}.tmp.COSMIC.coding.noncoding.dbsnp.dnNSFP.snpEFF > ${outdir}/${sample}.multi.snpEFF

########### Part 2: Variant calling using VarScan #############

#######################################
#         VarScan Variant Call        #
# http://dkoboldt.github.io/varscan/  #
#######################################

#set environment variables for software and input file
samtools=/home/thwang/software/samtools-1.2/bin/samtools #Samtools
reference=/home/thwang/reference/hg19.fa #human reference data. We use hg19 version not hg38
java=/data/bootcamp/seqprg/jre1.8.0_101/bin/java #Java
varscan=/home/thwang/software/VarScan.v2.3.9.jar #VarScan is software to call somatic mutation
snpEFF=/home/thwang/software/snpEff/snpEff.jar #snpEff is software to annotate variant calls based on VCF file
Nbam=/home/thwang/test/normal.bam #Tumor alignment bam file by bwa (https://sourceforge.net/projects/bio-bwa/)
Tbam=/home/thwang/test/tumor.bam #Normal alignment bam file by bwa (https://sourceforge.net/projects/bio-bwa/)

#generate mpileup using samtools. Runing time ~ 5min for each sample
$samtools mpileup -f $reference $Nbam > normal.mpileup
$samtools mpileup -f $reference $Tbam > tumor.mpileup

#calling somatic mutations using varscan, the output would be somatic.snp.vcf and somatic.indel.vcf. Runing time ~ 5min
$java -jar $varscan somatic normal.mpileup tumor.mpileup somatic --output-vcf

#annotate vcf file using snpEFF. Runing time ~ 2min
$java -jar $snpEFF hg19 somatic.indel.vcf > somatic.indel.ann.vcf
$java -jar $snpEFF hg19 somatic.snp.vcf > somatic.snp.ann.vcf
### when testing in the account of thwang, the time was right, but when I was using my account, the time would be very long....(more than 10min for steps which suppose to use 5min)​
