#!/bin/bash
# MID analysis using fgbio
# S Sandhu 160911

set -e
set -x
#gunzip *gz
#for f in *_R1_001.fastq.gz
#18-482-2_S6_L001_R1_001.fastq.gz  HD-0001_S2_L001_R1_001.fastq.gz  HD-001_S3_L001_R1_001.fastq.gz  HD-005_S4_L001_R1_001.fastq.gz  HD-WT_S1_L001_R1_001.fastq.gz
for f in HD-001_S3_L001_R1_001.fastq.gz HD-005_S4_L001_R1_001.fastq.gz HD-WT_S1_L001_R1_001.fastq.gz
do
    fq1=$f
    fq2=${f%_R1*}_R2_001.fastq.gz
    OUT=${fq1%%_L001*}

    #anno=/tools/snpEff
    #fgbio=/home/sandhu/tools/fgbio/target/scala-2.11/fgbio-0.1.4-SNAPSHOT.jar
    #picard=/seq/tools/picard/build/libs/picard.jar
    #trim=/tools/trimmomatic36

    dbnsfp=/data/snpEff/GRCh37.75/dbNSFP2.9.txt.gz
    cosmic=/data/snpEff/GRCh37.75/CosmicCodingMuts.vcf
    cosmicNonCod=/data/snpEff/GRCh37.75/CosmicNonCodingVariants.vcf
    clinvar=/data/snpEff/GRCh37.75/clinvar_20170516.vcf
    gatk_data=/data/gatk

    REF="/data/refgenomes/hg19.with.mt.fasta"
    BEDDIR="/data/bedfiles"
    #GATK="java -Xmx64g -jar /tools/GenomeAnalysisTK-2014.4-2-g9ad6aa8/GenomeAnalysisTK.jar"
    MASTER="${BEDDIR}/18-2132_EGFR_MID_Masterfile.txt"
    BED8="${BEDDIR}/egfrMID_8col.bed"  ## for better variant calling with Vardict in amplicon mode
    # make a (nonmerged) bed from  Masterfile

    awk '{print $1,$2,$3,$4}' OFS="\t" $MASTER > nonmerged_targets.bed

    # make a merged bedfile from nonmerged
    bedtools merge -nms -i nonmerged_targets.bed | sed 's/;.*//' > merged_targets.tmp
    awk '{print $1,$2,$3,"+",$4}' OFS="\t" merged_targets.tmp > merged_targets.bed

    ### overlap-removed non-merged target BED ###
    # get all overlapping regions for target bedfile

    bedtools intersect \
        -a nonmerged_targets.bed -b nonmerged_targets.bed |
        awk 'a[$1FS$2FS$3]++' OFS="\t" > overlapped_regions.bed
    bedtools subtract \
        -a nonmerged_targets.bed -b overlapped_regions.bed \
     > nonmerged_noolaps_targets.tmp
    awk '{print $1,$2,$3,"+",$4}' OFS="\t" nonmerged_noolaps_targets.tmp > nonmerged_noolaps_targets.bed

    #mv overlapped_regions.bed bed
    BED=merged_targets.bed
    NOOLAPBED=nonmerged_noolaps_targets.bed


    ########## PRE-MID analysis & POST-MID with fgbio ##############
    #generate index sequence by keep 1st 5 bases of each read
    trimmomatic SE -threads 20 $fq2 $OUT.R2.MID.fq.gz \
     CROP:10

    # MID fastq file
    mid=$OUT.R2.MID.fq.gz

    ######## crop 1st 10bp MID
    trimmomatic SE -threads 20 $fq2 $OUT.trimd.R2.fq.gz \
        HEADCROP:10

    ######## Adapter Trim, & crop 1st 5bp MID
    #java -jar $trim/trimmomatic-0.36.jar PE -threads 20 $fq1 $fq2 $OUT.trimd.R1.fq.gz \
     #   $OUT_unpairedR1.fq.gz $OUT.trimd.R2.fq.gz $OUT_unpairedR2.fq.gz \
      #  ILLUMINACLIP:${trim}/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:10


    #echo "***********bwa align 2s to hg19 and collect picard metrics********"
    #align using bwa with verbosity 1 (to print errors only), mark sec alignments and add RG
    bwa mem $REF $fq1 $OUT.trimd.R2.fq.gz -M -t 16 -v 1 > $OUT.sam

    # Sort, Add RG for nopclip bam
    picard AddOrReplaceReadGroups I=$OUT.sam O=$OUT.preMID.bam \
     SO=coordinate VALIDATION_STRINGENCY=STRICT RGID=Accel RGLB=MID RGSM=$OUT RGPL=Illumina \
     RGPU=Miseq CREATE_INDEX=TRUE

    ## try remove a from bwa to avoid * error with gatk
    # query sort for primerclip
    picard SortSam I=$OUT.sam O=$OUT.qsort.sam SO=queryname

    #PrimerCLIP before MID
    #/home/sandhu/.local/bin/primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam
    primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam

    # Sort
    picard AddOrReplaceReadGroups I=preM${OUT}.preMID.pclip.sam O=$OUT.preMID.pclip.bam \
     SO=coordinate VALIDATION_STRINGENCY=STRICT RGID=Accel RGLB=MID RGSM=$OUT RGPL=Illumina \
     RGPU=Miseq CREATE_INDEX=TRUE

    fgbio -Xms500m -Xmx64g AnnotateBamWithUmis -i $OUT.preMID.bam -f $mid -o $OUT.bamAnnotdWumi.bam >& log.${f%_L001*}.fgbioAnnoBam.txt

    #RevertSam to sanitise
    picard RevertSam I=$OUT.bamAnnotdWumi.bam O=$OUT.sanitised.bam \
        SANITIZE=true REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false

    #SetMateInfo to bring MQ tags
    fgbio -Xmx64g -Xms500m  SetMateInformation -i $OUT.sanitised.bam -o $OUT.sanitised.setmateinfo.bam
    #java -jar $fgbio SetMateInformation -i $OUT.pclip.sanitised.bam -o $OUT.pclip.sanitised.setmateinfo.bam

    #sort by queryname
    picard SortSam I=$OUT.sanitised.setmateinfo.bam \
        O=$OUT.sanitised.setmateinfo.qsort.bam SO=queryname
    #java -jar $picard SortSam I=$OUT.pclip.sanitised.setmateinfo.bam \
    #    O=$OUT.pclip.sanitised.setmateinfo.qsort.bam SO=queryname

    #Group Reads by UMI
    fgbio -Xms500m -Xmx64g GroupReadsByUmi -s adjacency --edits 1 \
        -i $OUT.sanitised.setmateinfo.qsort.bam -o $OUT.groupdbyMID.bam >& log.${f%_L001*}.fgbioGroupbyMID.txt
    #java -jar $fgbio GroupReadsByUmi -s adjacency --edits 1 \
    #    -i $OUT.pclip.sanitised.setmateinfo.qsort.bam -o $OUT.pclip.groupdbyMID.bam

    #Make consensus N molecules to make consensus Out is coordinate sorted
    fgbio -Xms500m -Xmx64g CallMolecularConsensusReads -M 3 \
        -i $OUT.groupdbyMID.bam -o $OUT.M3.consensus.bam -r $OUT.M3.notused4consensus.bam >& log.${f%_L001*}.fgbioCMCR-M3.txt

    fgbio -Xms500m -Xmx64g CallMolecularConsensusReads -M 1 \
        -i $OUT.groupdbyMID.bam -o $OUT.consensus.bam -r $OUT.notused4consensus.bam >& log.${f%_L001*}.fgbioCMCR.txt

    #java -jar $fgbio CallMolecularConsensusReads -M 3 \
    #    -i $OUT.pclip.groupdbyMID.bam -o $OUT.M3.pclip.consensus.bam -r $OUT.M3.pclip.notused4consensus.bam

    #bamtofastq
    picard SamToFastq I=$OUT.M3.consensus.bam \
        F=$OUT.consensus.R1.fq F2=$OUT.consensus.R2.fq FU=$OUT.unpaired.fq
    #java -Xmx128g -jar $picard SamToFastq I=$OUT.M3.pclip.consensus.bam \
    #    F=$OUT.pclip.consensus.R1.fq F2=$OUT.pclip.consensus.R2.fq FU=$OUT.unpaired.fq

    #Realign
    bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M > ${OUT}_consensus.sam
    #bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M -I 100,50,50,500 > ${OUT}_consensus.sam

    ## qsort for primer-clip
    picard SortSam I=${OUT}_consensus.sam \
        O=$OUT.consensus.qsort.sam SO=queryname
    #java -jar $picard SortSam I=${OUT}_pclip.consensus.sam \
    #    O=$OUT.pclip.consensus.qsort.sam SO=queryname

    #PrimerCLIP
    primerclip $MASTER $OUT.consensus.qsort.sam ${OUT}.fgbio.pclip.sam

    #sort and Add read groups
    picard AddOrReplaceReadGroups I=${OUT}.fgbio.pclip.sam \
        O=$OUT.fgbio.bam RGID=MID-Amp RGLB=$OUT RGSM=NA12878 RGPL=Illumina \
        RGPU=MiSeq SO=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=STRICT

    # make BED file for picard
    cd=$PWD
    cp $BED $cd/BED.bed
    samtools view -H $OUT.fgbio.bam > header.txt
    cat header.txt BED.bed > BED.picard.bed
    pBED=BED.picard.bed
    # Collect TargetedPCRMetrics from the TWO BAMs #######
    picard CollectTargetedPcrMetrics I=$OUT.preMID.bam O=$OUT.preMID.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF PER_TARGET_COVERAGE=$OUT.preMID.perTargetCov.txt

    picard CollectTargetedPcrMetrics I=$OUT.fgbio.bam O=$OUT.fgbio.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF PER_TARGET_COVERAGE=$OUT.fgbio.perTargetCov.txt


    #BEDtools
    coverageBed -abam $OUT.fgbio.bam -b $BED > $OUT.fgbio.cov
    coverageBed -abam $OUT.fgbio.bam -b $BED -d > $OUT.fgbio.covd

    coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED \
        > $OUT.noolap.fgbio.cov
    coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED \
        -d > $OUT.noolap.fgbio.covd

    # pre MID cov covd files with noolap BED
    coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED \
        > $OUT.noolap.preMID.cov
    coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED \
        -d > $OUT.noolap.preMID.covd

    ############# BQSR fgbio bam      ##########
    #Create target interval for Indelrealigner
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx32g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T RealignerTargetCreator -R $REF \
        -known $gatk_data/GRCh37.75/1000G_phase1.indels.b37.vcf \
        -known $gatk_data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
        -I $OUT.fgbio.bam -L $BED -o ${OUT}.forIndelRealigner.intervals

    #GATK Indel Realignment
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g -Xms500m  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T IndelRealigner -R $REF -I $OUT.fgbio.bam \
        -known $gatk_data/GRCh37.75/1000G_phase1.indels.b37.vcf \
        -known $gatk_data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
        --targetIntervals ${OUT}.forIndelRealigner.intervals \
        -L $BED -o ${OUT}.realigned.bam

    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T BaseRecalibrator -R $REF -I $OUT.fgbio.bam \
        --knownSites $gatk_data/GRCh37.75/dbsnp_138.b37.vcf -nct 12 \
        --knownSites $gatk_data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
        --knownSites $gatk_data/GRCh37.75/hapmap_3.3.b37.vcf \
        -L $BED -o ${OUT}.recal_data.table

    #generate Recalibrated bam
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T PrintReads -R $REF -I ${OUT}.realigned.bam \
        -BQSR ${OUT}.recal_data.table -o ${OUT}.fgbio.bqsrCal.bam

    #########################################
    #### BQSR pre_MID bam####################
    #Create target interval for Indelrealigner
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T RealignerTargetCreator -R $REF \
        -known $gatk_data/GRCh37.75/1000G_phase1.indels.b37.vcf \
        -known $gatk_data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
        -I $OUT.preMID.bam -L $BED -o ${OUT}.forIndelRealigner.intervals

    ##GATK Indel Realignment
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T IndelRealigner -R $REF -I $OUT.preMID.bam \
        -known $gatk_data/GRCh37.75/1000G_phase1.indels.b37.vcf \
        -known $gatk_data/GRCh37.75/Mills_and_1000G_gold_standard.indels.b37.vcf \
        --targetIntervals ${OUT}.forIndelRealigner.intervals \
        -L $BED -o ${OUT}.preMID.realigned.bam

    ## GATK BQSR
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T BaseRecalibrator -R $REF -I $OUT.preMID.bam \
        --knownSites $gatk_data/GRCh37.75/dbsnp_138.b37.vcf -nct 12 \
        --knownSites $gatk_data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
        --knownSites $gatk_data/GRCh37.75/hapmap_3.3.b37.vcf \
        -L $BED -o ${OUT}.recal_data.table

    ##generate Recalibrated bam
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T PrintReads -R $REF -I ${OUT}.preMID.realigned.bam \
        -BQSR ${OUT}.recal_data.table -o ${OUT}.preMID.bqsrCal.bam

    ################################################################
    ###############################################################

    #low freq variant calls

    lofreq call-parallel --pp-threads 16 --call-indels -f $REF ${OUT}.fgbio.bqsrCal.bam -l $BED \
        -o $OUT.fgbio.bqsr.lf.vcf
    lofreq call-parallel --pp-threads 16 --call-indels -f $REF ${OUT}.preMID.bqsrCal.bam -l $BED \
        -o $OUT.preMID.bqsr.lf.vcf

    lofreq call-parallel --pp-threads 16 --call-indels -f $REF ${OUT}.fgbio.bam -l $BED \
        -o $OUT.fgbio.NObqsr.lf.vcf

    ########      VARIANT CALLING      ##################################
    #lofreq vcf merge overlapping indels
    /usr/local/bin/	 $OUT.fgbio.bqsr.lf.vcf \
        > $OUT.fgbio.bqsr.lf.mergeIndels.vcf
    /usr/local/bin/lofreq2_indel_ovlp.py $OUT.fgbio.NObqsr.lf.vcf \
        > $OUT.fgbio.NObqsr.lf.mergeIndels.vcf

    #### Vardict
    vardict -G $REF -f 0.001 -N $OUT -b $OUT.fgbio.bam -c 1 -S 2 -E 3 -g 4 $BED8 \
     | teststrandbias.R | var2vcf_valid.pl -N $OUT -E -f 0.01 > ${OUT}.fgbio.vardict.vcf

    vardict -G $REF -f 0.001 -N $OUT -b $OUT.fgbio.bqsrCal.bam -c 1 -S 2 -E 3 -g 4 $BED8 \
     | teststrandbias.R | var2vcf_valid.pl -N $OUT -E -f 0.01 > ${OUT}.fgbio.bqsr.vardict.vcf

    #germline var calling GATK HC
    /usr/java/jdk1.8.0_181-amd64/bin/java -jar -Xmx16g  /usr/opt/gatk-3.8/GenomeAnalysisTK-3.8.jar -T HaplotypeCaller -I $OUT.fgbio.bqsrCal.bam -R $REF -L $BED -o $OUT.fgbio.bqsr.gatkHC.vcf
    #$GATK -T HaplotypeCaller -I $OUT.preMID.bqsrCal.bam -R $REF -L $BED -o $OUT.preMID.bqsr.gatkHC.vcf

    #java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.preMID.filter.vcf > $OUT.preMID.filter.cosmic.vcf
    SnpSift annotate $cosmic $OUT.fgbio.bqsr.lf.mergeIndels.vcf \
        > $OUT.fgbio.lf.cosmic.vcf
    SnpSift annotate $cosmic $OUT.fgbio.bqsr.gatkHC.vcf \
        > $OUT.fgbio.gatkHC.cosmic.vcf

    #java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.preMID.filter.cosmic.vcf > $OUT.preMID.filter.cosmic.clinvar.vcf
    SnpSift annotate $clinvar $OUT.fgbio.lf.cosmic.vcf \
        > $OUT.fgbio.lf.cosmic.clinvar.vcf
    SnpSift annotate $clinvar $OUT.fgbio.gatkHC.cosmic.vcf \
        > $OUT.fgbio.gatkHC.cosmic.clinvar.vcf

    # dbNSFP annotation for polyphen and pathogenic score
    #java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.preMID.filter.cosmic.clinvar.vcf > $OUT.preMID.filter.cosmic.clinvar.dbnsfp.vcf
    #java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.lf.cosmic.clinvar.vcf \
    #    > $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf
    #java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.gatkHC.cosmic.clinvar.vcf \
    #    > $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf

    #final vars from clinvar, cosmic &dbnsfp

    #java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf \
    # CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
    #    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.lf.finalvars.txt

    #java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf \
    # CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
    #    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.gatkHC.finalvars.txt

done

######  clean up ######
rm *.MID.fq.gz
#rm *fp.vcf
rm *var.vcf
rm *cosmic.vcf
#rm *sus.bam
#rm *MID.bam
rm *.sanitise*.bam
rm *notused4consensus.bam
rm *.realigned.ba[mi]
rm *table

#############################################################
##########   summarize coverage results and reporting ########
for f in *covd
do
    awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' $f \
        > $f_covd.tmp
    awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($7>=b)n++}END{print m,(n/FNR*100.0)}' \
     OFS="\t" $f_covd.tmp $f \
        > ${f%.covd}_covMetrics.txt
done
    #awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($6>=b)n++}END{print m,b,(n/FNR*100.0)}'
    # Summarize PCR metrics
    # new PICARD rs2
for f in *targetPCRmetrics.txt
do
        awk -v n=${f%%.target*} 'NR==1{print n,$5,$11,$14*100,$22*100.0,($17+$18)/$7*100.0}' \
            OFS="\t" <(head -n8 $f | tail -n1) \
            > ${f%%.txt}_summary.txt
        f2=${f%%.target*}_covMetrics.txt
        paste ${f%%.txt}_summary.txt $f2 > ${f2%%_cov*}_combined_cov_metrics.txt

done

echo "SampleID    Total_Reads    #UQ_Reads_Aligned    %UQ_Reads_Aligned    %Bases_OnTarget_Aligned     %Bases_OnTarget_Total    Mean_Coverage    %Coverage_Uniformity" \
    > final_metrics_report.txt
cat *_combined_cov_metrics.txt >> final_metrics_report.txt

#rm *_summary.txt
#rm *_covMetrics.txt
#rm *cov_metrics.txt
######################################################################
######### MID summary outputs   #####################################
## pre and post MID coverage
printf "Chr\tStart\tEND\tAmpliconID\tLength\n" > rowheader

    #make the header columnsfile
for f in *olap.preMID.cov
do
    awk '{print $1,$2,$3,$5,$7}' OFS="\t" $f > colheader
done

    cat rowheader colheader > cov_table_header

    ## parse coverage numbers
for f in *olap.preMID.cov
do
    #g=${f%.no*}.noolap.fgbio.cov
    echo "${f%.no*}_preMID" > ${f}.preCov.tmp
    awk '{print $6}' $f >> ${f}.preCov.tmp

    echo "${f%.no*}_postMID" > ${f}.postCov.tmp
    awk '{print $6}' ${f%.no*}.noolap.fgbio.cov >> ${f}.postCov.tmp

    paste ${f}.preCov.tmp ${f}.postCov.tmp > ${f}.pre-postCov

done

paste cov_table_header *pre-postCov > pre_post-MID_coverage.txt

# MID family size distribution
# Average fam size
# Variant binning in various AF buckets pre-postMID

for f in *groupdbyMID.bam
do
    samtools view $f |grep -o "RX:Z.*." | sed 's/RX:Z://g' | sort | uniq -c | sort -k1,1n | sed -e 's/^\s*//g' -e 's/\s/\t/g' > ${f%.group*}.mid_fam_counts.txt

    cut -f1 ${f%.group*}.mid_fam_counts.txt > ${f%.group*}.counts2hist
    cat ${f%.group*}.counts2hist | sort | uniq -c |sort -k2,2n | sed -e 's/^\s*//g' -e 's/\s/\t/g' > ${f%.group*}.freq.hist.txt

    #printf "Average (overall) Family Size\n"
    echo ${f%.group*} > ${f%.group*}.avgFamSize
    awk '{sum+=$1}END{print sum/NR}' ${f%.group*}.mid_fam_counts.txt >> ${f%.group*}.avgFamSize

    #printf "Average (>2 members) Family Size\n"
    awk '{if($1>3) print}' ${f%.group*}.mid_fam_counts.txt| awk '{sum+=$1}END{print sum/NR}' > ${f%.group*}.realFamSize
    cat ${f%.group*}.avgFamSize ${f%.group*}.realFamSize > ${f%.group*}.avg_real_famSize.txt

done
#rm *tmp

  printf "Sample\nAverage (overall) Family Size\nAverage (>2 members) Family Size\n" > fam_size_header

  paste fam_size_header *famSize.txt > MID_FamilySize_Stats_report.txt

  ## variant AF binning
  #vcf=$1

for f in *.bqsr.lf.mergeIndels.vcf
do

    ### Classify variants into 6 bins based on Allele frequency pre & post MID to show noise clean
    #### due to false positive reduction

    ## total variants preMID
    echo ${f%.bqsr*} > ${f%.bqsr}.totalvars.tmp
    grep -vc ^"#" $f >> ${f%.bqsr}.totalvars.tmp

    ## total variants postMID/fgbio
    #echo ${f%.bqsr*} > ${f%.bqsr}.totalvars.tmp
    #grep -vc ^"#" $f >> ${f%.bqsr}.totalvars.tmp

    ## bin variants in AF <0.1
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if($9<0.00099) n++}END {print n}' > ${f%.bqsr}.af0_lt0.1p.tmp


    ## bin variants in AF < 0.5%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.00099) && ($9<0.0049)) n++}END {print n}' > ${f%.bqsr}.af1_0.1-.5p.tmp

    ### bin variants with AF between 0.5 and 1%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0049)&&($9<0.0099)) n++}END {print n}' > ${f%.bqsr}.af2_0.5-1p.tmp

    ### bin variants with AF between 1% and 5%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0099)&&($9<0.0499)) n++}END {print n}' > ${f%.bqsr}.af3_1-5p.tmp

    ### bin variants with AF between 5% and 10%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0499)&&($9<0.0999)) n++}END {print n}' > ${f%.bqsr}.af4_5-10p.tmp

    ### bin variants with AF between 10% and 50%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0999)&&($9<0.4999)) n++}END {print n}' > ${f%.bqsr}.af5_10-50p.tmp

    ### bin variants with AF >50%
    grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if($9>=0.4999) n++}END {print n}' > ${f%.bqsr}.af6_gt50p.tmp

    cat ${f%.bqsr}.totalvars.tmp ${f%.bqsr}.af*tmp > ${f%.bqsr}_varbins
done

#### merge from all samples for final report
printf "Sample\nTotal_Vars\nVars_AF <0.1\nVars_AF 0.1-0.5\nVars_AF 0.5-1%%\nVars_AF 1-5%%\nVars_AF 5-10%%\nVars_AF 10-50%%\nVars_AF >= 50%%\n" > bin_header

paste bin_header *_varbins > variant_AFbins_prePOST-MID.txt

#clean up
#rm *tmp
#rm *bins
#rm *header
#rm *Size
#rm *Size.txt
#rm *counts.txt
#rm *pre-postCov

#{print $1":"$2,$9}' > afTable2plot2.txt
