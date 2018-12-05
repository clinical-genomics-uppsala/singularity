# Requires:
* GenomeAnalysisTK-3.8.jar
* jdk-8u181-linux-x64.rpm
* singularity

# Build
´´´
singularity build accel_amplicon_hs.simg accel_amplicon_hs.def accel_amplicon_hs.def
´´´

# Run
'''
singularity run  --app swift_hs -B /projects/wp1/nobackup/ngs/utveckling/analys/2018/20181012_AccAmpEGFR-HS_Test1_SWIFT -B gatk:/data/gatk  -B bedfiles:/data/bedfiles -B refgenomes/:/data/refgenomes accel_amplicon_hs.simg
'''
