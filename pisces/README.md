Illuminas variant caller.

# Build
'''
sudo singularity build pisces.simage  pisces.def
'''

# Run
'''
singularity run -B /projects --app pisces pisces.simage  -Bam test.bam  -G hg19 -i manifest_ext  -MinBQ 30 -MinVF 0.0001 -SSFilter False -MinMQ 30 -MaxVQ 100 -MinDepthFilter 500 -MinVQ 0 -VQFilter 20 -gVCF True -ReportNoCalls True -CallMNVs True -MaxMNVLength 3 -MaxGapBetweenMNV 1 -Collapse True -ReportRcCounts True -ReportTsCounts True -ThreadByChr True -outFolder test
'''
