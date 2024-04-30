#!/bin/bash
featureCounts -p -T 10 -a genes/chrX.gtf -o featurescount_out/all_featureCountx.txt hisat_out/ERR*.bam


