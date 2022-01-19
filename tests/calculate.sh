#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# For bam file we do the md5sum
module load samtools/1.9
samtools view Xenoclassify_filtered.bam | md5sum

# For json file we do the md5sum
find . -name *.json | xargs md5sum
