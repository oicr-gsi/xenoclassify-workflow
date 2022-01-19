#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# For bam file we do the md5sum

echo "transcriptome .bam file:"
find . -name *.bam | xargs md5sum

# For json file we do the md5sum

echo "transcriptome .bam file:"
find . -name *.json | xargs md5sum
