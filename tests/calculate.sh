#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# For bam file we do the md5sum
module load samtools/1.14
# For json file we do the md5sum
find . -name '*.json' | xargs md5sum
JUNCTIONS="./Xenoclassify.Chimeric.out.junction"
# Optional files from STAR
if [[ -f $JUNCTIONS ]];then
  find . -name '*.ReadsPerGene.out.tab' | xargs md5sum
  sort -V $JUNCTIONS | grep -v ^# | md5sum
  samtools view ./Xenoclassify.Aligned.sortedByCoord.out.bam | md5sum
else
  samtools view ./Xenoclassify_filtered.bam | md5sum
fi
