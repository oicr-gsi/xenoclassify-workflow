#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# For bam file we do the md5sum
module load samtools/1.14 jq/1.6
# For json file we do the md5sum
for f in $(find . -name '*.json');do cat $f | jq -S | md5sum;done
JUNCTIONS="./Xenoclassify.Chimeric.out.junction"
# Optional files from STAR
if [[ -f $JUNCTIONS ]];then
  find . -name '*.ReadsPerGene.out.tab' | xargs md5sum
  sort -V $JUNCTIONS | grep -v ^# | md5sum
  samtools view ./Xenoclassify.Aligned.sortedByCoord.out.bam | md5sum
else
  samtools view ./Xenoclassify_filtered.bam | md5sum
fi
