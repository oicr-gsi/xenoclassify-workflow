#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

qrsh -l h_vmem=8G -cwd -now n "\
module load samtools/1.9 2>/dev/null; \
find . -regex '.*\.bam$' \
-exec sh -c \" samtools flagstat {} | tr '\n' '\t'; echo \" \; \
| sort | uniq | tr '\t' '\n'"