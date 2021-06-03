#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load samtools/0.1.19 2>/dev/null
find . -regex '.*\.bam$' | tail -n 1 | xargs samtools flagstat  | tr '\t' '\n'; echo " " | sort | uniq | tr '\t' '\n'
