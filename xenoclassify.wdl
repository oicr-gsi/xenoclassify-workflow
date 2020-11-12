version 1.0

# imports workflows for the top portion of WGSPipeline
import "imports/bwaMem.wdl" as bwaMem

workflow xenoClassify {
input {
        File fastqR1
	File? fastqR2
	String refHost  = "$MM10_BWA_INDEX_ROOT/mm10.fa"
	String refGraft = "$HG19_BWA_INDEX_ROOT/hg19_random.fa"
        String rG = "'@RG\\tID:TEST-RUN_XENO\\tLB:XENOTEST\\tPL:ILLUMINA\\tPU:TEST-RUN_XENO\\tSM:TEST_XENOTEST_X'"
        String bwaMemModules = "bwa/0.7.17 samtools/1.9 hg19-bwa-index/0.7.17 mm10-bwa-index/0.7.17"
        String outputFileNamePrefix = ""
}

String outputPrefix = if outputFileNamePrefix=="" then basename(fastqR1, '.fastq.gz') else outputFileNamePrefix

call bwaMem.bwaMem as generateHostBam {
  input:
    fastqR1 = fastqR1, 
    fastqR2 = fastqR2, 
    runBwaMem_bwaRef = refHost, 
    runBwaMem_modules = bwaMemModules,
    readGroups = rG,
    outputFileNamePrefix = "host"
}

call bwaMem.bwaMem as generateGraftBam {
  input:
    fastqR1 = fastqR1,
    fastqR2 = fastqR2,
    runBwaMem_bwaRef = refGraft,
    runBwaMem_modules = bwaMemModules,
    readGroups = rG,
    outputFileNamePrefix = "graft"
}

call sortBam as sortHostBam { input: inBam = generateHostBam.bwaMemBam }
call sortBam as sortGraftBam { input: inBam = generateGraftBam.bwaMemBam }

call classify { input: hostBam = sortHostBam.sortedBam, graftBam = sortGraftBam.sortedBam, outputPrefix = outputPrefix }
call filterHost { input: xenoClassifyBam = classify.xenoClassifyBam, outputPrefix = outputPrefix }


output {
  File filteredResults = filterHost.outputBam
  File filteredResultsIndex = filterHost.outputBai
  File jsonReport = classify.jsonReport
}

parameter_meta {
  fastqR1: "fastq file for read 1"
  fastqR2: "fastq file for read 2"
  refHost: "The reference Host genome to align the sample with by BWA"
  refGraft: "The reference Graft genome to align the sample with by BWA"
  rG: "Read group string"
  bwaMemModules: "modules for bwaMem sub-workflow"
  outputFileNamePrefix: "Output file name prefix"
}

 
meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Xenoclassify 1.0"
  dependencies: [
    {
      name: "bwa/0.7.12",
      url: "https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
    },
    {
      name: "samtools/1.9",
      url: "https://github.com/samtools/samtools/archive/1.9.tar.gz"
    },
    {
      name: "xenoclassify/1.0",
      url: "https://github.com/oicr-gsi/xenoclassify/archive/1.1.tar.gz"
    }
  ]
  output_meta: {
    filteredResults: "bam file without host (most commonly mouse) reads",
    filteredResults: "index file for file without host reads",
    jsonReport: "a simple stats file with counts for differently tagged reads" 
  }
}

}


# =============================================
# run BWAmem as a subworkflow, and then
# TASK 1 of 3: sort bam
# =============================================
task sortBam {
input {
	File inBam
	Int jobMemory  = 10
        String? tmpDir
        String modules = "samtools/1.9"
        Int timeout = 72
}

command <<<
 samtools sort -n ~{inBam} ~{'-T ' + tmpDir} -o ~{basename(inBam, '.bam')}_sorted.bam
>>>

parameter_meta {
 inBam: "Input .bam file"
 tmpDir: "Optionally supply tmpDir for writing chunk bam files for sorting"
 jobMemory: "Memory allocated to sort task"
 modules: "Names and versions of modules needed for sorting"
 timeout: "Timeout for this task in hours"
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File sortedBam = "~{basename(inBam, '.bam')}_sorted.bam"
}
}

# ===================================
#  TASK 2 of 3: run xenoclassify bam
# ===================================
task classify {
input {
        File hostBam
        File graftBam
        String outputPrefix
	String modules = "samtools/1.9 xenoclassify/1.0"
	Int jobMemory = 10
        Int neitherThreshold = 20
        Int tolerance = 5
        Int difference = 5
        Int timeout = 72
}

command <<<
 set -euo pipefail
 python3 $XENOCLASSIFY_ROOT/bin/xenoclassify/xenoclassify.py -H ~{hostBam} -G ~{graftBam} -O . -b -p ~{outputPrefix} \
                                                             -n ~{neitherThreshold} -t ~{tolerance} -d ~{difference}
 python3<<CODE
 import json
 import os
 import re
 json_name = "~{outputPrefix}_tagReport.json"
 jsonDict = {}

 command = "samtools view ~{outputPrefix}_output.bam | awk \'{print \$NF}\' | sort | uniq -c"
 counts = os.popen(command).read().splitlines() 

 for line in counts:
   if line.find('CL:Z:') == 0:
         continue
   line = line.rstrip()
   line = re.sub('CL:Z:', '', line)
   tmp = line.split()
   jsonDict[tmp[1]] = tmp[0]

 with open(json_name, 'w') as json_file:
   json.dump(jsonDict, json_file)
 CODE
>>>

parameter_meta {
 hostBam:  "Input host .bam file"
 graftBam: "Input graft .bam file"
 outputPrefix: "Output prefix"
 neitherThreshold: "Threshold for score below which the reads are classified as 'neither'"
 tolerance: "Tolerance around the mean of alignment scores for a set of reads classified as 'both'"
 difference: "Difference between the sum of host and graft alignment scores for a set of reads classified as 'both'"
 jobMemory: "Memory allocated to classify task"
 timeout: "Timeout for this task in hours"
 modules: "Names and versions of modules needed for classification"
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File xenoClassifyBam  = "~{outputPrefix}_output.bam"
  File jsonReport = "~{outputPrefix}_tagReport.json"
}
}

# ================================
#  TASK 3 of 3: filter bam
# ================================
task filterHost {
input {
        File xenoClassifyBam
        String outputPrefix = "OUTPUT"
        String modules = "samtools/1.9"
        String? tmpDir
        Array[String] filterTags = ["host"]
        Int jobMemory = 5
        Int timeout = 72
}

parameter_meta {
 xenoClassifyBam: "Classified .bam file"
 tmpDir: "Optionally supply tmpDir for writing chunk bam files for sorting"
 outputPrefix: "Prefix for making filtered bam name"
 modules: "Names and versions of modules needed for filtering"
 filterTags: "Filter reads with these tags"
 jobMemory: "Memory allocated to filtering task"
 timeout: "Timeout for this task in hours"
}

command <<<
  set -euo pipefail
  python3<<CODE
  import os
  inputTags =  "~{sep=' ' filterTags}"
  tags = inputTags.split()
  
  command = "samtools view -h ~{xenoClassifyBam}"
  for t in tags:
    command = command + " | grep -v \'CL:Z:" + t + "\'"
  
  command = command + " | samtools sort -O bam ~{'-T ' + tmpDir} -o ~{outputPrefix}_filtered.bam -"
  os.system(command)
  CODE
  samtools index ~{outputPrefix}_filtered.bam ~{outputPrefix}_filtered.bai
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputBam = "${outputPrefix}_filtered.bam"
  File outputBai = "${outputPrefix}_filtered.bai"
}

}

