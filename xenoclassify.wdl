version 1.0

# imports workflows for the top portion of WGSPipeline
import "imports/pull_bwamem2.wdl" as bwaMem
import "imports/pull_star.wdl" as star

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

workflow xenoClassify {
input {
  Array[InputGroup] inputs
  String reference
  String libraryDesign
  String outputFileNamePrefix = ""
}

String outputPrefix = outputFileNamePrefix
# We support only single-lane data for WG
File fastqR1 = select_first(inputs).fastqR1
File fastqR2 = select_first(inputs).fastqR2
String rG = select_first(inputs).readGroup

if (libraryDesign == "WG" || libraryDesign == "EX" || libraryDesign == "TS") {

 call bwaMem.bwamem2 as generateHostBamWG {
   input:
     fastqR1 = fastqR1, 
     fastqR2 = fastqR2, 
     runBwamem2_readGroups = rG,
     reference = "mm10",
     outputFileNamePrefix = "host"
 }

 call bwaMem.bwamem2 as generateGraftBamWG {
   input:
     fastqR1 = fastqR1,
     fastqR2 = fastqR2,
     runBwamem2_readGroups = rG,
     reference = reference,
     outputFileNamePrefix = "graft"
 }
 call sortBam as sortHostBamWG { input: inBam = select_first([generateHostBamWG.bwamem2Bam, generateHostBamWT.starBam]) }
 call sortBam as sortGraftBamWG { input: inBam = select_first([generateGraftBamWG.bwamem2Bam, generateGraftBamWT.starBam]) }
 call classify as classifyWG { input: hostBam = sortHostBamWG.sortedBam, graftBam = sortGraftBamWG.sortedBam, outputPrefix = outputFileNamePrefix }
 call filterHost as filterHostWG { input: xenoClassifyBam = classifyWG.xenoClassifyBam, outputPrefix = outputFileNamePrefix }
 call mergeReports as mergeReportsWG { input: inputReports = [classifyWG.jsonReport], inputRgs = [rG], outputPrefix = outputFileNamePrefix }
} 

if (libraryDesign == "WT" || libraryDesign == "MR") {
  scatter(inp in inputs) {
    call star.star as generateHostBamWT {
      input:
        inputGroups = [inp],
        reference = "mm10",
        outputFileNamePrefix = "host"
    }
    call sortBam as sortHostBamWT { input: inBam = generateHostBamWT.starBam }
    call star.star as generateGraftBamWT {
      input:
        inputGroups = [inp],
        reference = reference,
        outputFileNamePrefix = "graft"
    }
    call sortBam as sortGraftBamWT { input: inBam = generateGraftBamWT.starBam }
    
    call classify as classifyWT { input: hostBam = sortHostBamWT.sortedBam, graftBam = sortGraftBamWT.sortedBam, outputPrefix = outputFileNamePrefix }
    call filterHost as filterHostWT { input: xenoClassifyBam = classifyWT.xenoClassifyBam, outputPrefix = outputFileNamePrefix }
    call makeFastq { input: inputBam = filterHostWT.outputBam, rG = inp.readGroup, outputPrefix = outputFileNamePrefix }
  }
  # Also, re-align filtered data with star and merge reports from the classify task
  call star.star as generateFinalBamWT {
    input:
      inputGroups = makeFastq.fastqData,
      reference = reference,
      outputFileNamePrefix = outputFileNamePrefix
  }
  call mergeReports as mergeReportsWT { input: inputReports = classifyWT.jsonReport, inputRgs = makeFastq.readGroup, outputPrefix = outputFileNamePrefix }
}

output {
  File filteredResults = select_first([generateFinalBamWT.starBam, filterHostWG.outputBam]) 
  File filteredResultsIndex = select_first([generateFinalBamWT.starIndex, filterHostWG.outputBai])
  File? starChimeric = generateFinalBamWT.starChimeric
  File? transcriptomeBam = generateFinalBamWT.transcriptomeBam
  File? geneReadFile = generateFinalBamWT.geneReadFile
  File jsonReport = select_first([mergeReportsWT.jsonReport, mergeReportsWG.jsonReport])
}

parameter_meta {
  inputs: "Array of fastq files for read 1 and 2 along with rG string"
  libraryDesign: "Supported library design acronym. We support WG, EX, TS, WT and MR. Default is WG"
  reference: "The reference of Graft to align the data with by either STAR or BWA"
  outputFileNamePrefix: "Output file name prefix"
}

 
meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Xenoclassify 1.3: This Seqware workflow classifies short-read sequencing data generated from xenograft samples using [XenoClassify](https://github.com/oicr-gsi/xenoclassify).\n\n ![Xenoclassify, how it works](docs/xenoclassify_wf.png)\n"
  dependencies: [
    {
      name: "bwa/0.7.12",
      url: "https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
    },
    {
      name: "star/2.7.6a",
      url: "https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz"
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
    filteredResultsIndex: "index file for file without host reads",
    starChimeric: "Chimeric Graft junctions, provisioned for WT data only",
    transcriptomeBam: "transcriptomeBam is a file produced for Graft WT data only",
    geneReadFile: ".tab file with Graft gene read outs, only for WT data",
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

# ====================================
#   Optional for WT: Make fastq files
# ====================================
task makeFastq {
input {
  Int jobMemory = 24
  Int overhead = 6
  Int timeout = 20
  File inputBam
  String outputPrefix
  String rG
  String picardParams = "VALIDATION_STRINGENCY=LENIENT"
  String modules = "samtools/1.9 picard/2.21.2"
}

Int javaMemory = jobMemory - overhead

command <<<
 set -euo pipefail
 unset _JAVA_OPTIONS
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar SamToFastq I=~{inputBam} F=FILTERED_1.fastq F2=FILTERED_2.fastq ~{picardParams}
 gzip -c FILTERED_1.fastq > ~{outputPrefix}_part_1.fastq.gz
 gzip -c FILTERED_2.fastq > ~{outputPrefix}_part_2.fastq.gz

>>>

parameter_meta {
 inputBam: "Input bam file, BWA-aligned reads"
 outputPrefix: "Output prefix for the result file"
 jobMemory: "Memory allocated to the task."
 rG: "Read group string for passing to the downstream process"
 overhead: "Ovrerhead for calculating heap memory, difference between total and Java-allocated memory"
 picardParams: "Additional parameters for picard SamToFastq, Default is VALIDATION_STRINGENCY=LENIENT"
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  InputGroup fastqData = {"fastqR1": "~{outputPrefix}_part_1.fastq.gz", "fastqR2": "~{outputPrefix}_part_2.fastq.gz", "ReadGroup": "~{rG}"}
  String readGroup = "~{rG}"
}
}

# =========================================
#   Optional for WT: Merge classify reports
# =========================================
task mergeReports {
  input {
    Array[File] inputReports
    Array[String] inputRgs
    String outputPrefix
    String modules = ""
    Int jobMemory = 4
    Int timeout = 4
  }

  parameter_meta {
    inputReports: "Array of input JSON files"
    inputRgs: "Array of RG strings"
    outputPrefix: "Output prefix for the result file"
    jobMemory: "Memory for the task, in gigabytes"
    modules: "Environment modules for the task"
    timeout: "Timeout for the task, in hours"
  }

 command <<<
   python <<CODE
   import json
   import re

   r = "~{sep=' ' inputReports}"
   inputJsons = r.split()
   inputRgs = "~{sep=' ' inputRgs}"
   
   data = {}

   def jsonRead(fileName):
       with open(fileName, "r") as f:
           jsonText = f.readlines()
           jsonText = "".join(jsonText)
           jsonText = jsonText.strip()
       return json.loads(jsonText)

   matches = re.findall('(?<=[ID]:)([\S]*)', inputRgs)

   for j in range(len(inputJsons)):
       if matches[j]:
           data[matches[j]] = jsonRead(inputJsons[j])

   metrics_file = "~{outputPrefix}_tagReport.json"
   with open(metrics_file, "w") as m:
       m.write(json.dumps(data, indent=2))
 
   CODE
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File jsonReport = "~{outputPrefix}_tagReport.json"
  }
}


