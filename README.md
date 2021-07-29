# xenoClassify

Xenoclassify 1.0: This Seqware workflow classifies short-read sequencing data generated from xenograft samples using [XenoClassify](https://github.com/oicr-gsi/xenoclassify).

 ![Xenoclassify, how it works](docs/xenoclassify_wf.png)


## Overview

## Dependencies

* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/1.9.tar.gz)
* [xenoclassify 1.0](https://github.com/oicr-gsi/xenoclassify/archive/1.1.tar.gz)


## Usage

### Cromwell
```
java -jar cromwell.jar run xenoclassify.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqR1`|File|fastq file for read 1


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|fastq file for read 2
`refHost`|String|"$MM10_BWA_INDEX_ROOT/mm10.fa"|The reference Host genome to align the sample with by BWA
`refGraft`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa"|The reference Graft genome to align the sample with by BWA
`rG`|String|"'@RG\\tID:TEST-RUN_XENO\\tLB:XENOTEST\\tPL:ILLUMINA\\tPU:TEST-RUN_XENO\\tSM:TEST_XENOTEST_X'"|Read group string
`bwaMemModules`|String|"bwa/0.7.17 samtools/1.9 hg19-bwa-index/0.7.17 mm10-bwa-index/0.7.17"|modules for bwaMem sub-workflow
`outputFileNamePrefix`|String|""|Output file name prefix


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateHostBam.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`generateHostBam.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`generateHostBam.indexBam_timeout`|Int|48|Hours before task timeout
`generateHostBam.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`generateHostBam.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateHostBam.bamMerge_timeout`|Int|72|Hours before task timeout
`generateHostBam.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`generateHostBam.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`generateHostBam.runBwaMem_timeout`|Int|96|Hours before task timeout
`generateHostBam.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`generateHostBam.runBwaMem_threads`|Int|8|Requested CPU threads
`generateHostBam.runBwaMem_addParam`|String?|None|Additional BWA parameters
`generateHostBam.adapterTrimming_timeout`|Int|48|Hours before task timeout
`generateHostBam.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`generateHostBam.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`generateHostBam.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`generateHostBam.slicerR2_timeout`|Int|48|Hours before task timeout
`generateHostBam.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`generateHostBam.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`generateHostBam.slicerR1_timeout`|Int|48|Hours before task timeout
`generateHostBam.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`generateHostBam.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`generateHostBam.countChunkSize_timeout`|Int|48|Hours before task timeout
`generateHostBam.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`generateHostBam.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`generateHostBam.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`generateHostBam.trimMinLength`|Int|1|minimum length of reads to keep [1]
`generateHostBam.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`generateHostBam.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`generateHostBam.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`generateGraftBam.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`generateGraftBam.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`generateGraftBam.indexBam_timeout`|Int|48|Hours before task timeout
`generateGraftBam.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`generateGraftBam.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateGraftBam.bamMerge_timeout`|Int|72|Hours before task timeout
`generateGraftBam.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`generateGraftBam.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`generateGraftBam.runBwaMem_timeout`|Int|96|Hours before task timeout
`generateGraftBam.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`generateGraftBam.runBwaMem_threads`|Int|8|Requested CPU threads
`generateGraftBam.runBwaMem_addParam`|String?|None|Additional BWA parameters
`generateGraftBam.adapterTrimming_timeout`|Int|48|Hours before task timeout
`generateGraftBam.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBam.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`generateGraftBam.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`generateGraftBam.slicerR2_timeout`|Int|48|Hours before task timeout
`generateGraftBam.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBam.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`generateGraftBam.slicerR1_timeout`|Int|48|Hours before task timeout
`generateGraftBam.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBam.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`generateGraftBam.countChunkSize_timeout`|Int|48|Hours before task timeout
`generateGraftBam.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBam.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`generateGraftBam.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`generateGraftBam.trimMinLength`|Int|1|minimum length of reads to keep [1]
`generateGraftBam.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`generateGraftBam.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`generateGraftBam.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`sortHostBam.jobMemory`|Int|10|Memory allocated to sort task
`sortHostBam.tmpDir`|String?|None|Optionally supply tmpDir for writing chunk bam files for sorting
`sortHostBam.modules`|String|"samtools/1.9"|Names and versions of modules needed for sorting
`sortHostBam.timeout`|Int|72|Timeout for this task in hours
`sortGraftBam.jobMemory`|Int|10|Memory allocated to sort task
`sortGraftBam.tmpDir`|String?|None|Optionally supply tmpDir for writing chunk bam files for sorting
`sortGraftBam.modules`|String|"samtools/1.9"|Names and versions of modules needed for sorting
`sortGraftBam.timeout`|Int|72|Timeout for this task in hours
`classify.modules`|String|"samtools/1.9 xenoclassify/1.0"|Names and versions of modules needed for classification
`classify.jobMemory`|Int|10|Memory allocated to classify task
`classify.neitherThreshold`|Int|20|Threshold for score below which the reads are classified as 'neither'
`classify.tolerance`|Int|5|Tolerance around the mean of alignment scores for a set of reads classified as 'both'
`classify.difference`|Int|5|Difference between the sum of host and graft alignment scores for a set of reads classified as 'both'
`classify.timeout`|Int|72|Timeout for this task in hours
`filterHost.modules`|String|"samtools/1.9"|Names and versions of modules needed for filtering
`filterHost.tmpDir`|String?|None|Optionally supply tmpDir for writing chunk bam files for sorting
`filterHost.filterTags`|Array[String]|["host"]|Filter reads with these tags
`filterHost.jobMemory`|Int|5|Memory allocated to filtering task
`filterHost.timeout`|Int|72|Timeout for this task in hours


### Outputs

Output | Type | Description
---|---|---
`filteredResults`|File|bam file without host (most commonly mouse) reads
`filteredResultsIndex`|File|index file for file without host reads
`jsonReport`|File|a simple stats file with counts for differently tagged reads


## Commands
 This section lists command(s) run by WORKFLOW workflow
 
 * Running WORKFLOW
 
 Xenoclassify aligns data to host and graft genomes using imported bwaMem workflow and then classify reads depending on their alignment scores.
 
 Sort bam files by read name
 
 ```
  samtools sort -n INPUT_BAM -T TMP_DIR -o BAM_BASENAME_sorted.bam
 
 ```
 
 Classify reads with xenoclassify.py script:
 
 ```
  
  xenoclassify.py is run on name-sorted bams from host and graft
 
  python3 $XENOCLASSIFY_ROOT/bin/xenoclassify/xenoclassify.py -H HOST_BAM -G GRAFT_BAM -O . -b -p PREFIX
                                                              -n NEITHER_THRESHOLD -t TOLERANCE -d DIFFERENCE
  ...
 
  samtools view PREFIX_output.bam | awk \'{print \$NF}\' | sort | uniq -c"
  counts = os.popen(command).read().splitlines() 
 
  Following Python code looks for host/graft (+ neither, both) classification and writes a summary into .json file
 
  for line in counts:
    if line.find('CL:Z:') == 0:
          continue
    line = line.rstrip()
    line = re.sub('CL:Z:', '', line)
    tmp = line.split()
    jsonDict[tmp[1]] = tmp[0]
 
  with open(json_name, 'w') as json_file:
    json.dump(jsonDict, json_file)
  
 ```
 
 Filter reads matching the supplied classification(s) and provision filtered .bam along with its index
 
 ```
  
   The following python code splits the supplied tags and then 
   removes all matching reads from the bam file  
   
   ...
 
   inputTags =  "~{sep=' ' filterTags}"
   tags = inputTags.split()
   
   command = "samtools view -h ~{xenoClassifyBam}"
   for t in tags:
     command = command + " | grep -v \'CL:Z:" + t + "\'"
   
   command = command + " | samtools sort -O bam ~{'-T ' + tmpDir} -o ~{outputPrefix}_filtered.bam -"
   os.system(command)
   
   samtools index ~{outputPrefix}_filtered.bam ~{outputPrefix}_filtered.bai
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
