# xenoClassify

Xenoclassify 1.2: This Seqware workflow classifies short-read sequencing data generated from xenograft samples using [XenoClassify](https://github.com/oicr-gsi/xenoclassify).

 ![Xenoclassify, how it works](docs/xenoclassify_wf.png)


## Overview

## Dependencies

* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [star 2.7.6a](https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz)
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
`fastqR2`|File|fastq file for read 2


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`refHost`|String|"$MM10_BWA_INDEX_ROOT/mm10.fa"|The reference Host genome to align the sample with by either STAR or BWA
`refGraft`|String|"$HG19_BWA_INDEX_ROOT/hg19_random.fa"|The reference Graft genome to align the sample with by either STAR or BWA
`libraryDesign`|String|"WG"|Supported library design acronym. We support WG, EX, TS, WT and MR. Default is WG
`rG`|String|"'@RG\\tID:TEST-RUN_XENO\\tLB:XENOTEST\\tPL:ILLUMINA\\tPU:TEST-RUN_XENO\\tSM:TEST_XENOTEST_X'"|Read group string
`alignerModules`|String|"bwa/0.7.17 samtools/1.9 hg19-bwa-index/0.7.17 mm10-bwa-index/0.7.17"|modules for the aligner sub-workflow
`outputFileNamePrefix`|String|""|Output file name prefix


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`generateHostBamWG.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`generateHostBamWG.indexBam_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`generateHostBamWG.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateHostBamWG.bamMerge_timeout`|Int|72|Hours before task timeout
`generateHostBamWG.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`generateHostBamWG.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`generateHostBamWG.runBwaMem_timeout`|Int|96|Hours before task timeout
`generateHostBamWG.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`generateHostBamWG.runBwaMem_threads`|Int|8|Requested CPU threads
`generateHostBamWG.runBwaMem_addParam`|String?|None|Additional BWA parameters
`generateHostBamWG.adapterTrimming_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`generateHostBamWG.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`generateHostBamWG.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`generateHostBamWG.slicerR2_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`generateHostBamWG.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`generateHostBamWG.slicerR1_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`generateHostBamWG.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`generateHostBamWG.countChunkSize_timeout`|Int|48|Hours before task timeout
`generateHostBamWG.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`generateHostBamWG.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`generateHostBamWG.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`generateHostBamWG.trimMinLength`|Int|1|minimum length of reads to keep [1]
`generateHostBamWG.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`generateHostBamWG.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`generateHostBamWG.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`generateGraftBamWG.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`generateGraftBamWG.indexBam_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`generateGraftBamWG.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateGraftBamWG.bamMerge_timeout`|Int|72|Hours before task timeout
`generateGraftBamWG.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`generateGraftBamWG.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`generateGraftBamWG.runBwaMem_timeout`|Int|96|Hours before task timeout
`generateGraftBamWG.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`generateGraftBamWG.runBwaMem_threads`|Int|8|Requested CPU threads
`generateGraftBamWG.runBwaMem_addParam`|String?|None|Additional BWA parameters
`generateGraftBamWG.adapterTrimming_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBamWG.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`generateGraftBamWG.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`generateGraftBamWG.slicerR2_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBamWG.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`generateGraftBamWG.slicerR1_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBamWG.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`generateGraftBamWG.countChunkSize_timeout`|Int|48|Hours before task timeout
`generateGraftBamWG.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`generateGraftBamWG.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`generateGraftBamWG.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`generateGraftBamWG.trimMinLength`|Int|1|minimum length of reads to keep [1]
`generateGraftBamWG.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`generateGraftBamWG.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`generateGraftBamWG.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`generateHostBamWT.indexBam_timeout`|Int|48|hours before task timeout
`generateHostBamWT.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`generateHostBamWT.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateHostBamWT.runStar_timeout`|Int|72|hours before task timeout
`generateHostBamWT.runStar_jobMemory`|Int|64|Memory allocated for this job
`generateHostBamWT.runStar_threads`|Int|6|Requested CPU threads
`generateHostBamWT.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`generateHostBamWT.runStar_peOvNbasesMin`|Int|12|minimum number of overlap bases to trigger mates merging and realignment
`generateHostBamWT.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`generateHostBamWT.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`generateHostBamWT.runStar_chimMulmapNmax`|Int|20|maximum number of chimeric multi-alignments
`generateHostBamWT.runStar_chimScoJunNonGTAG`|Int|-4|penalty for a non-GTAG chimeric junction
`generateHostBamWT.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`generateHostBamWT.runStar_alignIntMax`|Int|100000|maximum intron size
`generateHostBamWT.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`generateHostBamWT.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`generateHostBamWT.runStar_chimJunOvMin`|Int|12|minimum overhang for a chimeric junction
`generateHostBamWT.runStar_chimSegmin`|Int|12|minimum length of chimeric segment length
`generateHostBamWT.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`generateHostBamWT.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`generateHostBamWT.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`generateHostBamWT.runStar_addParam`|String?|None|Additional STAR parameters
`generateHostBamWT.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`generateHostBamWT.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`generateHostBamWT.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`generateHostBamWT.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
`generateGraftBamWT.indexBam_timeout`|Int|48|hours before task timeout
`generateGraftBamWT.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`generateGraftBamWT.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateGraftBamWT.runStar_timeout`|Int|72|hours before task timeout
`generateGraftBamWT.runStar_jobMemory`|Int|64|Memory allocated for this job
`generateGraftBamWT.runStar_threads`|Int|6|Requested CPU threads
`generateGraftBamWT.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`generateGraftBamWT.runStar_peOvNbasesMin`|Int|12|minimum number of overlap bases to trigger mates merging and realignment
`generateGraftBamWT.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`generateGraftBamWT.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`generateGraftBamWT.runStar_chimMulmapNmax`|Int|20|maximum number of chimeric multi-alignments
`generateGraftBamWT.runStar_chimScoJunNonGTAG`|Int|-4|penalty for a non-GTAG chimeric junction
`generateGraftBamWT.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`generateGraftBamWT.runStar_alignIntMax`|Int|100000|maximum intron size
`generateGraftBamWT.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`generateGraftBamWT.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`generateGraftBamWT.runStar_chimJunOvMin`|Int|12|minimum overhang for a chimeric junction
`generateGraftBamWT.runStar_chimSegmin`|Int|12|minimum length of chimeric segment length
`generateGraftBamWT.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`generateGraftBamWT.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`generateGraftBamWT.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`generateGraftBamWT.runStar_addParam`|String?|None|Additional STAR parameters
`generateGraftBamWT.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`generateGraftBamWT.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`generateGraftBamWT.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`generateGraftBamWT.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
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
`makeFastq.jobMemory`|Int|24|Memory allocated to the task.
`makeFastq.overhead`|Int|6|Ovrerhead for calculating heap memory, difference between total and Java-allocated memory
`makeFastq.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`makeFastq.picardParams`|String|"VALIDATION_STRINGENCY=LENIENT"|Additional parameters for picard SamToFastq, Default is VALIDATION_STRINGENCY=LENIENT
`makeFastq.modules`|String|"samtools/1.9 picard/2.21.2"|Names and versions of required modules.
`generateFinalBamWT.indexBam_timeout`|Int|48|hours before task timeout
`generateFinalBamWT.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`generateFinalBamWT.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`generateFinalBamWT.runStar_timeout`|Int|72|hours before task timeout
`generateFinalBamWT.runStar_jobMemory`|Int|64|Memory allocated for this job
`generateFinalBamWT.runStar_threads`|Int|6|Requested CPU threads
`generateFinalBamWT.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`generateFinalBamWT.runStar_peOvNbasesMin`|Int|12|minimum number of overlap bases to trigger mates merging and realignment
`generateFinalBamWT.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`generateFinalBamWT.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`generateFinalBamWT.runStar_chimMulmapNmax`|Int|20|maximum number of chimeric multi-alignments
`generateFinalBamWT.runStar_chimScoJunNonGTAG`|Int|-4|penalty for a non-GTAG chimeric junction
`generateFinalBamWT.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`generateFinalBamWT.runStar_alignIntMax`|Int|100000|maximum intron size
`generateFinalBamWT.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`generateFinalBamWT.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`generateFinalBamWT.runStar_chimJunOvMin`|Int|12|minimum overhang for a chimeric junction
`generateFinalBamWT.runStar_chimSegmin`|Int|12|minimum length of chimeric segment length
`generateFinalBamWT.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`generateFinalBamWT.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`generateFinalBamWT.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`generateFinalBamWT.runStar_addParam`|String?|None|Additional STAR parameters
`generateFinalBamWT.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`generateFinalBamWT.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`generateFinalBamWT.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`generateFinalBamWT.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file


### Outputs

Output | Type | Description
---|---|---
`filteredResults`|File|bam file without host (most commonly mouse) reads
`filteredResultsIndex`|File|index file for file without host reads
`starChimeric`|File?|Chimeric Graft junctions, provisioned for WT data only
`transcriptomeBam`|File?|transcriptomeBam is a file produced for Graft WT data only
`geneReadFile`|File?|.tab file with Graft gene read outs, only for WT data
`jsonReport`|File|a simple stats file with counts for differently tagged reads


## Commands

 This section lists command(s) run by Xenoclassify workflow
 
 * Running Xenoclassify workflow
 
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
   
   command = "samtools view -h XENOCLASSIFY_BAM"
   for t in tags:
     command = command + " | grep -v \'CL:Z:" + t + "\'"
   
   command = command + " | samtools sort -O bam -T  TMP_DIR -o OUTPUT_PREFIX_filtered.bam -"
   os.system(command)
   
   samtools index OUTPUT_PREFIX_filtered.bam OUTPUT_PREFIX_filtered.bai
 
 ```
 
 Extract reads from filtered file into fastq format with picard:
 
 ```
  set -euo pipefail
  unset _JAVA_OPTIONS
  java -Xmx32G -jar picard.jar SamToFastq I=FILTERED_BAM F=FILTERED_1.fastq F2=FILTERED_2.fastq ADDITIONAL_PARAMETERS
  gzip -c FILTERED_1.fastq > OUTPUT_PREFIX_part_1.fastq.gz
  gzip -c FILTERED_2.fastq > OUTPUT_PREFIX_part_2.fastq.gz
 
 ```
 
 In addition, we run second pass STAR alignments with reads from the filtered bam extracted into fastq

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
