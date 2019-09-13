version 1.0

workflow XenoclassifyWorkflow {
input {
        File fastqR1
	File fastqR2
	String refHost
	String refGraft
        String outputPrefix
}

call generateBam as generateHostBam { input: fastqR1=fastqR1, fastqR2=fastqR2,  refGenome=refHost, prefix="host" }
call generateBam as generateGraftBam { input: fastqR1=fastqR1, fastqR2=fastqR2, refGenome=refGraft, prefix="graft" }
call sortBam as sortHostBam { input: inBam=generateHostBam.outputBam }
call sortBam as sortGraftBam { input: inBam=generateGraftBam.outputBam }

call classify { input: hostBam = sortHostBam.sortedBam, graftBam = sortGraftBam.sortedBam, outPrefix = outputPrefix }
call filterHost { input: xenoClassifyBam = classify.xenoClassifyBam, filteredBam = "~{outputPrefix}_filtered.bam" }


output {
  File filteredResults = filterHost.outputBam
}
}

# ================================
#  TASK 1 of 4: generate bam
# ================================
task generateBam {
input {
	File   fastqR1
 	File   fastqR2
        String refGenome
        String prefix
        Int?   jobMemory = 20
        String? modules = "bwa/0.7.17 samtools/0.1.19 hg19-bwa-index/0.7.17 mm10-bwa-index/0.7.17"
}

command <<<
 OUTBAM=$(echo ~{basename(fastqR1)} | sed s/_R.*/_~{prefix}.bam/)
 bwa mem -t 8 -M ~{refGenome} ~{fastqR1} ~{fastqR2} | samtools view -Sb - > $OUTBAM
 echo $OUTBAM
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputBam = read_string(stdout())
}

}

# ================================
#  TASK 2 of 4: sort bam
# ================================/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.12/hg19_random.fa
task sortBam {
input {
	File inBam
	Int? jobMemory = 10
        String? modules = "samtools/0.1.19"
}

command <<<
 samtools sort -n ~{inBam} ~{basename(inBam, '.bam')}_sorted
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File sortedBam = "~{basename(inBam, '.bam')}_sorted.bam"
}
}

# ===================================
#  TASK 3 of 4: run xenoclassify bam
# ===================================
task classify {
input {
        File hostBam
        File graftBam
        String outPrefix
	String? modules = "xenoclassify/1.0"
	Int? jobMemory = 10
        Int? neitherThreshold = 20
        Int? tolerance = 5
        Int? difference = 5
}

command <<<
 python3 $XENOCLASSIFY_ROOT/bin/xenoclassify/xenoclassify.py -H ~{hostBam} -G ~{graftBam} -O . -b -p ~{outPrefix} \
                                                             -n ~{neitherThreshold} -t ~{tolerance} -d ~{difference}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File xenoClassifyBam  = "~{outPrefix}_output.bam"
}
}

# ================================
#  TASK 4 of 4: filter bam
# ================================
task filterHost {
input {
        File xenoClassifyBam
        String filteredBam
        String? modules = "samtools/0.1.19"
        Int? jobMemory = 5

}

command <<<
  samtools view -h ~{xenoClassifyBam} | grep -v 'CL:Z:mouse' | samtools view -Sh - -b > ~{filteredBam}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputBam = "${filteredBam}"
}

}


