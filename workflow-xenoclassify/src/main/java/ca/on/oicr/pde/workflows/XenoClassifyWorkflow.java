package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class XenoClassifyWorkflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;

    // Input Data
    private String fastqR1;
    private String fastqR2;
    
    // Params BWA
    private String hostBamPrefix;
    private String graftBamPrefix;
    private String hostRefFasta;
    private String graftRefFasta;
    
    // Params Xenoclassify
    private Integer neitherThreshold;
    private Integer tolerance;
    private Integer difference;
    private String outputPrefix;
    private String fastqOutput;
    private String bamOutput;
    
    //Tools
    private String bwaMem;
    private String xenoClassify;
    private String samtools;

    //Memory allocation
    private Integer  bwaMemMem;
    private Integer xenoClassifyMem;
    private String javaMem = "-Xmx8g";

    private boolean manualOutput;
    private String queue;

    // meta-types
    private final static String FASTQ_GZ_METATYPE = "application/fastq-gz";
    private final static String BAM_METATYPE = "application/bam";

    private void init() {
        try {
            //dir
            dataDir = "data/";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            fastqR1 = getProperty("fastq_read_1");
            fastqR2 = getProperty("fastq_read_2");
            
            //genome references
            hostRefFasta = getProperty("host_ref_fasta");
            graftRefFasta = getProperty("graft_ref_fasta");
            
            //bam file names
            hostBamPrefix = getProperty("host_bam_prefix");
            graftBamPrefix = getProperty("graft_bam_prefix");
            
            //XenoClassify params
            neitherThreshold = Integer.parseInt(getOptionalProperty("neither_threshold","20"));
            tolerance = Integer.parseInt(getOptionalProperty("tolerance","5"));
            difference = Integer.parseInt(getOptionalProperty("difference","5"));
            outputPrefix = getOptionalProperty("output_prefix","");

            //tools
            bwaMem = getProperty("bwa_mem");
            xenoClassify = getProperty("xenoClassify");
            samtools = getProperty("samtools");

            //java = getProperty("java");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            bwaMemMem = Integer.parseInt(getProperty("bwa_mem"));
            xenoClassifyMem = Integer.parseInt(getProperty("xenoClassify_mem"));
            
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file_R1 = this.createFile("fastqR1");
        file_R1.setSourcePath(fastqR1);
        file_R1.setType(FASTQ_GZ_METATYPE);
        file_R1.setIsInput(true);
        
        SqwFile file_R2 = this.createFile("fastqR2");
        file_R2.setSourcePath(fastqR2);
        file_R2.setType(FASTQ_GZ_METATYPE);
        file_R2.setIsInput(true);
        
        //SqwFile for fasta files?
        
        return this.getFiles();
    }
    
    @Override
    public void buildWorkflow() {
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output/";
        String inputTumourBam = getFiles().get("tumor").getProvisionedPath();
        String tumourPileupFile = this.dataDir + this.outputFilenamePrefix + ".tumour.mpileup";
        Job tumourPileup = generateMpileup(inputTumourBam, tumourPileupFile);
        parentJob = tumourPileup;
//        if (this.normalBam != null) {
//            // general normal pileup
//            String inputNormalBam = getFiles().get("normal").getProvisionedPath();
//            String normalPileupFile = this.dataDir + this.outputFilenamePrefix + ".normal.mpileup";
//            Job normalPileup = generateMpileup(inputNormalBam, normalPileupFile);
//            normalPileup.addParent(parentJob);
//            parentJob = normalPileup;
//            
//            // SNP calls
//            Job somaticCall = runVarScanSomatic(tumourPileupFile, normalPileupFile);
//            somaticCall.addParent(parentJob);
//            
//            // Provision VCF
//            String provOutVCF = this.dataDir + this.outputFilenamePrefix + ".varscan.somatic.vcf";
//            SqwFile outVCF = createOutputFile(provOutVCF, VCF_METATYPE, this.manualOutput);
//            outVCF.getAnnotations().put("somatic VCF", "varscan");
//            somaticCall.addFile(outVCF);
//            
//            //copy number
//            Job copyNumber = runVarScanCopyNumber(tumourPileupFile, normalPileupFile);
//            copyNumber.addParent(parentJob);
//            
//            // Provision Copynumber file
//            String provCopyNumber = this.dataDir + this.outputFilenamePrefix + ".copynumber";
//            SqwFile outCopyNumber = createOutputFile(provCopyNumber, COPYNUMBER_METATYPE, this.manualOutput);
//            outCopyNumber.getAnnotations().put("copy number", "varscan");
//            copyNumber.addFile(outCopyNumber);
//        }
//        else {
        Job germlineCall = runVarScanSingleSampleMode(tumourPileupFile);
        germlineCall.addParent(parentJob);
        parentJob = germlineCall;
        String provOutVCF = this.dataDir + this.outputFilenamePrefix + ".varscan.germline.vcf";
        
        // bgzip and tabix index
        Job tabixIndexVCF = tabixIndex(provOutVCF);
        tabixIndexVCF.addParent(parentJob);
        parentJob = tabixIndexVCF;
        
        
        // Provision VCF
        String gzProvOutVCF = provOutVCF + ".gz";
        SqwFile outVCF = createOutputFile(gzProvOutVCF, VCF_GZ_METATYPE, this.manualOutput);
        outVCF.getAnnotations().put("single sample VCF", "varscan");
        parentJob.addFile(outVCF);
        
        String tbiOut = gzProvOutVCF + ".tbi"; 
        SqwFile indexOut = createOutputFile(tbiOut, TBI_METATYPE, this.manualOutput);
        indexOut.getAnnotations().put("index file for VCF", "varscan");
        parentJob.addFile(indexOut);
//        }
        
        
    }
    
    private Job generateBam(String ref_genome) {
        Job ssGenerateBam = getWorkflow().createBashJob("generate_bam");
        Command cmd = ssGenerateBam.getCommand();
        cmd.addArgument(this.bwaMem);
        cmd.addArgument("-t 8");
        cmd.addArgument("-M "+ ref_genome);
        cmd.addArgument(this.fastqR1 + " " + this.fastqR1 + " |");
        cmd.addArgument(this.samtools + " view -Sb - > $SWID.bam" + ";");
        ssGenerateBam.getMaxMemory();
        ssGenerateBam.getQueue();
        return ssGenerateBam;
    }
    
    private Job sortBam(String prefix) {
        String sortedBam = prefix + ".bam";
        String finalBam = prefix + "_sorted.bam";
        Job ssSortBam = getWorkflow().createBashJob("generate_bam");
        Command cmd = ssSortBam.getCommand();
        cmd.addArgument(this.samtools +" sort -o " + sortedBam + " -n " + finalBam);
        ssSortBam.getMaxMemory();
        ssSortBam.getQueue();
        return ssSortBam;
    }
    
    private Job runVarScanSingleSampleMode(String tumourPileupFile) {
        Job ssSNP = getWorkflow().createBashJob("varscan_germline");
        Command cmd = ssSNP.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar " + this.varscan);
        cmd.addArgument("mpileup2cns");
        cmd.addArgument(tumourPileupFile);
        cmd.addArgument("--outout-vcf 1");
        cmd.addArgument("--variants 1");
        cmd.addArgument("--min-var-freq " + this.minVarFreq.toString());
        cmd.addArgument("--min-coverage " + this.minCovTumour.toString());
        cmd.addArgument("> " + this.dataDir + this.outputFilenamePrefix+".varscan.germline.vcf");
        ssSNP.setMaxMemory(Integer.toString(varscanMem * 1024));
        ssSNP.setQueue(this.queue);
        return ssSNP;
    }

//    private Job runVarScanSomatic(String tumourPileupFile, String normalPileupFile){
//        Job somaticSNP = getWorkflow().createBashJob("varscan_Somatic");
//        Command cmd = somaticSNP.getCommand();
//        cmd.addArgument(this.java);
//        cmd.addArgument(this.javaMem);
//        cmd.addArgument("-jar "+this.varscan);
//        cmd.addArgument("somatic");
//        cmd.addArgument(tumourPileupFile);
//        cmd.addArgument(normalPileupFile);
//        cmd.addArgument("--output-vcf 1");
//        cmd.addArgument("--min-var-freq " + this.minVarFreq.toString());
//        cmd.addArgument("--min-coverage-tumour " + this.minCovTumour.toString());
//        cmd.addArgument("--min-coverage-normal " + this.minCovNormal.toString());
//        cmd.addArgument("> " + this.dataDir + this.outputFilenamePrefix+".varscan.somatic.vcf");
//        somaticSNP.setMaxMemory(Integer.toString(varscanMem * 1024));
//        somaticSNP.setQueue(this.queue);
//        return somaticSNP;
//    }

//    private Job runVarScanCopyNumber(String tumourPileupFile, String normalPileupFile) {
//        Job vsCNA = getWorkflow().createBashJob("varscan_CNA");
//        Command cmd = vsCNA.getCommand();
////        cmd.addArgument("cd " + this.outDir + ";");
//        cmd.addArgument(this.java);
//        cmd.addArgument(this.javaMem);
//        cmd.addArgument("-jar " + this.varscan);
//        cmd.addArgument("copynumber");
//        cmd.addArgument(tumourPileupFile);
//        cmd.addArgument(normalPileupFile);
//        cmd.addArgument(this.dataDir + this.outputFilenamePrefix);
//        cmd.addArgument("--min-base-qual "+this.minBaseQual);
//        cmd.addArgument("--min-map-qual "+this.minMapQ);
//        cmd.addArgument("--min-coverage "+this.minCov);
//        vsCNA.setMaxMemory(Integer.toString(varscanMem * 1024));
//        vsCNA.setQueue(this.queue);
//        return vsCNA;
//    }
    
    // tabix index
    private Job tabixIndex(String vcfFile){
        Job tabixIndexVCF = getWorkflow().createBashJob("tabix_index_vcf");
        Command cmd = tabixIndexVCF.getCommand();
        cmd.addArgument(this.bgzip);
        cmd.addArgument(vcfFile + ";");
        cmd.addArgument(this.tabix + " -p vcf "+ vcfFile+".gz");
        tabixIndexVCF.setMaxMemory(Integer.toString(varscanMem * 1024));
        tabixIndexVCF.setQueue(this.queue);
        return tabixIndexVCF;
    }
    
//    private Job vcfPostProcessSingleSample(String inVCFgz){
//        Job vcfPostProcess = getWorkflow().createBashJob("post_process_VCF");
//        Command cmd = vcfPostProcess.getCommand();
//        cmd.addArgument("export PATH=$PATH:"+this.vcftools+":"+this.bcftools);
//        cmd.addArgument(this.vcftools+"vcf-merge " + inVCFgz + " " + inVCFgz +";");
//        String tumorNameId = this.outputFilenamePrefix + ".tumor";
//        String normalNameId = this.outputFilenamePrefix + ".blood";
//        cmd.addArgument("cat "+ tumorNameId + " " + normalNameId + ">" + this.outputFilenamePrefix+".reheader.txt"+";");
//        cmd.addArgument(this.bcftools+"/bcftools");
//        cmd.addArgument("reheader");
//        cmd.addArgument("--samples");
//        cmd.addArgument(this.outputFilenamePrefix+".reheader.txt");
//        cmd.addArgument(this.outputFilenamePrefix+".temp.vcf");
//        vcfPostProcess.setMaxMemory(Integer.toString(varscanMem * 1024));
//        vcfPostProcess.setQueue(this.queue);
//        return vcfPostProcess;
//    }
}
