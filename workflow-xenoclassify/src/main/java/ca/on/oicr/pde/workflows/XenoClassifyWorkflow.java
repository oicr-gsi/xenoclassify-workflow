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
    private String samtools;
    private String xenoClassify;

    //Memory allocation
    private Integer  bwaMem;
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
            fastqOutput = getOptionalProperty("fastq_output","");
            bamOutput = getOptionalProperty("bam_output","");

            //tools
            samtools = getProperty("samtools");
            xenoClassify = getProperty("xenoClassify");

            //java = getProperty("java");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            bwaMem = Integer.parseInt(getProperty("bwa_mem"));
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
        SqwFile file0 = this.createFile("tumor");
        file0.setSourcePath(tumorBam);
        file0.setType("application/bam");
        file0.setIsInput(true);
        if (normalBam != null) { // check for missing matched normals
            SqwFile file1 = this.createFile("normal");
            file1.setSourcePath(normalBam);
            file1.setType("application/bam");
            file1.setIsInput(true);
        }
        return this.getFiles();
    }
    
//    @Override
//    public void buildWorkflow() {
//        //Runs the script you submit to the Job
////        Job job = myScriptJob();
////
////        // You shouldn't need to modify this
////        // If your script does not submit jobs to the cluster or 'blocks', meaning
////        // that it pauses execution until completion, then you do not need the monitoring
////        // step. Remove the following line and amend the output files so that they 
////        // are defined for the job above. e.g. defineOutputFiles(job);
////        Job monitor = monitorSgeJobs(job);
//        String inputs = this.inputDir + "/" + " ";
//        String output1 = this.outputDir + "/" + " ";
//        Job parentJob = null;
//        Job t1 = test1();
//        t1.addParent(parentJob);
//        
//        Job t2 = test2(inputs);
//        t2.addParent(t1);
//        
//        Job t3 = test3();
//        t3.addParent(t1);
//        
//	try {
//            defineOutputFiles(monitor);
//	}
//	catch (Exception e) {
//	    e.printStackTrace();
//	    System.exit(1);
//	}
//    }
//
//    /**
//     * Creates a job that wraps a user-provided script. Modify this method if
//     * necessary.
//     *
//     * @return the Job that describes
////     */
////    private Job myScriptJob() {
//        
////        //Set up the job
////        Job job1 = newJob("RunScript");
////        job1.setMaxMemory(getProperty("my_script_mem_mb"));
//
//        //Create the command
////        String command = String.format("%s %s >> output", 
////                getProperty("my_script"), getProperty("my_script_parameters"));
////        job1.getCommand().addArgument(command);
////	job1.setQueue(getOptionalProperty("queue", ""));
////        return job1;
//    }
//    
//    private Job test1 (String sortedBam, String finalBam) {
//        Job test1 = new Job("test_part1");
//        Command cmd = test1.getCommand();
//        cmd.addArgument("bwa mem");
//        cmd.addArgument("-t 8");
//        cmd.addArgument("-M "+this.RefGenome);
//        cmd.addArgument(this.fastq1 + " " this.fastq2 +" |");
//        cmd.addArgument(this.samtools + " view -Sb - > $SWID.bam" + ";");
//        cmd.addArgument(this.samtool+" sort -o " +sortedBam + " -n " + finalBam);
//        test1.getMaxMemory();
//        test1.getQueue();
//        return test1;
//    }

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

    
    private Job generateMpileup(String inBam, String outPileup) {
        Job ssmpileup = getWorkflow().createBashJob("generate_mpileup");
        Command cmd = ssmpileup.getCommand();
        cmd.addArgument(this.samtools);
        cmd.addArgument("mpileup -B");
        cmd.addArgument("-f " + this.refFasta);
        cmd.addArgument(inBam);
        cmd.addArgument(">" + outPileup);
        ssmpileup.setMaxMemory(Integer.toString(this.varscanMem * 1024));
        ssmpileup.setQueue(getOptionalProperty("queue", ""));
        return ssmpileup;
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
