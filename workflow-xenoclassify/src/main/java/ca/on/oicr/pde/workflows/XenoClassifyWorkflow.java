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
    private String python;

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
            outDir = getProperty("out_dir");

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
   
        return this.getFiles();
    }
    
    @Override
    public void buildWorkflow() {
        
//        String fastqR1 = getFiles().get("fastqR1").getProvisionedPath();
//        String fastqR2 = getFiles().get("fastqR2").getProvisionedPath();
        
        Job parentJob = null;
        
        //align fastqs to human and mouse reference genome
        String hostBam = this.dataDir + this.hostBamPrefix + ".bam";
        String graftBam = this.dataDir + this.graftBamPrefix + ".bam";
        Job generateHostBam = generateBam(this.hostRefFasta, hostBam);
        Job generateGraftBam = generateBam(this.graftRefFasta, graftBam);
        generateHostBam.addParent(parentJob);
        generateGraftBam.addParent(parentJob);
        
        String sortedHostBam = this.dataDir + this.hostBamPrefix + "_sorted.bam";
        String sortedGraftBam = this.dataDir + this.graftBamPrefix + "_sorted.bam";
        Job sortHost = sortBam(hostBam, sortedHostBam);
        Job sortGraft = sortBam(graftBam, sortedGraftBam);
        sortHost.addParent(generateHostBam);
        sortGraft.addParent(generateGraftBam);
        
        Job classifyXeno = classify(sortedHostBam, sortedGraftBam);
        classifyXeno.addParent(sortHost);
        classifyXeno.addParent(sortGraft);
        
        //Provision out output bam created by XenoClassify
        String xenoClassifyBam = this.dataDir + this.outputPrefix + ".bam";
        SqwFile xenoClassifyBamFile = createOutputFile(xenoClassifyBam, BAM_METATYPE, this.manualOutput);
        xenoClassifyBamFile.getAnnotations().put("output BAM", "xenoclassify");
        classifyXeno.addFile(xenoClassifyBamFile);
        
        String filteredBam = this.dataDir + this.outputPrefix + "_filtered.bam";
        Job filter = filterHost(xenoClassifyBam, filteredBam);
        
        //Provision out filtered bam file
        SqwFile filteredBamFile = createOutputFile(filteredBam, BAM_METATYPE, this.manualOutput);
        filter.addFile(filteredBamFile);
    }
    
    private Job generateBam(String ref_genome, String outBam) {
        Job ssGenerateBam = getWorkflow().createBashJob("generate_bam");
        Command cmd = ssGenerateBam.getCommand();
        cmd.addArgument(this.bwaMem);
        cmd.addArgument("-t 8");
        cmd.addArgument("-M "+ ref_genome);
        cmd.addArgument(this.fastqR1 + " " + this.fastqR1 + " |");
        cmd.addArgument(this.samtools + " view -Sb - > " + outBam + ";");
        ssGenerateBam.getMaxMemory();
        ssGenerateBam.getQueue();
        return ssGenerateBam;
    }
    
    private Job sortBam(String inBam, String outBam) {
        Job ssSortBam = getWorkflow().createBashJob("generate_bam");
        Command cmd = ssSortBam.getCommand();
        cmd.addArgument(this.samtools + " sort -o " + outBam + " -n " + inBam);
        ssSortBam.getMaxMemory();
        ssSortBam.getQueue();
        return ssSortBam;
    }
    
    private Job classify(String hostBam, String graftBam) {
        Job ssClassify = getWorkflow().createBashJob("classify");
        Command cmd = ssClassify.getCommand();
        cmd.addArgument(this.python + " " + this.xenoClassify + " -M ");
        cmd.addArgument(hostBam + " -H " + graftBam + " -O " + this.dataDir);
        cmd.addArgument(" -b " + " -p " + this.outputPrefix);
        ssClassify.getMaxMemory();
        ssClassify.getQueue();
        return ssClassify;
    }
    
    // job for removing mouse reads from bam files
    private Job filterHost(String xenoClassifyBam, String filteredBam) {
        Job ssFilter = getWorkflow().createBashJob("filter");
        Command cmd = ssFilter.getCommand();
        cmd.addArgument("grep -v 'CL:Z:mouse' " + xenoClassifyBam + " > " + filteredBam); 
        return ssFilter;
    }
    
}
