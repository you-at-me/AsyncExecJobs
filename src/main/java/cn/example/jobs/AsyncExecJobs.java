package cn.example.jobs;

import edu.sysu.pmglab.commandParser.CommandGroup;
import edu.sysu.pmglab.commandParser.CommandOption;
import edu.sysu.pmglab.commandParser.CommandOptions;
import edu.sysu.pmglab.commandParser.CommandParser;
import edu.sysu.pmglab.commandParser.exception.ParameterException;
import edu.sysu.pmglab.commandParser.types.*;
import edu.sysu.pmglab.commandParser.usage.DefaultStyleUsage;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.*;

import static cn.example.jobs.AsyncExecJobs.Builder.*;
import static edu.sysu.pmglab.commandParser.CommandItem.HELP;
import static edu.sysu.pmglab.commandParser.CommandItem.HIDDEN;

/**
 * created by Carl on 2023/03/28
 */
public class AsyncExecJobs {

    // private static final com.sun.org.slf4j.internal.Logger log = LoggerFactory.getLogger(AsyncExecJobs.class);
    // private static final Logger logger = Logger.getLogger("cn.example.jobs.AsyncExecJobs");

    private final Logger logger = SimpleConsoleFormatter.installFormatter(Logger.getLogger(AsyncExecJobs.class.getSimpleName()));

    private final String SAM_FILE_PATH = "/sam/";
    private final String BAM_FILE_PATH = "/bam/";
    private final String DUPLICATED_FILE_PATH = "/duplicated/";
    private final String GVCF_FILE_PATH = "/gvcf/";
    private final String PLINK_FILE_PATH = "/plink/";
    private final String GCTA_FILE_PATH = "/gcta/";
    private final String LOG_PATH = "/log/";
    private final String PCA_PATH = "/pca/";

    private final String eigvenval_file_name = "all_raw_snp.gcta_pca.eigenval";
    private final String eigvenc_file_name = "all_snp.plink.eigenvec";
    private final Integer bam_sort_choice = 0;

    /**
     * 样本数量阈值，超过该值为大量样本
     */
    private final short SAMPLE_THRESHOLD = 100;

    /**
     * 最终聚合的gvcf 文件
     */
    private final String ALL_COMBINED_GVCF = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_combined.gvcf";
    private final String ALL_RAW_CVF_PATH = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_raw.vcf";
    private final String ALL_RAW_SNP_VCF_PATH = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_raw_snp.vcf";
    private final String ALL_RAW_SNP_HFA_PATH = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_raw_snp_hardfiltereannotated.vcf";
    private final String ALL_RAW_SNP_HF_PATH = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_raw_snp_hardfiltered.vcf";

    /**
     * 目录集合列表
     */
    private final List<String> list = Arrays.asList(SAM_FILE_PATH, BAM_FILE_PATH, DUPLICATED_FILE_PATH, GVCF_FILE_PATH, PLINK_FILE_PATH, GCTA_FILE_PATH, LOG_PATH, PCA_PATH);

    private final BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(JOB);

    private final ThreadPoolExecutor poolExecutor = new ThreadPoolExecutor(
            THREAD_NUM,  // 线程池的核心线程数量，默认4
            8,  // 线程池的最大线程数，也就是核心线程和临时线程的总和，这里默认是8
            6, // 临时线程的存活时间数
            TimeUnit.SECONDS,  // 临时线程的存活时间单位
            blockingQueue, // 线程任务等待的阻塞队列，设置多少个阻塞的线程数量，等待执行的线程任务就只能最多是多少个，超过的线程数将被通过下述策略去处理那些超过的线程任务。
            // Executors.defaultThreadFactory(),  // 创建线程的工厂类型，一般走默认就行，如果希望改名字也可自定义
            r -> {
                SecurityManager s = System.getSecurityManager();
                ThreadGroup group = (s != null) ? s.getThreadGroup() : Thread.currentThread().getThreadGroup();
                String namePrefix = "prefix-thread-pool" + new AtomicInteger(1).getAndIncrement() + "-thread-";
                Thread t = new Thread(group, r, namePrefix +  new AtomicInteger(1).getAndIncrement(), 0);
                if (t.isDaemon()) t.setDaemon(false);
                if (t.getPriority() != Thread.NORM_PRIORITY) t.setPriority(Thread.NORM_PRIORITY);
                return t;
            },
            new ThreadPoolExecutor.DiscardOldestPolicy() // 当线程提交的任务超过了线程池的最大线程数和阻塞队列等待的线程数总和时，该是以何种策略去处理这些线程。一共可选择的策略有四种，这里设置的是把这个被拒绝的任务替换掉在阻塞队列中等待时间最长的那个线程任务，然后由线程池重试execute执行这个任务，如果这个时候核心线程和临时线程仍在忙碌，则就走线程池的那一套执行流程，也就是放入到阻塞队列当中。一般不太建议使用的是：CallerRunsPolicy执行策略，因为这种策略在对待被拒绝的任务是，直接在调用execute方法的线程中运行(run)这个被拒绝的任务，而一般情况下我们都是在主线程在通过线程池执行这个任务的，所以被拒绝的任务最终会被主线程执行，如果这个任务执行很慢，会严重影响到整个程序的执行性能。
    );

    {
        this.poolExecutor.allowCoreThreadTimeOut(true);
    }

    private static class AsyncJobsCommandParser {
        private static final CommandParser PARSER = new CommandParser(false);

        private final CommandOptions options;
        private final CommandOption<?> help;
        private final CommandOption<File> sourceFilePath;
        private final CommandOption<File> fastsFilePath;
        private final CommandOption<String> stage;
        private final CommandOption<Integer> threadNum;
        private final CommandOption<Integer> job;
        private final CommandOption<Float> mm;
        private final CommandOption<Float> mf;
        private final CommandOption<Long> s1;
        private final CommandOption<Long> wp1;
        private final CommandOption<Integer> wps1;
        private final CommandOption<Integer> n1;

        AsyncJobsCommandParser(String... args) {
            this.options = PARSER.parse(args);
            this.help = new CommandOption<>("--help", this.options);
            this.sourceFilePath = new CommandOption<>("--sp", this.options);
            this.fastsFilePath = new CommandOption<>("--fasta", this.options);
            this.stage = new CommandOption<>("--stage", this.options);
            this.threadNum = new CommandOption<>("--threadnum", this.options);
            this.job = new CommandOption<>("--job", this.options);
            this.mm = new CommandOption<>("--mm", this.options);
            this.mf = new CommandOption<>("--mf", this.options);
            this.s1 = new CommandOption<>("--s1", this.options);
            this.wp1 = new CommandOption<>("--wp1", this.options);
            this.wps1 = new CommandOption<>("--wps1", this.options);
            this.n1 = new CommandOption<>("--n1", this.options);
        }

        private static AsyncJobsCommandParser parse(String... args) {
            return new AsyncJobsCommandParser(args);
        }

        private AsyncJobsCommandParser parseFile(File argsFile) throws IOException {
            return new AsyncJobsCommandParser(CommandParser.readFromFile((edu.sysu.pmglab.container.File) argsFile));
        }

        /**
         * 将每一对参数以换行的形式展示，而每一对参数命令与值用空格隔开，例如:
         * <p>-SP ../aa \ <p>
         * -fasta-file-path ab
         */
        private static CommandOptions getOptionsByParseArgs(String... args) {
            return PARSER.parse(args);
        }

        /**
         * Get CommandParser, toString() is the same to usage()
         */
        private static CommandParser getParser() {
            return PARSER;
        }

        /**
         * Get the usage of CommandParser
         */
        private static String usage() {
            return PARSER.toString();
        }

        /**
         * Get CommandOptions
         */
        @Deprecated
        private CommandOptions getOptions() {
            return this.options;
        }

        static {
            // 设置程序名 (默认值: <main class>)。
            PARSER.setProgramName("java -jar asyncExecJobs.jar");
            // 例如解析：bgzip compress <file> --level 5 -t 4 -o ~/test.gz, 当 offset = 3 时，传入的下列参数将跳过前 3 个参数，解析 "--level 5 -t 4 -o ~/test.gz"
            PARSER.offset(0);
            // 是否为调试模式 (默认值: false)。
            PARSER.debug(false);
            // 是否将 @ 识别为取地址符 (地址对应的文件内容作为参数传入) (默认值: true)。
            PARSER.usingAt(true);
            // 设置最大匹配参数项个数 (默认值: -1)。传入的参数项达到最大个数时，后续的参数不再解析，而是作为最后一个匹配的参数项的值。例如：对于指令: bgzip compress <file> decompress <file>; maxMatchedItems = 1 时，下列参数只匹配 "bgzip"，剩下的参数 "compress <file> decompress <file>" 则作为 "bgzip" 参数项的值
            PARSER.setMaxMatchedNum(-1);
            // 当没有指令被传入时，是否自动添加 help 指令 (默认值: false)。
            PARSER.setAutoHelp(false);
            // 设置参数命令展示的文档格式 (默认值: DefaultStyleUsage.UNIX_TYPE_1)。
            PARSER.setUsageStyle(DefaultStyleUsage.UNIX_TYPE_1);

            CommandGroup group = PARSER.addCommandGroup("Options");
            group.register(IType.NONE, "-h", "--help", "-help").addOptions(HELP, HIDDEN);
            group.register(STRING.VALUE, "-STAGE", "--stage").defaultTo(STAGE).setDescription("choose step to run,-STAGE 1 is mapping+snp calling+vcf generate,-STAGE 2 is vcf filter,-STAGE 3 is pca+admixture+Phylogenetic analyse to divide populations,-STAGE 4 is pca+admixture+Phylogenetic analyse+LDdecay+snp denisity+Genetic diversity datas+selective sweep analysis etc.");
            group.register(INTEGER.VALUE, "-J", "--job").defaultTo(JOB).setDescription("the number of jobs to submit, default=8");
            group.register(INTEGER.VALUE, "-T", "--threadnum").defaultTo(THREAD_NUM).setDescription("the number of thread, default=4");
            // group.register((IType) FILE.validateWith(true, false, true), "-SP", "--sp").defaultTo(DEFAULT_SOURCE_FILE_PATH);
            group.register(FILE.VALUE, "-SP", "--sp").defaultTo(SOURCE_FILE_PATH);
            group.register(FILE.VALUE, "-fasta-file-path", "--fasta").defaultTo(FASTA_FILE_PATH);
            group.register(INTEGER.VALUE, "-N", "--n1").defaultTo(8).setDescription("times of bootstraps, default=8");
            group.register(LONG.VALUE, "-s", "--s1").defaultTo(100000).setDescription("the window size of snp density, default=10w");
            group.register(LONG.VALUE, "-window-pi", "--wp1").defaultTo(5000).setDescription("the window size of Pi/Fst/xp-clr, default=5000");
            group.register(INTEGER.VALUE, "-window-pi-step", "--wps1").defaultTo(2000).setDescription("the window step of Pi/Fst/xp-clr, default=2000");
            group.register(FLOAT.VALUE, "-MF", "--mf").defaultTo(0.05).setDescription("the Minor Allele Frequency to filter snp, default=0.05");
            // 染色体例子 2,3,8,24
            group.register(FLOAT.VALUE, "-MM", "--mm").defaultTo(0.2).setDescription("the Max-missing rate, default=0.2");
            group.register(IType.NONE, "-CHR", "--chr").setDescription("Chromosomes splited with ',' e.g -CHR 1,2,3,4,5");
            group.register(IType.NONE, "-K", "--k").setDescription("your belief of the number of ancestral populations");
        }
    }

    static class Builder {
        static String SOURCE_FILE_PATH = "/home/user/module/resource/LargeYellowCroaker_RAD/clean/test";
        static String FASTA_FILE_PATH = SOURCE_FILE_PATH + "/GCF_000972845.2/GCF_000972845.2_L_crocea_2.0_genomic.fna";
        static Integer THREAD_NUM = 4;
        static String STAGE = "ALL";
        static Integer JOB = 8;
        static Float MM = 0.2f;
        static Float MF = 0.05f;
        static Long S1 = 100000L;
        static Long WP1 = 5000L;
        static Integer WPS1 = 2000;

        public Builder() {
        }

        @SuppressWarnings("all")
        private Builder(String sourceFilePath, String fastaFilePath, Integer threadNum, String stage, Integer job, Float mm, Float mf, Long s1, Long wp1, Integer wps1) {
            SOURCE_FILE_PATH = sourceFilePath;
            FASTA_FILE_PATH = fastaFilePath;
            THREAD_NUM = threadNum;
            STAGE = stage;
            JOB = job;
            MM = mm;
            MF = mf;
            S1 = s1;
            WP1 = wp1;
            WPS1 = wps1;
        }

        private static AsyncExecJobs.Builder instance() {
            return new AsyncExecJobs.Builder();
        }

        private AsyncExecJobs.Builder setSourceFilePath(String sourceFilePath) {
            SOURCE_FILE_PATH = sourceFilePath;
            return this;
        }

        private AsyncExecJobs.Builder setFastaFilePath(String fastaFilePath) {
            if (!Objects.equals(fastaFilePath, "")) {
                FASTA_FILE_PATH = fastaFilePath;
            }
            return this;
        }

        private AsyncExecJobs.Builder setThreadNum(Integer threadNum) {
            THREAD_NUM = Objects.nonNull(threadNum) ? threadNum : THREAD_NUM;
            return this;
        }

        private AsyncExecJobs.Builder setStage(String stage) {
            STAGE = Objects.equals(stage, "") ? STAGE : stage;
            return this;
        }

        private AsyncExecJobs.Builder setJob(Integer job) {
            JOB = Objects.isNull(job) ? JOB : job;
            return this;
        }

        private AsyncExecJobs.Builder setMm(Float mm) {
            MM = Objects.isNull(mm) ? MM : mm;
            return this;
        }

        private AsyncExecJobs.Builder setMf(Float mf) {
            MF = Objects.isNull(mf) ? MF : mf;
            return this;
        }

        @SuppressWarnings("all")
        private AsyncExecJobs.Builder setWps1(Integer wps1) {
            WPS1 = Objects.isNull(wps1) ? WPS1 : wps1;
            return this;
        }
    }

    /**
     * 为{@link java.util.logging.Logger}实现自定义的日志输出，可以输出IDEA自动识别源码位置的日志格式，方便调试
     */
    private static class SimpleConsoleFormatter extends Formatter {

        @Override
        public String format(LogRecord record) {
            String message = formatMessage(record);
            String throwable = "";
            if (Objects.nonNull(record.getThrown())) {
                StringWriter sw = new StringWriter();
                PrintWriter pw = new PrintWriter(sw);
                pw.println();
                record.getThrown().printStackTrace(pw);
                pw.close();
                throwable = "\n" + sw;
            }
            String dateTime = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(record.getMillis());
            StackTraceElement stackTrace = Thread.currentThread().getStackTrace()[8];
            return String.format("%s [%s] (%s:%d) %s%s\n", dateTime, Thread.currentThread().getName(), stackTrace.getFileName(), stackTrace.getLineNumber(), message, throwable);
        }

        /**
         * 将{@link SimpleConsoleFormatter}实例指定为{@link Logger}的输出格式
         */
        public static Logger installFormatter(Logger logger) {
            if (Objects.nonNull(logger)) {
                // 禁用原输出handler，否则会输出两次
                logger.setUseParentHandlers(false);
                ConsoleHandler consoleHandler = new ConsoleHandler();
                consoleHandler.setFormatter(new SimpleConsoleFormatter());
                logger.addHandler(consoleHandler);
            }
            return logger;
        }
    }

    /**
     * 为{@link java.util.logging.Logger}自定义输出日志格式信息到指定的文件当中，这种方法仅限于将格式的日志信息打印到文本，却不能打印到控制台，控制台使用的依旧是原始的日志格式信息
     */
    @Deprecated
    private void formatLogToFile() {
        FileHandler fileHandler = null;
        try {
            fileHandler = new FileHandler("", true);
        } catch (IOException e) {
            e.printStackTrace();
        }
        assert fileHandler != null;
        fileHandler.setFormatter(new Formatter() {
            private final SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

            public String format(LogRecord record) {
                StringBuilder sb = new StringBuilder();
                String dataFormat = this.sdf.format(record.getMillis());
                sb.append(dataFormat).append(" ");
                sb.append("level:").append(record.getLevel()).append(" ");
                sb.append(record.getMessage()).append("\n");

                return sb.toString();
            }
        });
        logger.addHandler(fileHandler);
    }

    private void checkCommandParser(String[] args) {
        try {
            if (args.length == 0) {
                System.out.println(AsyncJobsCommandParser.usage());
            }
            AsyncJobsCommandParser parser = AsyncJobsCommandParser.parse(args);
            if (parser.help.isPassedIn) {
                System.out.println(AsyncJobsCommandParser.usage());
                AsyncJobsCommandParser.getParser().iterator().forEachRemaining(System.out::println);
                System.exit(0);
            }

            // test
            if (parser.sourceFilePath.isPassedIn) {
                // System.out.println(parser.sourceFilePath.value); // 输出对于值的对象
                // System.out.println("=============");
                // System.out.println(parser.sourceFilePath.value.getName()); // aa
                // System.out.println(parser.sourceFilePath.value.getPath()); // ..\..\aa
                // System.out.println(parser.sourceFilePath.value.getParent()); // ..\..
                // // G:\Workplaces\Java\IdeaExercise\aa
                // System.out.println(parser.sourceFilePath.value.getCanonicalPath());
                System.out.println(parser.sourceFilePath.value.getCanonicalFile());
                // // G:\Workplaces\Java\IdeaExercise\Works\jdk8-juc-thread\..\..\aa
                // System.out.println(parser.sourceFilePath.value.getAbsoluteFile());
                // System.out.println("=============");
                // // 对象地址参数: edu.sysu.pmglab.commandParser.CommandOption@6d5380c2
                // System.out.println(parser.sourceFilePath);
                // System.out.println(AsyncJobsCommandParser.getOptionsByParseArgs(args));
            }

            Builder.instance().setSourceFilePath(parser.sourceFilePath.value.getCanonicalPath()).setFastaFilePath(parser.fastsFilePath.value.getCanonicalPath()).setThreadNum(parser.threadNum.value).setStage(parser.stage.value).setJob(parser.job.value).setMm(parser.mm.value).setMf(parser.mf.value).setWps1(parser.wps1.value);

        } catch (Exception | Error e) {
            e.printStackTrace();
            logger.log(Level.WARNING, String.format("%s", e.getMessage()));
            // logger.warning(String.format("%s", e.getMessage()));
            if (e instanceof ParameterException) {
                logger.info(AsyncJobsCommandParser.usage());
            }
            System.exit(0);
        }
    }

    private void init() {
        if (!new File(SOURCE_FILE_PATH).exists()) {
            System.out.printf("the %s is not exist", SOURCE_FILE_PATH);
            System.exit(0);
        }
        List<String> fileDirPath = new ArrayList<>(list.size());
        list.forEach(path -> fileDirPath.add(SOURCE_FILE_PATH + path));
        for (String fileDir : fileDirPath) {
            File file = new File(fileDir);
            if (!file.exists()) {
                file.mkdir();
                System.out.printf("=====the %s is not exist, already created %n", fileDir);
            }
        }
        checkFile(SOURCE_FILE_PATH);
    }

    /**
     * 检查文件是否符合fqgz命名规则，若有不符合的文件直接退出整个程序的运行
     */
    private void checkFile(String sourceFilePath) {
        System.out.println("========= start check *.fq.gz file =========");
        List<String> invalidFqgzFiles = new ArrayList<>();
        Objects.requireNonNull(listFilesWithEnd(sourceFilePath, ".gz")).forEach(filePath -> {
            if (!filePath.endsWith("1.fq.gz") && !filePath.endsWith("2.fq.gz")) {
                invalidFqgzFiles.add(filePath);
            }
        });
        if (invalidFqgzFiles.size() > 0) {
            System.out.printf("file pattern check failed, check:%s %n", invalidFqgzFiles);
            System.exit(0);
        }
    }

    private void createFaIndex() {
        String faLogPath = SOURCE_FILE_PATH + LOG_PATH + "fa_index.log";
        String shell = String.format("bwa index -a bwtsw %s >%s 2>&1 %n", FASTA_FILE_PATH, faLogPath);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "==== start create index ===="), poolExecutor).whenComplete((v, e) -> {
            if (Objects.isNull(e) && v == 0) {
                System.out.printf("generate index finished, check:%s %n%n", faLogPath);
                this.transferFqgzFilesToSam();
            }
        }).exceptionally(throwable -> {
            System.out.printf("generate index failed, check:%s %n", faLogPath);
            System.out.printf("exception information:%s %n", throwable.getMessage());
            return null;
        });
    }

    private void transferFqgzFilesToSam() {
        System.out.println("========= start transfer *.fq.gz file to sam =========");
        List<String> gzList = listFilesWithEnd(SOURCE_FILE_PATH, ".1.fq.gz");
        if (gzList.size() < 1) {
            System.out.println(".1.fq.gz files not found!");
            System.exit(0);
        }
        CountDownLatch latch = new CountDownLatch(gzList.size());
        gzList.forEach(filePath -> {
            String targetPrefix = getPrefixName(filePath, ".1.fq.gz");
            String samFilePath = SOURCE_FILE_PATH + SAM_FILE_PATH + targetPrefix + ".sam";
            String firstFqFilePath = SOURCE_FILE_PATH + "/" + targetPrefix + ".1.fq.gz";
            String secondFqFilePath = SOURCE_FILE_PATH + "/" + targetPrefix + ".2.fq.gz";
            String fgToSamLogPath = SOURCE_FILE_PATH + LOG_PATH + targetPrefix + ".fqToSam.log";
            String shell = String.format("bwa mem -t %s -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\" %s %s %s -o %s >%s 2>&1 %n", THREAD_NUM, firstFqFilePath, secondFqFilePath, FASTA_FILE_PATH, firstFqFilePath, secondFqFilePath, samFilePath, fgToSamLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "==== start transfer " + targetPrefix + " to sam ===="), poolExecutor).whenComplete((v, e) -> {
                if (Objects.isNull(e) && v == 0) {
                    System.out.printf("==== the %s file has finished transferring to sam ! ==== %n%n", targetPrefix);
                    latch.countDown();
                } else {
                    System.out.printf("==== the %s file failed transferring to sam", targetPrefix);
                    System.out.printf("==== the %s file has finished transferring to sam ! ==== %n%n", targetPrefix);
                }
            }).exceptionally(throwable -> {
                System.out.printf("the %s file failed to transfer to sam %n", targetPrefix);
                System.out.printf("exception information:%s %n", throwable.getMessage());
                return null;
            });
        });
        latchWait(latch);
        System.out.println(("========= transfer all fq.gz to sam finished! ========="));
    }

    private void validateSam() {
        System.out.println("====== start validate sam files =======");
        String shell = "gatk ValidateSamFile --INPUT %s >%s 2>&1 %n";
        List<String> samFileLists = listFilesWithEnd(SOURCE_FILE_PATH + SAM_FILE_PATH, ".sam");
        CountDownLatch latch = new CountDownLatch(samFileLists.size());
        for (String filePath : samFileLists) {
            String[] split = filePath.split("/");
            String prefixName = split[split.length - 1];
            String validateSamLogPath = SOURCE_FILE_PATH + LOG_PATH + getPrefixName(filePath, ".sam") + ".validateSam.log";
            CompletableFuture.supplyAsync(() -> executeAsyncShell(String.format(shell, filePath, validateSamLogPath), "==== start validate " + prefixName + " ===="), poolExecutor).whenComplete((v, e) -> {
                if (v != 0) {
                    System.out.printf("==== %s validate failed ==== %n", prefixName);
                    System.exit(0);
                }
                System.out.printf("==== finish validate %s ==== %n", prefixName);
                latch.countDown();
            });
        }
        latchWait(latch);
        System.out.println("====== all sam files validate finished ======");
    }

    private void samToBam() {
        System.out.println("====== start transfer sam to bam  =======");
        String originShell = "gatk SortSam -I %s -O %s --TMP_DIR %s -SO coordinate >%s 2>&1 %n";
        AtomicInteger samFileCount = new AtomicInteger();
        List<String> samFileLists = listFilesWithEnd(SOURCE_FILE_PATH + SAM_FILE_PATH, ".sam");
        CountDownLatch latch = new CountDownLatch(samFileLists.size());
        samFileLists.forEach(filePath -> {
            if (!filePath.endsWith(".sam")) {
                System.exit(0);
            }
            //Using atomic class cas to solve the counting problem under multi-threading
            samFileCount.getAndIncrement();
            String samPrefixName = getPrefixName(filePath, ".sam");
            String bamFilePath = SOURCE_FILE_PATH + BAM_FILE_PATH + samPrefixName + ".sort.bam";
            String tmpFilePath = SOURCE_FILE_PATH + BAM_FILE_PATH + samPrefixName + ".tmp";
            String samToBamLogPath = SOURCE_FILE_PATH + LOG_PATH + samPrefixName + ".samToBam.log";
            String shell = String.format(originShell, filePath, bamFilePath, tmpFilePath, samToBamLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, String.format("==== start transfer %s to bam ====", samPrefixName)), poolExecutor).whenComplete((v, e) -> {
                if (!Objects.isNull(e) || v != 0) {
                    System.out.printf("==== the %s 'sam to bam' failed ==== %n", samPrefixName);
                    System.exit(0);
                }
                System.out.printf("==== the %s 'sam to bam' has transferred ==== %n", samPrefixName);
                latch.countDown();
            });
        });
        latchWait(latch);
        List<String> sortBamLists = listFilesWithEnd(SOURCE_FILE_PATH + BAM_FILE_PATH, "sort.bam");
        // Verify the number of generated files
        if (sortBamLists.size() != samFileCount.get()) {
            System.out.println("======= No. of sam not equals No. of bam =======");
            System.exit(0);
        }
        System.out.println("====== all 'sam to bam' have been converted ======");
    }

    private void markDuplicates() {
        System.out.println("======= start mark duplicates =======");
        List<String> bamFileLists = listFilesWithEnd(SOURCE_FILE_PATH + BAM_FILE_PATH, ".sort.bam");
        String originShell = "gatk MarkDuplicates -I %s -O %s.deduplicated -M %s.deduplicated.metrics >%s 2>&1 %n";
        CountDownLatch latch = new CountDownLatch(bamFileLists.size());
        bamFileLists.forEach(filePath -> {
            String bamPrefixName = getPrefixName(filePath, ".sort.bam");
            String bamPrefixNamePath = SOURCE_FILE_PATH + DUPLICATED_FILE_PATH + bamPrefixName;
            String markDuplicatesLogPath = SOURCE_FILE_PATH + LOG_PATH + bamPrefixName + ".duplicate.log";
            String Shell = String.format(originShell, filePath, bamPrefixNamePath, bamPrefixNamePath, markDuplicatesLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(Shell, String.format("==== start markDuplicates %s ====", bamPrefixName)), poolExecutor).whenComplete((v, e) -> {
                if (v != 0 || !Objects.isNull(e)) {
                    System.out.printf("%s failed markDuplicates %n", bamPrefixName);
                    System.exit(0);
                }
                System.out.printf("%s finished markDuplicates %n", bamPrefixName);
                latch.countDown();
            });
        });
        latchWait(latch);
        System.out.println("====== all .sort.bam files finished markDuplicates ======");
    }

    private void createBaiIndex() {
        System.out.println("======= start create deduplicated and fa index =======");
        List<String> deduplicatedLists = listFilesWithEnd(SOURCE_FILE_PATH + DUPLICATED_FILE_PATH, ".deduplicated");
        String deduplicatedIndex = "samtools index %s %n";
        CountDownLatch latch = new CountDownLatch(deduplicatedLists.size() + 2);
        deduplicatedLists.forEach(filePath -> {
            String filePrefixName = getPrefixName(filePath, "");
            CompletableFuture.supplyAsync(() -> executeAsyncShell(String.format(deduplicatedIndex, filePath), String.format("==== start create %s index ====", filePrefixName)), poolExecutor).whenComplete((v, e) -> {
                if (v != 0 || !Objects.isNull(e)) {
                    System.out.printf("failed %s index %n", filePrefixName);
                    System.exit(0);
                }
                System.out.printf("finished %s index %n", filePrefixName);
                latch.countDown();
            });
        });
        createDicAndFaiIndex(latch);
        System.out.println("====== all deduplicated index and fa index finished ! ======");
    }

    private void createDicAndFaiIndex(CountDownLatch latch) {
        /* 以下是针对fa文件进行.dic和.fai索引文件的构建，速度很快 */
        String faidxIndex = String.format("samtools faidx %s %n", FASTA_FILE_PATH);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(faidxIndex, "==== start create faidx index ===="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed create faidx index");
                System.exit(0);
            }
            System.out.println("finished create faidx index");
            latch.countDown();
        });

        String createSequenceDictionary = SOURCE_FILE_PATH + LOG_PATH + "CreateSequenceDictionary.log";
        String dicShell = String.format("gatk CreateSequenceDictionary -R %s >%s 2>&1 %n", FASTA_FILE_PATH, createSequenceDictionary);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(dicShell, "==== start create sequence dictionary index ===="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed create sequence dictionary index");
                System.exit(0);
            }
            System.out.println("finished create sequence dictionary index");
            latch.countDown();
        });
        latchWait(latch);
    }

    private void generateGvcfFiles() {
        System.out.println("======= start generate gvcf files =======");
        List<String> deduplicatedLists = listFilesWithEnd(SOURCE_FILE_PATH + DUPLICATED_FILE_PATH, ".deduplicated");
        String shell = "gatk HaplotypeCaller --pcr-indel-model CONSERVATIVE -ERC GVCF -R %s -I %s -O %s >%s 2>&1 %n";
        CountDownLatch latch = new CountDownLatch(deduplicatedLists.size());
        deduplicatedLists.forEach(filePath -> {
            String prefixName = getPrefixName(filePath, ".deduplicated");
            String gvcfFilepath = SOURCE_FILE_PATH + GVCF_FILE_PATH + prefixName + ".gvcf";
            String gvcfLogPath = SOURCE_FILE_PATH + LOG_PATH + prefixName + ".gvcf.log";
            String gvcfShell = String.format(shell, FASTA_FILE_PATH, filePath, gvcfFilepath, gvcfLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(gvcfShell, String.format("==== start generate %s gvcf file ====", prefixName)), poolExecutor).whenComplete((v, e) -> {
                if (v != 0 || !Objects.isNull(e)) {
                    System.out.printf("failed generate %s gvcf file %n", prefixName);
                    System.exit(0);
                }
                System.out.printf("finished generate %s gvcf file %n", prefixName);
                latch.countDown();
            });
        });
        latchWait(latch);
        System.out.println("======= all gvcf files have been generated ! =======");
    }

    private void mergeGvcfFiles() {
        List<String> gvcfLists = listFilesWithEnd(SOURCE_FILE_PATH + GVCF_FILE_PATH, ".gvcf");
        if (gvcfLists.size() > SAMPLE_THRESHOLD) {
            System.out.println("too many files");
            System.exit(0);
        }
        StringBuilder shellBuilder = new StringBuilder(String.format("gatk CombineGVCFs -R %s", FASTA_FILE_PATH));
        for (String filePath : gvcfLists) {
            shellBuilder.append(String.format(" -V %s", filePath));
        }
        String shell = shellBuilder.toString();
        shell = String.format(shell + " -O " + ALL_COMBINED_GVCF + " %n");
        String finalShell = shell;
        CompletableFuture.supplyAsync(() -> executeAsyncShell(finalShell, "==== start merge gvcf files ===="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed merge gvcf file");
                System.exit(0);
            }
            System.out.println("finished merge gvcf file");
            System.out.println("======= all gvcf files have been merged ! =======");
            gvcfToVcf();
        });
    }

    private void gvcfToVcf() {
        String shell = String.format("gatk GenotypeGVCFs -R %s -V %s -O %s %n", FASTA_FILE_PATH, ALL_COMBINED_GVCF, ALL_RAW_CVF_PATH);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "==== start transfer gvcf to vcf files ===="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed merge gvcf file");
                System.exit(0);
            }
            System.out.println("finished merge gvcf file");
            System.out.println("======= all gvcf to vcf files have been transferred ! =======");
            System.out.println("======= the first step over ！！！ =======");
            // stepOne();
        });
    }

    private void stepOne() {
        System.out.println("======= the second step is running ！ =======");
        String allRawSnpLogPath = SOURCE_FILE_PATH + LOG_PATH + "all_raw_snp.vcf.log";
        String shell = String.format("gatk SelectVariants -V %s -select-type SNP -O %s >%s 2>&1 %n", ALL_RAW_CVF_PATH, ALL_RAW_SNP_VCF_PATH, allRawSnpLogPath);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "===== start generate all_raw_snp.vcf ====="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed generate all_raw_snp.vcf");
                System.exit(0);
            }
            System.out.println("finished generate all_raw_snp.vcf");
            stepTwo();
        });
    }

    private void stepTwo() {
        String hardFilteredEAnnotatedLogPath = SOURCE_FILE_PATH + LOG_PATH + "all_raw_snp_hardfiltereannotated.vcf.log";
        String shell = String.format("gatk VariantFiltration -V %s -filter 'QD < 2.0' --filter-name 'QD2' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O %s >%s 2>&1 %n", ALL_RAW_SNP_VCF_PATH, ALL_RAW_SNP_HFA_PATH, hardFilteredEAnnotatedLogPath);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "===== start generate all_raw_snp_hardfiltereannotated.vcf ====="), poolExecutor).whenComplete((v, e) -> {
            if (v != 0 || !Objects.isNull(e)) {
                System.out.println("failed generate all_raw_snp_hardfiltereannotated.vcf");
                System.exit(0);
            }
            System.out.println("finished generate all_raw_snp_hardfiltereannotated.vcf");

            String hardFilteredLogPath = SOURCE_FILE_PATH + LOG_PATH + "all_raw_snp_hardfiltered.vcf.log";
            String shell2 = String.format("gatk SelectVariants --exclude-filtered true -V %s -O %s >%s 2>&1 %n", ALL_RAW_SNP_HFA_PATH, ALL_RAW_SNP_HF_PATH, hardFilteredLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(shell2, "===== start generate all_raw_snp_hardfiltered.vcf ====="), poolExecutor).whenComplete((x, y) -> {
                if (x != 0 || !Objects.isNull(y)) {
                    System.out.println("failed generate all_raw_snp_hardfiltered.vcf");
                    System.exit(0);
                }
                System.out.println("finished generate all_raw_snp_hardfiltered.vcf");
                stepThree();
            });
        });
    }

    private void stepThree() {
        String ALL_SNP_VCF_PATH = SOURCE_FILE_PATH + GVCF_FILE_PATH + "all_snp.vcf";
        String softFilterLogPath = SOURCE_FILE_PATH + LOG_PATH + "all_snp.softfilter.log";
        String shell = String.format("vcftools --vcf %s --max-missing %s --maf %s --mac 3 --recode --recode-INFO-all --out %s >%s 2>&1 %n", ALL_RAW_SNP_HF_PATH, MM, MF, ALL_SNP_VCF_PATH, softFilterLogPath);
        CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "===== start generate all_snp.vcf ====="), poolExecutor).whenComplete((x, y) -> {
            if (x != 0 || !Objects.isNull(y)) {
                System.out.println("failed generate all_snp.vcf");
                System.exit(0);
            }
            System.out.println("finished generate all_snp.vcf");
            System.out.println("======= the second step over ！！！ =======");
        });
    }


    private String getPrefixName(String filePath, String suffix) {
        String[] split = filePath.substring(0, filePath.length() - suffix.length()).split("/");
        return split[split.length - 1];
    }

    private void latchWait(CountDownLatch latch) {
        if (latch.getCount() != 0) {
            try {
                latch.await(); // Block the code after calling the await method
            } catch (InterruptedException e) {
                e.printStackTrace();
                System.out.printf("latch await exception:%s %n", e.getMessage());
                System.exit(0);
            }
        }
    }

    /**
     * 返回对应文件夹路径path下以end_with结尾的文件集合列表
     */
    private List<String> listFilesWithEnd(String path, String endWith) {
        File dirObj = new File(path);
        if (!dirObj.isDirectory()) {
            System.out.printf("listFilesWithEnd is empty, %s is not a directory %n", path);
            System.exit(0);
        }
        System.out.printf("==== start list files, path:%s, endWith : %s ==== %n", path, endWith);
        List<String> conformFileLists = new ArrayList<>();
        for (File file : Objects.requireNonNull(dirObj.listFiles())) {
            String filePath = file.getAbsolutePath();
            if (filePath.endsWith(endWith)) {
                conformFileLists.add(filePath);
            }
        }
        System.out.printf("==== listFilesWithEnd : %s %n%n", conformFileLists);
        return conformFileLists;
    }

    /**
     * 注解 @SneakyThrow将避免 javac 坚持要捕获或抛出方法主体中的语句声明它们生成的任何已检查异常。
     */
    @lombok.SneakyThrows
    private int executeAsyncShell(String shell, String execTitle) {
        String threadName = Thread.currentThread().getName();
        System.out.println(execTitle);
        System.out.printf("%s is executing ======== %s", threadName, shell);
        ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", shell);
        builder.redirectErrorStream(true);

        Process process = builder.start();
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line;
        while ((line = reader.readLine()) != null) {
            System.out.println(line);
        }
        // process.waitFor()的调用查看执行结果的返回值，当值不为0的时候表示脚本执行出现问题
        if (process.waitFor() != 0) {
            // throw new IOException();
            System.out.printf("%s executes shell command failed %n", execTitle);
            return -1;
        }
        {
            return 0;
        }
    }

    private void first() {
        init();
        createFaIndex();

        validateSam();
        samToBam();
        markDuplicates();
        createBaiIndex();
        generateGvcfFiles();
        mergeGvcfFiles();

    }

    private void second() {
        stepOne();

    }

    public void test(String[] args) {
        checkCommandParser(args);
        System.out.println(AsyncJobsCommandParser.getOptionsByParseArgs(args));
        logger.info("===== main run begin: =====");

        if (Objects.equals(STAGE, "ALL")) {
            first();
            second();
        } else {
            logger.info("error choice");
        }
    }

}
