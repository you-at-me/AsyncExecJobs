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

    // private static final Logger logger = LoggerFactory.getLogger(AsyncExecJobs.class);SimpleConsoleFormatter
    // private static final Logger logger = Logger.getLogger("cn.example.jobs.AsyncExecJobs");
    private final Logger logger = SimpleConsoleFormatter.installFormatter(Logger.getLogger(AsyncExecJobs.class.getSimpleName()));

    private final String SAM_FILE_PATH = "/sam/";
    private final String BAM_FILE_PATH = "/bam/";
    private final String SORTED_BAM_FILE_PATH = "/sorted_bam/";
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
    private final Integer SAMPLE_THRESHOLD = 100;

    /**
     * 最终聚合的gvcf 文件
     */
    private final String ALL_COMBINED_GVCF = "all_combined.gvcf";

    /**
     * 目录集合列表
     */
    private final List<String> list = Arrays.asList(SAM_FILE_PATH, BAM_FILE_PATH, SORTED_BAM_FILE_PATH, DUPLICATED_FILE_PATH, GVCF_FILE_PATH, PLINK_FILE_PATH, GCTA_FILE_PATH, LOG_PATH, PCA_PATH);

    private final ThreadPoolExecutor poolExecutor = new ThreadPoolExecutor(
            THREAD_NUM,  // 线程池的核心线程数量，默认4
            8,  // 线程池的最大线程数，也就是核心线程和临时线程的总和，这里默认是8
            60, // 临时线程的存活时间数
            TimeUnit.SECONDS,  // 临时线程的存活时间单位
            new ArrayBlockingQueue<>(JOB), // 线程任务等待的阻塞队列，设置多少个阻塞的线程数量，等待执行的线程任务就只能最多是多少个，超过的线程数将被通过下述策略去处理那些超过的线程任务。
            Executors.defaultThreadFactory(),  // 创建线程的工厂类型，一般走默认就行
            new ThreadPoolExecutor.CallerRunsPolicy() // 当线程提交的任务超过了线程池的最大线程数和阻塞队列等待的线程数总和时，该是以何种策略去处理这些线程。一个可选择的策略有四种，这里设置的那些超过的线程任务替换掉等待时间最长的那个线程任务
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

        AsyncJobsCommandParser(String... args) {
            this.options = PARSER.parse(args);
            this.help = new CommandOption<>("--help", this.options);
            this.sourceFilePath = new CommandOption<>("--sp", this.options);
            this.fastsFilePath = new CommandOption<>("--fasta", this.options);
            this.stage = new CommandOption<>("--stage", this.options);
            this.threadNum = new CommandOption<>("--threadnum", this.options);
            this.job = new CommandOption<>("--job", this.options);
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
            group.register(LONG.VALUE, "-window-pi", "--wp1").defaultTo(50000).setDescription("the window size of Pi/Fst/xp-clr, default=5000");
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

        public Builder() {
        }

        @SuppressWarnings("all")
        private Builder(String sourceFilePath, String fastaFilePath, Integer threadNum, String stage, Integer job) {
            SOURCE_FILE_PATH = sourceFilePath;
            FASTA_FILE_PATH = fastaFilePath;
            THREAD_NUM = threadNum;
            Builder.STAGE = stage;
            Builder.JOB = job;
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

        @SuppressWarnings("all")
        private AsyncExecJobs.Builder setJob(Integer job) {
            JOB = Objects.isNull(job) ? JOB : job;
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

            Builder.instance().setSourceFilePath(parser.sourceFilePath.value.getCanonicalPath()).setFastaFilePath(parser.fastsFilePath.value.getCanonicalPath()).setThreadNum(parser.threadNum.value).setStage(parser.stage.value).setJob(parser.job.value);

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
                boolean mkdir = file.mkdir();
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
            // 程序正常退出
            System.exit(0);
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

    private void createFaIndex() {
        String faLogPath = SOURCE_FILE_PATH + LOG_PATH + "bwa_index.log";
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
            String targetPrefix = filePath.substring(filePath.length() - 15, filePath.length() - 8);
            String samFilePath = SOURCE_FILE_PATH + SAM_FILE_PATH + targetPrefix + ".sam";
            String firstFqFilePath = SOURCE_FILE_PATH + "/" + targetPrefix + ".1.fq.gz";
            String secondFqFilePath = SOURCE_FILE_PATH + "/" + targetPrefix + ".2.fq.gz";
            String bwaLogPath = SOURCE_FILE_PATH + LOG_PATH + targetPrefix + ".log";
            String shell = String.format("bwa mem -t %s -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\" %s %s %s -o %s >%s 2>&1 %n", THREAD_NUM, firstFqFilePath, secondFqFilePath, FASTA_FILE_PATH, firstFqFilePath, secondFqFilePath, samFilePath, bwaLogPath);
            CompletableFuture.supplyAsync(() -> executeAsyncShell(shell, "==== start transfer " + targetPrefix + " to sam ===="), poolExecutor).whenComplete((v, e) -> {
                if (Objects.isNull(e) && v == 0) {
                    System.out.printf("==== the %s file has finished transferring to sam ! ==== %n%n", targetPrefix);
                    latch.countDown();
                }
            }).exceptionally(throwable -> {
                System.out.printf("the %s file failed to transfer to sam %n", targetPrefix);
                System.out.printf("exception information:%s %n", throwable.getMessage());
                return null;
            });
        });
        if (latch.getCount() != 0) { // 注意线程减少计数的count if判断的写法
            try {
                latch.await();
            } catch (InterruptedException e) {
                e.printStackTrace();
                System.out.printf("latch await exception:%s %n", e.getMessage());
            }
        }
        {
            System.out.println(("========= transfer all fq.gz to sam finished! ========="));
            this.validateSam();
        }
    }

    private void validateSam() {
        System.out.println("====== start validate sam files =======");
        List<String> invalidSamFiles = new ArrayList<>();
        String shell = "gatk ValidateSamFile --INPUT %s %n";
        List<String> samFileLists = listFilesWithEnd(SOURCE_FILE_PATH + SAM_FILE_PATH, ".sam");
        CyclicBarrier cyclicBarrier = new CyclicBarrier(samFileLists.size(), () -> {
            if (invalidSamFiles.size() > 0) {
                System.out.printf("sam file validate failed, check:%s %n", invalidSamFiles);
                System.exit(0);
            }
            {
                System.out.println("====== sam files validate finished ======");
                this.samToBam();
            }
        });
        samFileLists.forEach(filePath -> CompletableFuture.supplyAsync(() -> {
            String[] split = filePath.split("/");
            return executeAsyncShell(String.format(shell, filePath), "==== start validate " + split[split.length - 1] + " ====");
        }, poolExecutor).whenComplete((v, e) -> {
            if (!Objects.isNull(e) || v != 0) {
                invalidSamFiles.add(filePath);
                try {
                    cyclicBarrier.await();
                } catch (InterruptedException | BrokenBarrierException ex) {
                    ex.printStackTrace();
                }
            }
        }));
    }

    private void samToBam() {
        System.out.println("====== start transfer sam to bam  =======");
        String originShell = "gatk SortSam -I %s -O %s --TMP_DIR %s -SO coordinate";
        List<String> TransferSamToBamFiles = new ArrayList<>();
        AtomicInteger samFileCount = new AtomicInteger();
        List<String> samFileLists = listFilesWithEnd(SOURCE_FILE_PATH + SAM_FILE_PATH, ".sam");
        CountDownLatch latch = new CountDownLatch(samFileLists.size());
        samFileLists.forEach(filePath -> {
            if (!filePath.endsWith(".sam")) {
                System.exit(0);
            }
            samFileCount.getAndIncrement();
            String[] split = filePath.substring(0, filePath.length() - ".sam".length()).split("/");
            String samFileNamePrefix = split[split.length - 1];
            String bamFilePath = SOURCE_FILE_PATH + BAM_FILE_PATH + samFileNamePrefix + ".sort.bam";
            String tmpFilePath = SOURCE_FILE_PATH + BAM_FILE_PATH + samFileNamePrefix + ".tmp";
            String shell = String.format(originShell, filePath, bamFilePath, tmpFilePath);
            CompletableFuture.supplyAsync(() -> {
                return executeAsyncShell(shell, "==== start transfer " + filePath.split("/")[1] + " to bam ====");
            }).whenComplete((v, e) -> {
                if (!Objects.isNull(e) || v != 0) {
                    TransferSamToBamFiles.add(filePath);
                }
                latch.countDown();
            });
        });
        if (latch.getCount() != 0) {
            try {
                latch.await();
            } catch (InterruptedException e) {
                e.printStackTrace();
                System.out.printf("latch await exception:%s %n", e.getMessage());
            }
        }
        {
            if (TransferSamToBamFiles.size() != 0) {
                System.out.printf("sorted to bam files failed, list:%s", TransferSamToBamFiles);
                System.exit(0);
            }
            // generate file failed
            List<String> sortBamLists = listFilesWithEnd(SOURCE_FILE_PATH + BAM_FILE_PATH, "sort.bam");
            // Verify the number of generated files
            if (sortBamLists.size() != samFileCount.get()) {
                System.out.println("======= No. of sam not equals No. of bam =======");
                System.exit(0);
            }
            System.out.println("====== sam to bam transfer finished ======");
        }

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
        if (process.waitFor() == 0) {
            return 0;
        }
        {
            // throw new IOException();
            System.out.printf("'%s' executes shell command failed %n", execTitle);
            return -1;
        }
    }

    private void first() {
        init();
        createFaIndex();
    }

    public void test(String[] args) {
        checkCommandParser(args);
        System.out.println(AsyncJobsCommandParser.getOptionsByParseArgs(args));
        logger.info("===== main run begin: =====");

        if (Objects.equals(STAGE, "ALL")) {
            first();
        } else {
            logger.info("error choice");
        }
    }


}
