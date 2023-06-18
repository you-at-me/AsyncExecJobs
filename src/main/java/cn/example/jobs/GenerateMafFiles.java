package cn.example.jobs;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.*;

/**
 * @author Carl
 * @since 2023-05-03
 */
public class GenerateMafFiles {

    private static final String FAS_FILE_PATH = "/home/xushutan/zyshu/resource/data/fas_Danio_rerio";
    private static final String OUT_FILE_PATH = "/home/xushutan/zyshu/resource/data/";

    private final BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(3);

    private final ThreadPoolExecutor poolExecutor = new ThreadPoolExecutor(
            16,  // 线程池的核心线程数量，默认16
            24,  // 线程池的最大线程数，也就是核心线程和临时线程的总和，这里默认是8
            6, // 临时线程的存活时间数
            TimeUnit.SECONDS,  // 临时线程的存活时间单位
            blockingQueue, // 线程任务等待的阻塞队列，设置多少个阻塞的线程数量，等待执行的线程任务就只能最多是多少个，超过的线程数将被通过下述策略去处理那些超过的线程任务。
            Executors.defaultThreadFactory(),  // 创建线程的工厂类型，一般走默认就行
            new ThreadPoolExecutor.CallerRunsPolicy() // 当线程提交的任务超过了线程池的最大线程数和阻塞队列等待的线程数总和时，该是以何种策略去处理这些线程。一个可选择的策略有四种，这里设置的那些超过的线程任务替换掉等待时间最长的那个线程任务
    );

    /**
     * 注解 @SneakyThrow将避免 javac 坚持要捕获或抛出方法主体中的语句声明它们生成的任何已检查异常。
     */
    @lombok.SneakyThrows
    private int executeAsyncShell(String shell, String execTitle) {
        System.out.println(execTitle);
        System.out.printf("%s is executing ======== %s", Thread.currentThread().getName(), shell);
        ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", shell);
        builder.redirectErrorStream(true);
        Process process = builder.start();

        if (process.waitFor() != 0) {
            System.out.printf("%s executes shell command failed %n", execTitle);
            return -1;
        }
        {
            return 0;
        }
    }

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

    private String getPrefixName(String filePath, String suffix) {
        String[] split = filePath.substring(0, filePath.length() - suffix.length()).split("/");
        return split[split.length - 1];
    }

    private void generateFiles() {
        List<String> faFiles = listFilesWithEnd(FAS_FILE_PATH, ".fa");
        faFiles.forEach(filePath -> poolExecutor.submit(() -> {
            String prefix = getPrefixName(filePath, ".fa");
            String outFilePath = OUT_FILE_PATH + prefix + ".axt";
            String shell = String.format("lastz %s Cyprinus_carpio.fa O=400 E=30 K=3000 L=3000 H=2200 T=1 --format=axt --ambiguous=n --ambiguous=iupac >%s", filePath, outFilePath);
            int code = executeAsyncShell(shell, filePath);
            if (code != 0) {
                System.out.printf("%s executes shell command failed %n", filePath);
                System.exit(0);
            } else {
                System.out.printf("%s executes shell successfully %n", filePath);
            }
        }));
    }

    public static void main(String[] args) {
        new GenerateMafFiles().generateFiles();
    }

}
