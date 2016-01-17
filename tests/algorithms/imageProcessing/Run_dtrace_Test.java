package algorithms.imageProcessing;

import algorithms.imageProcessing.features.CornersOfHouseTest;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * run dtrace on a class to see the classes it uses.
 * 
 * @author nichole
 */
public class Run_dtrace_Test extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    protected String processID = null;

    protected String CWD = null;

    protected Process dTraceProcess = null;

    protected Thread task = null;
    
    protected boolean dtraceStopped = false;
    
    // might need superuser
    protected boolean enable = false;
    
    
    public Run_dtrace_Test(String testName) {
        super(testName);
    }
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void test1() throws Exception {
        
        if (!enable) {
            return;
        }
                
        processID = getProcessID();

        CWD = getCWD();
        
        CountDownLatch procStartLatch = new CountDownLatch(1);
        CountDownLatch dTraceStartedLatch = new CountDownLatch(1);
        CountDownLatch procEndedLatch = new CountDownLatch(1);
        
        TRunner runner = new TRunner(procStartLatch, dTraceStartedLatch, procEndedLatch);
        runner.start();
        
        log.info("wait for procStartLatch countdown");
        procStartLatch.await();
        
        
        startDTrace(dTraceStartedLatch);
        
        log.info("wait for procEndedLatch countdown");
        procEndedLatch.await(); 
        
        log.info("procEndedLatch unlatched");
                                    
        if (!dtraceStopped) {
            // sometimes, have grepped the wrong PID and so there's an error with dtrace
            stopDTrace();
        }
    }
    
    protected class TRunner extends Thread {

        protected final CountDownLatch procStartLatch;
        protected final CountDownLatch procEndedLatch;
        protected final CountDownLatch dTraceStartedLatch;
        
        public TRunner(CountDownLatch procStartLatch, 
            CountDownLatch dTraceStartedLatch, CountDownLatch procEndedLatch) {
            this.procStartLatch = procStartLatch;
            this.procEndedLatch = procEndedLatch;
            this.dTraceStartedLatch = dTraceStartedLatch;
        }
        
        @Override
        public void run() {
            
            log.info("run");
            
            try {
                Thread.sleep(5000);
                
                log.info("release procStartLatch");
                procStartLatch.countDown();
                
                log.info("wait for dTraceStartedLatch");
                dTraceStartedLatch.await();
                
                
                log.info("dTraceStartedLatch unlatched");
                
                CornersOfHouseTest test = new CornersOfHouseTest(
                    "CornersOfHouseTest");
                test.testProcess();

                stopDTrace();
                procEndedLatch.countDown();
                
            } catch(Exception e) {
                log.severe(e.getMessage());
                if (procStartLatch.getCount() > 0) {
                    procStartLatch.countDown();
                }
                if (procEndedLatch.getCount() > 0) {
                    procEndedLatch.countDown();
                }
                if (dTraceStartedLatch.getCount() > 0) {
                    dTraceStartedLatch.countDown();
                }
            }
        }
    }
    
    protected String getProcessID() throws IOException, InterruptedException {

        String pid = "0";
        String cmd2 = "/bin/ps -aef | grep ExtendedDTraceProbes | grep Run_dtrace_Test | grep -v sudo | awk -F' ' '{print $2}'";
        List<String> commands = new ArrayList<String>();
        commands.add("/bin/sh");
        commands.add("-c");
        commands.add(cmd2);
        ProcessBuilder processBuilder2 = new ProcessBuilder(commands);
        Process proc2 = processBuilder2.start();
        int exit = proc2.waitFor();
        if (exit == 0) {
            String content = getContent(proc2.getInputStream());
            if (content != null && content.length() > 0) {
                String[] piditems = content.split("\n");
                if (piditems != null && piditems.length >= 1) {
                    pid = piditems[0];
                }
            }
        }
        log.info("PID=" + pid);
        return pid;
    }
    protected String getCWD() throws IOException, InterruptedException {

        String pwd = "";
        String cmd2 = "/bin/pwd";
        List<String> commands = new ArrayList<String>();
        commands.add("/bin/sh");
        commands.add("-c");
        commands.add(cmd2);
        ProcessBuilder processBuilder2 = new ProcessBuilder(commands);
        Process proc2 = processBuilder2.start();
        int exit = proc2.waitFor();
        if (exit == 0) {
            pwd = getContent(proc2.getInputStream()).trim();
        }
        log.info("CWD=" + pwd);
        return pwd;
    }
    
    protected void startDTrace(final CountDownLatch dTraceStartedLatch) 
        throws IOException, InterruptedException {

         final String path = CWD + "/dtrace";
        
         log.info("startDTrace");
        
         task = new Thread() {

            @Override
            public void run() {

                try {
                                        
                    List<String> commands = new ArrayList<String>();
                    //commands.add(path);
                    commands.add(path + "/j_flow_obj_alloc.d");
                    commands.add("-p");
                    commands.add(processID);

                    ProcessBuilder processBuilder2 = new ProcessBuilder(commands);
                    processBuilder2.directory(new File(path));
                    processBuilder2.inheritIO();
                    dTraceProcess = processBuilder2.start();

                    // allows the process to start
                    dTraceStartedLatch.countDown();
                    
                    int exit = dTraceProcess.waitFor();

                    log.info("exited dTraceProcess");
                    
                    if (exit == 0) {
                        String content = getContent(dTraceProcess.getInputStream());
                        log.info(content);
                    } else {
                        log.info("throw exception for dtrace exit=" + exit);
                        // throw exception that leads to latch countdown and hence stop of dtrace
                        throw new IOException("DTrace exit=" + exit);
                    }

                } catch (Exception e) {
                    if (e.getMessage() != null) {
                        log.severe("Error: " + e.getMessage());
                        if (dTraceStartedLatch.getCount() > 0) {
                            dTraceStartedLatch.countDown();
                        }
                    }
                }
            }

        };
         
        task.setPriority(Thread.MAX_PRIORITY);
        task.start();
    }

    protected static String getContent(InputStream input) throws IOException {
        if (input == null) {
            return null;
        }

        byte[] b = new byte[1024];

        int readBytes = 0;

        StringBuilder result = new StringBuilder();

        while ((readBytes = input.read(b)) >= 0) {
            result.append(new String(b, 0, readBytes));
        }
        return result.toString();
    }

    protected void stopDTrace() {
        
        log.info("stopDTrace");

        dTraceProcess.destroy();

        task.interrupt();
        
        dtraceStopped = true;
        
        log.info("DTrace stopped");
    }

    public static void main(String[] args) {
        try {
            Run_dtrace_Test test = new Run_dtrace_Test("Run_dtrace_Test");
            test.test1();
            int z = 1;
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
