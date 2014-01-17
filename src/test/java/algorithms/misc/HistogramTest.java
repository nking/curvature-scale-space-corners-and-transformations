package algorithms.misc;

import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class HistogramTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testPersistedHistograms() throws Exception {
        
        log.info("testPersistedHistograms");
        
        // tests to set up for regression tests for a range of datasets
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        File[] data = ResourceFinder.findFilesInTestResources("void_densities_");
        
        if (data != null) {
            
            for (File file : data) {
                String label = "";
                String fn = file.getName();
                Pattern p = Pattern.compile("void\\_densities\\_(.*?)\\.txt");
                Matcher m = p.matcher(fn);
                if (m.matches()) {
                    label = m.group(1);
                }
                
                ArrayPair pair = readTestFile(file);
                
                float[] a = pair.getX();
                float[] ae = pair.getY();

                int n = a.length;

                HistogramHolder hist = Histogram.defaultHistogramCreator(a, ae);

                plotter.addPlot(hist.getXHist(), hist.getYHistFloat(), new float[n], new float[n], label);
                plotter.writeFile2();

                boolean threwException = false;
                try {
                    Histogram.defaultHistogramCreator(null, ae);
                } catch (IllegalArgumentException e) {
                    threwException = true;
                }
                assertTrue(threwException);
                
                
                threwException = false;
                try {
                    Histogram.defaultHistogramCreator(a, null);
                } catch (IllegalArgumentException e) {
                    threwException = true;
                }
                assertTrue(threwException);
            }
        }
    }
    
    /**
     * Test of createHistogram method, of class Histogram.
     */
    public void testCreateHistogram_4args() {

        log.info("testCreateHistogram_4args");

        float[] aa = new float[]{1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5};
        int nBins = 5;

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        Histogram.createHistogram(aa, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 1));
            assertTrue(xHist[i] == (i + 1 + 0.5));
        }

        aa = new float[]{0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4};
        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(aa, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 1));
            assertTrue(xHist[i] == (i + 0.5));
        }

        aa = new float[]{0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.4f, 0.4f, 0.4f, 0.4f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};

        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(aa, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            float expected = (i + 1);
            float found = yHist[i];
            assertTrue(expected == found);

            expected = (float) ((i + 1 + 0.5)*0.1f);
            found = xHist[i];
            assertTrue( Math.abs(expected - found) < 0.01);
        }

        aa = new float[]{
            100, 100, 100, 100, 100,
            200, 200, 200, 200, 200, 200,
            300, 300, 300, 300, 300, 300, 300,
            400, 400, 400, 400, 400, 400, 400, 400,
            500, 500, 500, 500, 500, 500, 500, 500, 500,
        };

        float[] aae = new float[aa.length];
        for (int i = 0; i < aa.length; i++) {
            aae[i] = 0.1f*(float)Math.sqrt(aa[i]);
        }

        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(aa, nBins, xHist, yHist);
        for (int i = 0; i < yHist.length; i++) {
            float expected = (i + 5);
            float found = yHist[i];
            assertTrue(expected == found);

            expected = (float) ((i + 1 + 0.5)*100.f);
            found = xHist[i];
            assertTrue( Math.abs(expected - found) < 0.01);
        }

        float[] xHistErrorsOutput = new float[xHist.length];
        float[] yHistErrorsOutput = new float[xHist.length];

        Histogram.calulateHistogramBinErrors(xHist, yHist, aa, aae, xHistErrorsOutput, yHistErrorsOutput);
        for (int i = 0; i < yHist.length; i++) {
            float yh = yHist[i];
            float xh = xHist[i];
            float xhe = xHistErrorsOutput[i];
            float yhe = yHistErrorsOutput[i];
            //TODO:  add assert of rough expected values
            assertTrue(xhe < xh);
            assertTrue(yhe < yh);
        }

        
        boolean threwException = false;
        try {
            Histogram.createHistogram(null, nBins, xHistErrorsOutput, yHist);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(aa, nBins, null, yHist);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(null, nBins, xHistErrorsOutput, null);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
    }

    /**
     * Test of createHistogram method, of class Histogram.
     */
    public void testCreateHistogram_6args() {

        log.info("testCreateHistogram_6args");

        //  2   3   4   5   6
        //    0   1   2   3   4
        float[] aa = new float[]{2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6};
        int nBins = 5;

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        Histogram.createHistogram(aa, nBins, 2, 6, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 2));
            assertTrue(xHist[i] == (i + 2 + 0.5));
        }

        nBins = 4;
        xHist = new float[nBins];
        yHist = new int[nBins];
        Histogram.createHistogram(aa, nBins, 2, 5, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 2));
            assertTrue(xHist[i] == (i + 2 + 0.5));
        }
        
        boolean threwException = false;
        try {
            Histogram.createHistogram(null, nBins, 2, 5, xHist, yHist);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(aa, nBins, 2, 5, null, yHist);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(aa, nBins, 2, 5, xHist, null);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        

        threwException = false;
        try {
            Histogram.createHistogram(null, nBins, 2, 5, xHist, yHist, 2);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(aa, nBins, 2, 5, null, yHist, 2);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            Histogram.createHistogram(aa, nBins, 2, 5, xHist, null, 2);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
    }
    
    public void testCalculateHistogramWidthYLimitError() throws Exception {
        System.out.println("testCalculateHistogramWidthYLimitError()");
        
        /* 1.0|     *
         *    |        *
         * 0.5|  *     
         *    |           *
         *   0|_______________
         *    0  1  2  3  4
         */
        float[] xHist = new float[]{1.0f, 2.0f, 3.0f,  4.0f};
        float[] yHist = new float[]{0.5f, 1.0f, 0.75f, 0.25f};
        float[] xErrors = new float[]{0.1f, 0.1f, 0.1f, 0.1f};
        float[] yErrors = new float[]{0.1f, 0.1f, 0.1f, 0.1f};
        float yMaxFactor = 1;
        float tmp = Histogram.calculateHistogramWidthYLimitError(xHist, yHist, xErrors, yErrors, yMaxFactor);
        assertTrue(Math.abs(tmp - 0.141f) < 0.01f);
    }

    protected ArrayPair readTestFile(String fileName) throws Exception {
       
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        File file = new File(filePath);
        
        return readTestFile(file);
    }
    
    protected ArrayPair readTestFile(File file) throws Exception {
        
        float[] a = null;
        float[] ae = null;
        int nValues = 0;
        
        FileReader reader = null;
        BufferedReader in = null;

        try {
            int count = 0;
            reader = new FileReader(file);
            in = new BufferedReader(reader);

            String line = in.readLine();
            while (line != null) {
                if (count == 0) {
                    nValues = Integer.valueOf(line);
                    a = new float[nValues];
                    ae = new float[nValues];
                } else {
                    String[] values = line.split("\t");
                    a[count-1] = Float.valueOf(values[0]);
                    ae[count - 1] = Float.valueOf(values[1]);
                }
                line = in.readLine();
                count++;
            }
            
            ArrayPair r = new ArrayPair( Arrays.copyOf(a, nValues), Arrays.copyOf(ae, nValues));
            
            return r;
            
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (in != null) {
                in.close();
            }
        }
    }

    public void testReadWriteExternal() throws Exception {

        float[] aa = new float[]{
            2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11,
            12, 13, 14, 15
        };
        float[] aae = Errors.populateYErrorsBySqrt(aa);
        
        boolean threwException = false;
        try {
            HistogramHolder histogram = Histogram.createSimpleHistogram(null, aae);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        threwException = false;
        try {
            HistogramHolder histogram = Histogram.createSimpleHistogram(aa, null);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        threwException = false;
        try {
            HistogramHolder histogram = Histogram.createSimpleHistogram(10, null, aae);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        threwException = false;
        try {
            HistogramHolder histogram = Histogram.createSimpleHistogram(10, aa, null);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);

        
        HistogramHolder histogram = Histogram.createSimpleHistogram(aa, aae);


        PipedOutputStream pipedOut = null;
        PipedInputStream pipedIn = null;

        ObjectOutputStream oos = null;

        try {

            pipedIn = new PipedInputStream(1024);
            pipedOut = new PipedOutputStream(pipedIn);


            final CountDownLatch writeLatch = new CountDownLatch(1);
            final CountDownLatch doneLatch = new CountDownLatch(1);

            Reader reader = new Reader(writeLatch, doneLatch, pipedIn);
            reader.start();

            oos = new ObjectOutputStream(pipedOut);

            histogram.writeExternal(oos);

            pipedOut.close();
            pipedOut = null;
            writeLatch.countDown();

            doneLatch.await();

            HistogramHolder rHistogram = reader.getHistogramHolder();

            assertTrue( Arrays.equals(histogram.xHist, rHistogram.xHist));

            assertTrue( Arrays.equals(histogram.yHist, rHistogram.yHist));

            assertTrue( Arrays.equals(histogram.yHistFloat, rHistogram.yHistFloat));

            assertTrue( Arrays.equals(histogram.xErrors, rHistogram.xErrors));

            assertTrue( Arrays.equals(histogram.yErrors, rHistogram.yErrors));

        } finally {

            if (oos != null) {
                oos.close();
            }
            if (pipedOut != null) {
                pipedOut.close();
            }
            if (pipedIn != null) {
                pipedIn.close();
            }
        }

    }

    private class Reader extends Thread {

        private final CountDownLatch writeLatch;
        private final CountDownLatch doneLatch;

        protected HistogramHolder histogram = null;

        protected final PipedInputStream pipedIn;

        Reader(CountDownLatch writeLatch, CountDownLatch doneLatch, PipedInputStream pipedInputStream) {
            this.writeLatch = writeLatch;
            this.doneLatch = doneLatch;
            this.pipedIn = pipedInputStream;
        }
        @Override
        public void run() {

            ObjectInputStream in = null;

            try {
                writeLatch.await();

                in = new ObjectInputStream(pipedIn);

                histogram = new HistogramHolder();

                histogram.readExternal(in);

            } catch (Exception e) {
                fail(e.getMessage());
            } finally {
                if (in != null) {
                    try {
                        in.close();
                    } catch (IOException e) {}
                }
                doneLatch.countDown();
            }
        }

        public HistogramHolder getHistogramHolder() {
            return histogram;
        }
    }

    public void test0() throws Exception {
        Histogram h = new Histogram();
    }
}
