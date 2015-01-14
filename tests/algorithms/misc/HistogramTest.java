package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.logging.Logger;

import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class HistogramTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

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

        aa = new float[]{0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.4f, 0.4f, 0.4f, 
            0.4f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};

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

        Histogram.calulateHistogramBinErrors(xHist, yHist, aa, aae, 
            xHistErrorsOutput, yHistErrorsOutput);
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
        float[] aa = new float[]{2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 
            6, 6, 6, 6, 6};
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
        float tmp = Histogram.calculateHistogramWidthYLimitError(xHist, yHist, 
            xErrors, yErrors, yMaxFactor);
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
            
            ArrayPair r = new ArrayPair( Arrays.copyOf(a, nValues), 
                Arrays.copyOf(ae, nValues));
            
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

    public void testReduceHistogramToFirstPeak() throws Exception {
        
        String[] fileNames = new String[] {
            "Aggregation.txt", "Compound.txt", "Pathbased.txt" , "Spiral.txt",
            "D31.txt", "R15.txt" , "Jain.txt", "Flame.txt",
            "a1.txt", "a2.txt", "a3.txt"
        };
        int[] yPeakIdxes = new int[] {
            4, 5, 3, 7, 
            3, 3, 5, 11,
            3, 3, 3
        };
        int[] yPeakMinimaIdxes = new int[] {
            11, 11, 11, 9,
            10, 11, 11, 18,
            10, 10, 19
        };
        
        for (int i = 0; i < fileNames.length; i++) {
            
            String fileName = fileNames[i];
            
            HistogramHolder hist = getHistogram(fileName);
            
            int yPeakIdx = Histogram.findFirstPeakIndex(hist);
            
            int yPeakMinimaIdx = Histogram.findFirstMinimaFollowingPeak(hist, 
                yPeakIdx);
            
            assertTrue(yPeakIdx == yPeakIdxes[i]);
            //System.out.println(yPeakIdx + ", " + yPeakMinimaIdx + ":" + yPeakMinimaIdxes[i]);
            assertTrue(yPeakMinimaIdx == yPeakMinimaIdxes[i]);  
        }
    }

    public void test0() throws Exception {
        Histogram h = new Histogram();
    }
    
    protected HistogramHolder getHistogram(String fileName) {
        
        HistogramHolder hist = new HistogramHolder();
        
        if (fileName.equals("Aggregation.txt")) {
            
            hist.setXHist(new float[]{0.38890636f, 0.44124335f, 0.49358034f, 
                0.5459173f, 0.59825426f, 0.65059125f, 0.70292825f, 0.75526524f, 
                0.8076022f, 0.8599392f, 0.91227615f, 0.9646132f, 1.0169501f, 
                1.0692872f, 1.1216241f, 1.1739612f, 1.2262981f, 1.2786351f, 
                1.3309721f, 1.3833091f, 1.435646f, 1.4879831f, 1.54032f, 
                1.5926571f, 1.644994f});
            
            hist.setYHistFloat(new float[]{77f, 149f, 205f, 204f, 212f, 195f, 
                195f, 181f, 172f, 174f, 138f, 130f, 134f, 128f, 147f, 116f, 
                159f, 160f, 141f, 134f, 171f, 155f, 140f, 157f, 185f}
            );
            
            hist.setXErrors(new float[]{0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f, 0.026168495f, 
                0.026168495f, 0.026168495f, 0.026168495f}
            );
            
            hist.setYErrors(new float[]{0.009232345f, 0.0074687637f, 
                0.007094023f, 0.007819966f, 0.008376563f, 0.009491833f, 
                0.010212443f, 0.011338231f, 0.012392217f, 0.013110398f, 
                0.015541896f, 0.01688796f, 0.017545085f, 0.018921891f, 
                0.018456193f, 0.021673597f, 0.019400077f, 0.02021657f, 
                0.022346234f, 0.023710761f, 0.021814043f, 0.023816979f, 
                0.025850967f, 0.025245331f, 0.024015011f}
            );
            
        } else if (fileName.equals("Compound.txt")) {
            
            hist.setXHist(new float[]{0.42835927f, 0.59034604f, 0.75233275f, 
                0.91431946f, 1.0763062f, 1.2382929f, 1.4002798f, 1.5622666f, 
                1.7242532f, 1.88624f, 2.0482266f, 2.2102134f, 2.3722003f, 
                2.5341868f, 2.6961737f, 2.8581603f, 3.020147f, 3.182134f, 
                3.3441205f, 3.5061073f, 3.668094f, 3.8300807f, 3.9920676f, 
                4.154054f, 4.316041f}
            );
            
            hist.setYHistFloat(new float[]{86.f, 176.f, 253.f, 265.f, 255.f, 
                289.f, 250.f, 207.f, 166.f, 162.f, 142.f, 120.f, 176.f, 172.f, 
                115.f, 97.f, 48.f, 21.f, 6.f, 1.f, 3.f, 0.f, 0.f, 0.f, 0.f}
            );
            
            hist.setXErrors(new float[]{0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f, 0.080993384f, 
                0.080993384f, 0.080993384f, 0.080993384f}
            );
            
            hist.setYErrors(new float[]{0.011044767f, 0.010449851f, 0.0109357f, 
                0.012686985f, 0.014824005f, 0.015973134f, 0.018977497f, 
                0.022860596f, 0.027863028f, 0.03064927f, 0.0352371f, 
                0.041238185f, 0.036702387f, 0.039611485f, 0.050845973f, 
                0.058226686f, 0.08654159f, 0.13695203f, 0.26727128f, 
                0.69430524f, 0.41582954f, 0.0016532135f, 0.0016532135f, 
                0.0016532135f, 0.0016532135f}
            );
            
        } else if (fileName.equals("Pathbased.txt")) {
            
            hist.setXHist(new float[]{0.4257853f, 0.49367052f, 0.5615558f, 
                0.6294411f, 0.6973263f, 0.7652115f, 0.8330968f, 0.9009821f, 
                0.9688673f, 1.0367526f, 1.1046377f, 1.1725229f, 1.2404083f, 
                1.3082935f, 1.3761787f, 1.444064f, 1.5119492f, 1.5798346f, 
                1.6477197f, 1.7156049f, 1.7834903f, 1.8513755f, 1.9192606f, 
                1.987146f, 2.0550313f}
            );
            
            hist.setYHistFloat(new float[]{58.f, 106.f, 104.f, 114.f, 98.f, 
                86.f, 94.f, 86.f, 84.f, 83.f, 89.f, 82.f, 72.f, 67.f, 60.f, 
                64.f, 74.f, 76.f, 74.f, 58.f, 49.f, 50.f, 40.f, 26.f, 30.f}
            );
            
            hist.setXErrors(new float[]{0.03394261f, 0.03394261f, 0.03394261f, 
                0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 
                0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 
                0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 
                0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 0.03394261f, 
                0.03394261f, 0.03394261f}
            );
            
            hist.setYErrors(new float[]{0.011730701f, 0.010066399f, 
                0.01142342f, 0.012213565f, 0.014438428f, 0.01679357f, 
                0.01751225f, 0.019671658f, 0.021318225f, 0.022934532f, 
                0.023668472f, 0.026051551f, 0.029224776f, 0.031819463f, 
                0.035415683f, 0.036052197f, 0.03501306f, 0.0361575f, 
                0.038158353f, 0.044788245f, 0.050578f, 0.051721167f, 
                0.060159065f, 0.07719813f, 0.074333526f}
            );
            
        } else if (fileName.equals("Spiral.txt")) {
            
            hist.setXHist(new float[]{0.4140958f, 0.44548082f, 0.47686586f, 
                0.5082509f, 0.53963596f, 0.57102096f, 0.60240597f, 0.633791f, 
                0.66517603f, 0.69656104f, 0.7279461f, 0.7593311f, 0.7907161f, 
                0.8221012f, 0.8534862f, 0.8848712f, 0.91625625f, 0.94764125f, 
                0.9790263f, 1.0104113f, 1.0417963f, 1.0731814f, 1.1045663f, 
                1.1359514f, 1.1673365f}
            );
            
            hist.setYHistFloat(new float[]{12.f, 26.f, 31.f, 53.f, 37.f, 61.f, 
                92.f, 166.f, 134.f, 96.f, 110.f, 121.f, 134.f, 153.f, 196.f, 
                230.f, 320.f, 471.f, 713.f, 241.f, 27.f, 1.f, 3.f, 5.f, 1.f}
            );
            
            hist.setXErrors(new float[]{0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f, 0.015692517f, 
                0.015692517f, 0.015692517f, 0.015692517f}
            );
            
            hist.setYErrors(new float[]{0.024243813f, 0.017987486f, 
                0.017602937f, 0.014496157f, 0.01819318f, 0.015133856f, 
                0.013169752f, 0.010439589f, 0.012025907f, 0.014726643f, 
                0.014365254f, 0.0143139595f, 0.014187602f, 0.013855188f, 
                0.01275252f, 0.012232615f, 0.0108399335f, 0.00936545f, 
                0.007978725f, 0.013484718f, 0.039806105f, 0.21334535f, 
                0.12535408f, 0.10019355f, 0.22705628f}
            );
            
        } else if (fileName.equals("D31.txt")) {
            
            hist.setXHist(new float[]{0.4555705f, 0.563395f, 0.67121947f, 
                0.77904403f, 0.8868685f, 0.99469304f, 1.1025175f, 1.210342f, 
                1.3181666f, 1.425991f, 1.5338155f, 1.6416401f, 1.7494646f, 
                1.8572892f, 1.9651135f, 2.072938f, 2.1807625f, 2.288587f, 
                2.3964114f, 2.504236f, 2.6120605f, 2.719885f, 2.8277094f, 
                2.935534f, 3.0433586f, 3.1511831f, 3.2590077f, 3.366832f, 
                3.4746566f, 3.5824811f, 3.6903055f, 3.79813f, 3.9059546f, 
                4.013779f, 4.1216035f, 4.229428f, 4.3372526f, 4.445077f, 
                4.5529013f, 4.660726f}
            );
            
            hist.setYHistFloat( new float[]{972.f, 1907.f, 1934.f, 2349.f, 
                1934.f, 1580.f, 1330.f, 1068.f, 936.f, 857.f, 825.f, 882.f, 
                925.f, 929.f, 1008.f, 1016.f, 995.f, 949.f, 912.f, 908.f, 815.f, 
                807.f, 750.f, 683.f, 744.f, 643.f, 615.f, 572.f, 489.f, 408.f, 
                387.f, 358.f, 355.f, 296.f, 254.f, 280.f, 239.f, 199.f, 195.f, 
                178.f}
            );
            
            hist.setXErrors( new float[]{0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f, 0.053912267f, 0.053912267f, 
                0.053912267f, 0.053912267f}
            );
            
            hist.setYErrors(new float[]{0.003280208f, 0.0028552834f, 
                0.0033263417f, 0.0034723307f, 0.004279499f, 0.0052337823f, 
                0.0062627206f, 0.007604179f, 0.008795256f, 0.0099332305f, 
                0.010843717f, 0.011234753f, 0.011665408f, 0.012330402f, 
                0.0125227235f, 0.013146416f, 0.013953508f, 0.014978957f, 
                0.015964845f, 0.016715772f, 0.018358396f, 0.019203408f, 
                0.020682639f, 0.02247743f, 0.022326805f, 0.024848318f, 
                0.026240822f, 0.028098948f, 0.03132962f, 0.035332084f, 
                0.03732318f, 0.039911482f, 0.04121992f, 0.04639081f, 
                0.051403314f, 0.050190173f, 0.05567299f, 0.0625567f, 
                0.06467474f, 0.069301546f}
            );
            
        } else if (fileName.equals("R15.txt")) {
            
            hist.setXHist(new float[]{0.62060326f, 0.73796844f, 0.8553337f, 
                0.97269887f, 1.090064f, 1.2074293f, 1.3247944f, 1.4421597f, 
                1.5595249f, 1.67689f, 1.7942553f, 1.9116205f, 2.0289857f, 
                2.146351f, 2.2637162f, 2.3810813f, 2.4984467f, 2.6158118f, 
                2.733177f, 2.8505423f, 2.9679074f, 3.0852726f, 3.202638f, 
                3.320003f, 3.4373682f}
            );
            
            hist.setYHistFloat(new float[]{159.f, 95.f, 249.f, 397.f, 178.f, 
                96.f, 101.f, 134.f, 117.f, 102.f, 105.f, 62.f, 77.f, 68.f, 
                71.f, 70.f, 89f, 101f, 93f, 92f, 98f, 116f, 131f, 121f, 135f}
            );
            
            hist.setXErrors(new float[]{0.05868259f, 0.05868259f, 0.05868259f, 
                0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 
                0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 
                0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 
                0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 0.05868259f, 
                0.05868259f, 0.05868259f}
            );
            
            hist.setYErrors(new float[]{0.0107018985f, 0.015640663f, 
                0.011573987f, 0.010270698f, 0.016768077f, 0.024951924f, 
                0.02662559f, 0.025208334f, 0.029132104f, 0.033262078f, 
                0.035021987f, 0.04855885f, 0.04599282f, 0.051813576f, 
                0.053645186f, 0.056534085f, 0.052718095f, 0.051854666f, 
                0.05629769f, 0.059160825f, 0.059598666f, 0.056964062f, 
                0.055779833f, 0.05999243f, 0.058679536f}
            );
            
        } else if (fileName.equals("Jain.txt")) {
            
            hist.setXHist(new float[]{0.45988864f, 0.6571544f, 0.8544201f, 
                1.0516859f, 1.2489518f, 1.4462174f, 1.6434833f, 1.840749f, 
                2.0380147f, 2.2352805f, 2.4325461f, 2.6298118f, 2.8270779f, 
                3.0243435f, 3.221609f, 3.4188747f, 3.6161408f, 3.8134065f, 
                4.010672f, 4.207938f, 4.405204f, 4.6024694f, 4.799735f, 
                4.9970007f, 5.194267f}
            );
            
            hist.setYHistFloat(new float[]{228f, 221f, 190f, 217f, 224f, 253f, 
                214f, 155f, 151f, 173f, 90f, 75f, 57f, 40f, 28f, 30f, 25f, 5f, 
                13f, 2f, 7f, 0f, 4f, 1f, 0f}
            );
            
            hist.setXErrors(new float[]{0.09863287f, 0.09863287f, 0.09863287f, 
                0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 
                0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 
                0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 
                0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 0.09863287f, 
                0.09863287f, 0.09863287f}
            );
            
            hist.setYErrors(new float[]{0.009425132f, 0.011626688f, 0.015341333f, 
                0.017103666f, 0.019535251f, 0.020987337f, 0.024980692f, 
                0.031921536f, 0.035407007f, 0.036211606f, 0.053336043f, 
                0.062634915f, 0.07647853f, 0.09658514f, 0.121567294f, 
                0.1241457f, 0.14483534f, 0.34087843f, 0.21985154f, 0.58915895f, 
                0.3301306f, 0.0031047345f, 0.46462873f, 0.9887259f, 
                0.0031047345f}
            );
            
        } else if (fileName.equals("Flame.txt")) {
            
            hist.setXHist(new float[]{0.62069994f, 0.7673387f, 0.91397744f, 
                1.0606163f, 1.2072549f, 1.3538938f, 1.5005324f, 1.6471713f, 
                1.7938099f, 1.9404488f, 2.0870876f, 2.2337263f, 2.380365f, 
                2.5270038f, 2.6736426f, 2.820281f, 2.96692f, 3.1135588f, 
                3.2601976f, 3.406836f, 3.553475f, 3.7001138f, 3.8467526f, 
                3.993391f, 4.14003f}
            );
            
            hist.setYHistFloat(new float[]{23f, 88f, 93f, 115f, 132f, 109f, 
                117f, 115f, 128f, 132f, 123f, 138f, 122f, 102f, 59f, 32f, 14f, 
                8f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
            );
            
            hist.setXErrors(new float[]{0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f, 0.073319376f, 
                0.073319376f, 0.073319376f, 0.073319376f}
            );
            
            hist.setYErrors(new float[]{0.026985954f, 0.016974032f, 
                0.019325659f, 0.020193089f, 0.021366263f, 0.0261085f, 
                0.027905164f, 0.030773634f, 0.03189279f, 0.033850543f, 
                0.03759002f, 0.038075577f, 0.04302981f, 0.04993721f, 
                0.06874155f, 0.09792941f, 0.15489633f, 0.21386705f, 
                7.798004E-4f, 7.798004E-4f, 7.798004E-4f, 7.798004E-4f, 
                7.798004E-4f, 7.798004E-4f, 7.798004E-4f}
            );
            
        } else if (fileName.equals("a1.txt")) {
            
            hist.setXHist(new float[]{0.00935694f, 0.011871975f, 0.014387009f, 
                0.016902043f, 0.019417077f, 0.021932112f, 0.024447147f, 
                0.02696218f, 0.029477216f, 0.03199225f, 0.034507286f, 
                0.03702232f, 0.03953735f, 0.04205239f, 0.04456742f, 
                0.047082458f, 0.04959749f, 0.052112523f, 0.05462756f, 
                0.057142593f, 0.05965763f, 0.062172662f, 0.0646877f, 
                0.06720273f, 0.069717765f, 0.0722328f, 0.07474784f, 0.07726287f, 
                0.079777904f, 0.08229294f, 0.08480798f, 0.08732301f, 0.08983804f, 
                0.092353076f, 0.09486811f, 0.09738315f, 0.09989818f, 
                0.102413215f, 0.10492825f, 0.10744329f}
            );
            
            hist.setYHistFloat(new float[]{569f, 1060f, 1769f, 1846f, 1537f, 
                1595f, 1311f, 1052f, 896f, 814f, 786f, 826f, 852f, 879f, 962f, 
                950f, 960f, 871f, 964f, 897f, 939f, 866f, 871f, 748f, 751f, 
                700f, 649f, 653f, 615f, 552f, 488f, 461f, 418f, 391f, 354f, 
                329f, 339f, 261f, 253f, 198f}
            );
            
            hist.setXErrors(new float[]{0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f, 0.0012575174f, 0.0012575174f, 
                0.0012575174f, 0.0012575174f}
            );
            
            hist.setYErrors(new float[]{8.577941E-5f, 8.077057E-5f, 
                7.34335E-5f, 8.410244E-5f, 1.0420043E-4f, 1.1481581E-4f, 
                1.3960237E-4f, 1.702588E-4f, 2.005605E-4f, 2.2762375E-4f, 
                2.4968886E-4f, 2.6075472E-4f, 2.7399117E-4f, 2.8670713E-4f, 
                2.9037052E-4f, 3.0807292E-4f, 3.2241817E-4f, 3.5512922E-4f, 
                3.5383727E-4f, 3.831793E-4f, 3.9085434E-4f, 4.2366184E-4f, 
                4.3909895E-4f, 4.9153727E-4f, 5.089278E-4f, 5.4542144E-4f, 
                5.857012E-4f, 6.0360454E-4f, 6.4116344E-4f, 6.9816754E-4f, 
                7.647089E-4f, 8.095671E-4f, 8.737772E-4f, 9.286401E-4f, 
                0.0010016178f, 0.0010661648f, 0.0010778557f, 0.0012571953f, 
                0.0013089585f, 0.0015139922f}
            );
            
        } else if (fileName.equals("a2.txt")) {
            
            hist.setXHist(new float[]{0.009381481f, 0.011838005f, 0.014294529f, 
                0.016751053f, 0.019207576f, 0.0216641f, 0.024120625f, 
                0.026577149f, 0.029033672f, 0.031490196f, 0.03394672f, 
                0.036403242f, 0.038859766f, 0.04131629f, 0.043772813f, 
                0.046229336f, 0.04868586f, 0.051142383f, 0.053598907f, 
                0.05605543f, 0.058511954f, 0.060968477f, 0.063425004f, 
                0.06588153f, 0.06833806f, 0.07079458f, 0.073251106f, 
                0.07570763f, 0.07816415f, 0.08062068f, 0.0830772f, 0.08553372f, 
                0.08799025f, 0.09044677f, 0.09290329f, 0.09535982f, 0.09781634f, 
                0.100272864f, 0.10272939f, 0.10518591f}
            );
            
            hist.setYHistFloat(new float[]{1380f, 3243f, 3964f, 4475f, 3680f, 
                3111f, 2552f, 1904f, 1599f, 1372f, 1282f, 1372f, 1333f, 1475f, 
                1615f, 1602f, 1601f, 1534f, 1576f, 1583f, 1542f, 1540f, 1451f, 
                1413f, 1254f, 1253f, 1162f, 1125f, 1089f, 986f, 909f, 832f, 
                830f, 727f, 645f, 644f, 559f, 541f, 488f, 432f}
            );
            
            hist.setXErrors(new float[]{0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f, 0.0012282617f, 0.0012282617f, 
                0.0012282617f, 0.0012282617f}
            );
            
            hist.setYErrors(new float[]{5.6087247E-5f, 4.616087E-5f, 
                4.960274E-5f, 5.4017593E-5f, 6.713565E-5f, 8.123371E-5f, 
                9.867547E-5f, 1.2457524E-4f, 1.4775095E-4f, 1.7246323E-4f, 
                1.9198225E-4f, 1.9864693E-4f, 2.1487378E-4f, 2.1703863E-4f, 
                2.1978703E-4f, 2.3260189E-4f, 2.4515323E-4f, 2.6248163E-4f, 
                2.712859E-4f, 2.8283207E-4f, 2.9894907E-4f, 3.1133762E-4f, 
                3.3349742E-4f, 3.5061763E-4f, 3.8584252E-4f, 3.9951445E-4f, 
                4.2890297E-4f, 4.5039074E-4f, 4.7264723E-4f, 5.112962E-4f, 
                5.489404E-4f, 5.90328E-4f, 6.0785393E-4f, 6.6699093E-4f, 
                7.270518E-4f, 7.4637827E-4f, 8.2143315E-4f, 8.5589435E-4f, 
                9.224865E-4f, 0.0010042344f}
            );
            
        } else if (fileName.equals("a3.txt")) {
            
            hist.setXHist(new float[]{0.008774442f, 0.010016886f, 0.0112593295f, 
                0.0125017725f, 0.013744216f, 0.01498666f, 0.016229104f, 
                0.017471548f, 0.018713992f, 0.019956436f, 0.02119888f, 
                0.022441324f, 0.023683768f, 0.024926212f, 0.026168656f, 
                0.0274111f, 0.028653543f, 0.029895987f, 0.031138431f, 
                0.032380875f, 0.03362332f, 0.034865763f, 0.036108207f, 
                0.03735065f, 0.03859309f, 0.039835535f, 0.04107798f, 
                0.042320423f, 0.043562867f, 0.04480531f, 0.046047755f, 
                0.0472902f, 0.048532642f, 0.049775086f, 0.05101753f, 
                0.052259974f, 0.053502418f, 0.054744862f, 0.055987306f, 
                0.05722975f}
            );
            
            hist.setYHistFloat(new float[]{1315f, 2521f, 3032f, 3623f, 3402f, 
                3135f, 3103f, 3154f, 2961f, 2878f, 2512f, 2297f, 1967f, 1587f, 
                1431f, 1184f, 1171f, 1005f, 1002f, 906f, 912f, 950f, 990f, 976f, 
                923f, 1026f, 1013f, 1078f, 1155f, 1167f, 1167f, 1106f, 1116f, 
                1154f, 1070f, 1155f, 1098f, 1156f, 1139f, 1140f}
            );
            
            hist.setXErrors(new float[]{6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 6.2122196E-4f, 
                6.2122196E-4f, 6.2122196E-4f}
            );
            
            hist.setYErrors(new float[]{5.1306324E-5f, 4.2015632E-5f, 
                4.299305E-5f, 4.3372645E-5f, 4.8949325E-5f, 5.521418E-5f, 
                5.988406E-5f, 6.384715E-5f, 7.030728E-5f, 7.580886E-5f, 
                8.5873966E-5f, 9.475232E-5f, 1.07729364E-4f, 1.258811E-4f, 
                1.3884724E-4f, 1.5965596E-4f, 1.6760733E-4f, 1.8853376E-4f, 
                1.9659598E-4f, 2.1477323E-4f, 2.2228323E-4f, 2.2581349E-4f, 
                2.2895084E-4f, 2.3846551E-4f, 2.5329072E-4f, 2.4811437E-4f, 
                2.5725443E-4f, 2.569715E-4f, 2.555256E-4f, 2.6144116E-4f, 
                2.6852306E-4f, 2.8322372E-4f, 2.8930558E-4f, 2.9176503E-4f, 
                3.1023935E-4f, 3.0606487E-4f, 3.2122963E-4f, 3.2034193E-4f, 
                3.30011E-4f, 3.3713735E-4f}
            );
           
        }
        
        return hist;
    }
    
    public void testCreateBinaryHistogram() throws Exception {
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        img.setValue(3, 1);
        img.setValue(6, 1);
        img.setValue(9, 1);
        img.setValue(12, 1);
        
        int[] ns = Histogram.createBinaryHistogram(img);
        
        assertNotNull(ns);
        assertTrue(ns[0] == (img.getNPixels() - 4));
        assertTrue(ns[1] == 4);
    }

}
