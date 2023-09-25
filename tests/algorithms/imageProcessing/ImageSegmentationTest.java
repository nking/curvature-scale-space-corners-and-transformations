package algorithms.imageProcessing;

import algorithms.misc.Misc0;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TIntList;
import java.io.IOException;
import java.util.ArrayList;
import algorithms.imageProcessing.segmentation.*;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageSegmentationTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ImageSegmentationTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
    }
  
    private int count(List<Set<PairInt>> setList) {
        
        int c = 0;
        for (Set<PairInt> set : setList) {
            c += set.size();
        }
        
        return c;
    }

    private void plotFFT(GreyscaleImage crImg, String fileNameRoot) throws IOException {

        int bn = 1;//8
        float[] xPoints = new float[crImg.getWidth() / bn];
        float[] yPoints = new float[crImg.getWidth() / bn];
        float xmn = 0;
        float xmx = crImg.getWidth() / bn;
        float ymn = Float.MAX_VALUE;
        float ymx = Float.NEGATIVE_INFINITY;
        int row = 50;
        for (int i = 0; i < (crImg.getWidth() / bn) - 1; ++i) {
            int ii = bn * i;
            xPoints[i] = i;
            for (int k = 0; k < bn; ++k) {
                yPoints[i] += crImg.getValue(ii + k, row);
            }
            yPoints[i] /= bn;
            if (yPoints[i] < ymn) {
                ymn = yPoints[i];
            }
            if (yPoints[i] > ymx) {
                ymx = yPoints[i];
            }
        }

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(xmn, xmx, ymn, ymx,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
            "X fft_" + fileNameRoot);

        xPoints = new float[crImg.getHeight() / bn];
        yPoints = new float[crImg.getHeight() / bn];
        xmn = 0;
        xmx = crImg.getHeight() / bn;
        ymn = Float.MAX_VALUE;
        ymx = Float.NEGATIVE_INFINITY;
        int col = 50;
        for (int j = 0; j < (crImg.getHeight() / bn) - 1; ++j) {
            int jj = bn * j;
            xPoints[j] = j;
            for (int k = 0; k < bn; ++k) {
                yPoints[j] += crImg.getValue(jj + k, row);
            }
            yPoints[j] /= bn;
            if (yPoints[j] < ymn) {
                ymn = yPoints[j];
            }
            if (yPoints[j] > ymx) {
                ymx = yPoints[j];
            }
        }

        plotter.addPlot(xmn, xmx, ymn, ymx,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
            "y fft_" + fileNameRoot);

        plotter.writeFile(fileNameRoot + "_blob_cr_fft");

    }
   
}
