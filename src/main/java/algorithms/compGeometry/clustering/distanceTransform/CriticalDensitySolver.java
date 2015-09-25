package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CriticalDensitySolver {
    
    private boolean debug = false;
    
    public CriticalDensitySolver() {
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    float findCriticalDensity(int[][] distTrans, int nPoints, int width, 
        int height) {
        
        int w = distTrans.length;
        int h = distTrans[0].length;
        
        float[] values = new float[w*h];
        int count2 = 0;
        for (int i0 = 0; i0 < w; ++i0) {
            for (int j0 = 0; j0 < h; ++j0) {
                int v = distTrans[i0][j0];
                values[count2] = (float)(1./Math.sqrt(v));
                count2++;
            }
        }
        
        float[] vErrors = Errors.populateYErrorsBySqrt(values);

        List<HistogramHolder> histList = new ArrayList<HistogramHolder>();

        //TODO: revise this:
        float xl = 4.0f;
        int nb = 40;
        
        /*
        steps xl down as needed and changes nb to 20 when small
        */
        
        int hc = 0;
        
        while (true) {
            
            HistogramHolder hist = Histogram.createSimpleHistogram(
                0, xl, nb, values, vErrors);
        
            if (debug) {
                try {
                    hist.plotHistogram("clstr", "_cluster_" + hc);
                } catch (IOException ex) {
                    Logger.getLogger(CriticalDensitySolver.class.getName()).
                        log(Level.SEVERE, null, ex);
                }
            }
            
            histList.add(hist);
            
            // TODO: examine hist and decide whether to revise xl and nb
            if (true) {
                break;
            }
            
            hc++;
        }
        
        //TODO: examine the latest histogram in histList for first peak
        throw new UnsupportedOperationException("not yet implemented");
    }
    
}
