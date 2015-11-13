package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntensityFeaturesTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void test2() throws Exception {

        // 2
        int cellDim = IntensityFeatures.getDefaultCellDimForExtract();
        
        // 6
        int nCellsAcross = IntensityFeatures.getDefaultNCellsAcrossForExtract();
        int len = IntensityFeatures.getDefaultLengthForCellExtractOffsets();
    
        float[] xTrq0 = new float[len];
        float[] yTrq0 = new float[xTrq0.length];
        
        int rotationInDegrees = 45;
        
        IntensityFeatures.populateRotationOffsetsQ0(cellDim, nCellsAcross,
            rotationInDegrees, xTrq0, yTrq0);

        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        int n2 = cellDim * cellDim;
        int idx = 0;
        int oIdx = 0;
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < n2; ++i) {
                    float xOff = xTrq0[idx];
                    float yOff = yTrq0[idx];
                    sb.append(String.format("(%.1f,%.1f) ", xOff, yOff));
                    if (xOff==0 && yOff==0) {
                        log.info("center at output index=" + oIdx);
                    }
                    idx++;
                }
                oIdx++;
                log.info(String.format("%d   (%d, %d): %s", idx, dx, dy, sb.toString()));
            }
        }
    }
 
    public void test0() throws Exception {
        
        int cellDim = 2;
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        float[] output = new float[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];
        
        int count = 0;
        
        int xCenter = 10;
        int yCenter = 10;
        
        int w = 20;
        int h = 20;
        
        int rotation = 45;
        
        int tc = 0;
                
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                
                // --- calculate values for the cell ---
                boolean withinBounds = IntensityFeatures.transformCellCoordinates(
                    rotation, xCenter, yCenter, dx, dy, cellDim, w, h, xT, yT);
                
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < xT.length; ++i) {
                    sb.append(String.format("(%.1f,%.1f) ", xT[i] - xCenter, 
                        yT[i] - yCenter));
                }
                                
                log.info(String.format("prev %d %d (%d, %d): %s", count, tc, 
                    dx, dy, sb.toString()));
                
                tc += xT.length;
                
                count++;
            }
        }
    }
}
