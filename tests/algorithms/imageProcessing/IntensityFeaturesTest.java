package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntensityFeaturesTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
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
    
    public void testCalculate45DegreeOrientation() throws Exception {
        
        /*
               |3|3|2|1|1|
               |3|3|2|1|1|
               |4|4| |0|0|
               |5|5|6|7|7|
               |5|5|6|7|7|
        */
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        int blockHalfWidths = 5;
        
        boolean useNormalizedIntensities = true;
        
        GreyscaleImage img = new GreyscaleImage(25, 25, 
            GreyscaleImage.Type.Bits32FullRangeInt);
                
        boolean[] sign = new boolean[]{true, false};
        /*
               |3|3|2|1|1|
               |3|3|2|1|1|
               |4|4| |0|0|
               |5|5|6|7|7|
               |5|5|6|7|7|
        */
        int xc = 2;
        int yc = 2;
        for (boolean pos : sign) {
            
            for (int i = 1; i < 2; ++i) {
            
                // values are cached so create a new instance for each "img"
                IntensityFeatures features = new IntensityFeatures(blockHalfWidths,
                    useNormalizedIntensities, rOffsets);
        
                int v = i;
                if (!pos) {
                    v *= -1;
                }
                
                img.setValue(xc + 1, yc, v);
                img.setValue(xc + 2, yc, v);
                
                v = i + 1;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc + 1, yc + 1, v);
                img.setValue(xc + 2, yc + 1, v);
                img.setValue(xc + 1, yc + 2, v);
                img.setValue(xc + 2, yc + 2, v);
                
                v = i + 2;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc, yc + 1, v);
                img.setValue(xc, yc + 2, v);
                /*
                   |3|3|2|1|1|
                   |3|3|2|1|1|
                   |4|4| |0|0|
                   |5|5|6|7|7|
                   |5|5|6|7|7|
                */
                
                v = i + 3;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc + 1, v);
                img.setValue(xc - 2, yc + 1, v);
                img.setValue(xc - 1, yc + 2, v);
                img.setValue(xc - 2, yc + 2, v);
                
                v = i + 4;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc, v);
                img.setValue(xc - 2, yc, v);
                
                v = i + 5;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc - 1, yc - 1, v);
                img.setValue(xc - 2, yc - 1, v);
                img.setValue(xc - 1, yc - 2, v);
                img.setValue(xc - 2, yc - 2, v);
                /*
                   |3|3|2|1|1|
                   |3|3|2|1|1|
                   |4|4| |0|0|
                   |5|5|6|7|7|
                   |5|5|6|7|7|
                */
                
                v = i + 6;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc, yc - 1, v);
                img.setValue(xc, yc - 2, v);
                
                v = i + 7;
                if (v > 7) {
                    v = v - 8;
                }
                if (!pos) {
                    v *= -1;
                }
                img.setValue(xc + 1, yc - 1, v);
                img.setValue(xc + 2, yc - 1, v);
                img.setValue(xc + 1, yc - 2, v);
                img.setValue(xc + 2, yc - 2, v);
                
                features.calculateGradientWithGreyscale(img);
                
                int orientation = features.calculate45DegreeOrientation(2, 2);
                
                switch(i) {
                    case 0: {
                        if (pos) {
                            assertEquals(135, orientation);
                        } else {
                            assertEquals(315, orientation);
                        }
                        break;
                    }
                    case 1: {
                        if (pos) {
                            assertEquals(90, orientation);
                        } else {
                            assertEquals(270, orientation);
                        }
                        break;
                    }
                    case 2: {
                        if (pos) {
                            assertEquals(45, orientation);
                        } else {
                            assertEquals(225, orientation);
                        }
                        break;
                    }
                    case 3: {
                        if (pos) {
                            assertEquals(0, orientation);
                        } else {
                            assertEquals(180, orientation);
                        }
                        break;
                    }
                    case 4: {
                        if (pos) {
                            assertEquals(315, orientation);
                        } else {
                            assertEquals(135, orientation);
                        }
                        break;
                    }
                    case 5: {
                        if (pos) {
                            assertEquals(270, orientation);
                        } else {
                            assertEquals(90, orientation);
                        }
                        break;
                    }
                    case 6: {
                        if (pos) {
                            assertEquals(225, orientation);
                        } else {
                            assertEquals(45, orientation);
                        }
                        break;
                    } 
                    default: {
                        if (pos) {
                            assertEquals(180, orientation);
                        } else {
                            assertEquals(0, orientation);
                        }
                        break;
                    }
                }
            }
        }
        
    }
}
