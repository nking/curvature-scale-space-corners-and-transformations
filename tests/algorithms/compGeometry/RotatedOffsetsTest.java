package algorithms.compGeometry;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RotatedOffsetsTest extends TestCase {
    
    public RotatedOffsetsTest() {
    }
    
    public void test0() throws Exception {
        
        RotatedOffsets ro = RotatedOffsets.getInstance();
        
        for (int r = 0; r < 360; ++r) {
        
            int[] xExpected = getExpectedXOffsets(r);
            int[] yExpected = getExpectedYOffsets(r);
            
            int[] xOffsets = ro.getXOffsets(r);
            int[] yOffsets = ro.getYOffsets(r);
        
            assertTrue(xOffsets.length == xExpected.length);
            assertTrue(yOffsets.length == yExpected.length);
            
            
            for (int i = 0; i < xExpected.length; ++i) {
                int ae = xExpected[i];
                int a = xOffsets[i];
                if (ae != a) {
                    int z = 1;
                }
                int be = yExpected[i];
                int b = yOffsets[i];
                if (be != b) {
                    int z = 1;
                }
            }
            
            assertTrue(Arrays.equals(xExpected, xOffsets));
            assertTrue(Arrays.equals(yExpected, yOffsets));
        }
    }
    
    public int[] getExpectedXOffsets(int rotationInDegrees) {
        
        int cellDimension = 2;
        int range0 = 6;
    
        double rotationInRadians = rotationInDegrees * Math.PI/180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);
            
        int[] outputX = new int[144];
        
        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {
                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {
                        double xt = (((dx + dxc) * mc) + ((dy + dyc) * ms));
                        outputX[count] = (int)Math.round(xt);
                        count++;
                    }
                }
            }
        }
        
        return outputX;
    }
    
    public int[] getExpectedYOffsets(int rotationInDegrees) {
        
        int cellDimension = 2;
        int range0 = 6;
    
        double rotationInRadians = rotationInDegrees * Math.PI/180.;
        double mc = Math.cos(rotationInRadians);
        double ms = Math.sin(rotationInRadians);
            
        int[] outputY = new int[144];
        
        int count = 0;
        for (int dx = -range0; dx < range0; dx += cellDimension) {
            for (int dy = -range0; dy < range0; dy += cellDimension) {
                for (int dxc = 0; dxc < cellDimension; ++dxc) {
                    for (int dyc = 0; dyc < cellDimension; ++dyc) {
                        double yt = (-((dx + dxc) * ms)) + ((dy + dyc) * mc);
                        outputY[count] = (int)Math.round(yt);
                        count++;
                    }
                }
            }
        }
        
        return outputY;
    }
}
