package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PartialShapeMatcherTest extends TestCase {
    
    public PartialShapeMatcherTest() {
    }
 
    public void testDescriptots() throws Exception {
        
        PairIntArray p = getWineGlassShape();
        
        plot(p);
        
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        float[][] a = shapeMatcher.createDescriptorMatrix(p);
    
        /*
        for (int i = 0; i < a.length; ++i) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("row=%3d: ", i));
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(
                String.format("[%2d] %.4f,", j, a[i][j]));
            }
            System.out.println(sb.toString());
        }*/
    }
    
    protected PairIntArray getWineGlassShape() {
        
        PairIntArray p = new PairIntArray();
       p.add(140, 500 - 450);
       p.add(220, 500 - 410);
       p.add(225, 500 - 390);
       p.add(225, 500 - 325);
       p.add(190, 500 - 280);
       p.add(130, 500 - 240);
       p.add(130, 500 - 180);
       p.add(200, 500 - 180);
       p.add(250, 500 - 180);
       p.add(310, 500 - 180);
       p.add(380, 500 - 200);
       p.add(330, 500 - 260);
       p.add(280, 500 - 310);
       p.add(275, 500 - 380);
       p.add(310, 500 - 440);
       p.add(375, 500 - 460);
        p.add(320, 500 - 475);
        p.add(275, 500 - 475);
        p.add(200, 500 - 475);
        return p;
    }

    private void plot(PairIntArray p) throws Exception {

        float[] x = new float[p.getN()];
        float[] y = new float[p.getN()];
        
        for (int i = 0; i < x.length; ++i) {
            x[i] = p.getX(i);
            y[i] = p.getY(i);
        }
        
        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;
        
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        
        plot.addPlot(0, xMax, 0, yMax, 
            x, y, x, y, "");
        
        plot.writeFile(200);
    }
}
