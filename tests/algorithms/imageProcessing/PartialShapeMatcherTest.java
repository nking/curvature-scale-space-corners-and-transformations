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
 
    public void estDescriptors() throws Exception {
        
        PairIntArray p = getWineGlassShape();
        
        plot(p);
        
        PairIntArray q = p.copy();
        q.rotateLeft(q.getN() - 3);
        
        //TODO: add occlusion to q
        
        
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        /*
        float[][] a = shapeMatcher.createDescriptorMatrix(p);    
        for (int i = 0; i < a.length; ++i) {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("row=%3d: ", i));
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(
                String.format("[%2d] %.4f,", j, a[i][j]));
            }
            System.out.println(sb.toString());
        }*/
        
        shapeMatcher.match(p, q);
        
    }
    
    public void testDescriptors2() throws Exception {
        
        PairIntArray p = getScissors1();
        //plot(p);
        
        PairIntArray q = getScissors2();
        //plot(q);
        
        System.out.println("p.n=" + p.getN() 
            + " q.n=" + q.getN());
        
        q.rotateLeft(q.getN() - 3);
        
        //TODO: add occlusion to q
        
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        shapeMatcher.match(q, p);        
    }
    
    protected PairIntArray getScissors1() {
        
        PairIntArray p = new PairIntArray();
        p.add(95,58); p.add(102,52); p.add(108,42);
        p.add(110,35); p.add(118,28); p.add(128,25);
        p.add(135,25); p.add(142,27); p.add(152,32);
        p.add(150,40); p.add(142,48); p.add(135,52);
        p.add(128,58); p.add(119,59); p.add(110,60);
        p.add(102,62); p.add(95,65); p.add(100,72);
        p.add(108,80); p.add(115,87); p.add(122,95);
        p.add(128,104); p.add(130,111); p.add(130,120);
        p.add(125,120); p.add(115,122); p.add(108,122);
        p.add(100,112); p.add(95,105); p.add(97,97);
        p.add(97,89); p.add(95,80); p.add(89,72);
        p.add(80,70); p.add(71,70); p.add(65,70);
        p.add(57,70); p.add(49,70); p.add(40,70);
        p.add(32,70); p.add(23,66); p.add(20,65);
        p.add(29,63); p.add(37,63); p.add(45,63);
        p.add(52,60); p.add(60,60); p.add(70,60);
        p.add(61,55); p.add(52,51); p.add(45,45);
        p.add(38,40); p.add(28,35); p.add(38,32);
        p.add(46,35); p.add(55,40); p.add(62,42);
        p.add(70,48); p.add(78,52); p.add(86,60);
        
        // accidently entered y for x, so
        // reverse them here
        for (int i = 0; i < p.getN(); ++i) {
            int x = p.getX(i);
            int y = p.getY(i);
            p.set(i, y, x);
        }
        
        return p;
    }
    
    protected PairIntArray getScissors2() {
        
        PairIntArray p = new PairIntArray();
        p.add(80,105); p.add(75,112);
        p.add(65,117); p.add(58,121);
        p.add(52,130); p.add(51,138);
        p.add(50,147); p.add(55,155);
        p.add(60,160); p.add(66,156);
        p.add(72,150); p.add(78,141);
        p.add(82,133); p.add(85,125);
        p.add(85,116); p.add(90,110);
        p.add(92,100); p.add(100,95);
        p.add(108,91); p.add(116,88);
        p.add(125,85); p.add(132,82);
        p.add(141,78); p.add(149,70);
        p.add(155,62); p.add(150,60);
        p.add(142,52); p.add(135,52);
        p.add(125,53); p.add(118,58);
        p.add(112,65); p.add(110,75);
        p.add(105,82); p.add(97,87);
        p.add(97,78); p.add(97,70);
        p.add(97,60); p.add(97,52);
        p.add(97,45); p.add(95,36);
        p.add(90,30); p.add(90,40);
        p.add(90,48); p.add(90,55);
        p.add(88,62); p.add(87,70);
        p.add(70,93); p.add(62,92);
        p.add(52,95); p.add(45,95);
        p.add(37,97); p.add(28,98);
        p.add(20,100); p.add(27,103);
        p.add(37,107); p.add(47,107);
        p.add(52,107); p.add(62,105);
        p.add(68,104);
       
        return p;
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
