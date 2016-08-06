package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.util.Arrays;
import java.util.List;
import algorithms.imageProcessing.PartialShapeMatcher.Sequences;
import algorithms.imageProcessing.PartialShapeMatcher.Sequence;
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
 
    public void testDescriptors() throws Exception {
        
        PairIntArray p = getWineGlassShape();
        
        plot(p, 101);
        
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
    
    public void testSummedAreaTables() {
        
        /*
        2 9  2  7     9 11  18   16 22  36
        1 5  1  2     5  6   8    7 11  18
        0 2  3  5     2  5  10    2  5  10
          0  1  2      
        */
        
        float[][] a = new float[3][];
        a[0] = new float[]{2, 3, 5};
        a[1] = new float[]{5, 1, 2};
        a[2] = new float[]{9, 2, 7};
        
        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.applySummedAreaTableConversion(a);
        
        assertTrue(Arrays.equals(new float[]{2, 5, 10}, a[0]));
        assertTrue(Arrays.equals(new float[]{7, 11, 18}, a[1]));
        assertTrue(Arrays.equals(new float[]{16, 22, 36}, a[2]));
        
    }
    
    public void testDescriptors2() throws Exception {
        
        PairIntArray p = getScissors1();
        //plot(p, 200);
        
        PairIntArray q = getScissors2();
        //plot(q, 201);
        
        System.out.println("p.n=" + p.getN() 
            + " q.n=" + q.getN());
        
        //q.rotateLeft(q.getN() - 3);
                
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        // articulated:
        Sequences sequences = shapeMatcher.match(p, q);        
    
        /*
        unless improve the first image blade positions:
        expecting roughly:
         [junit] frac=0.7541, avgDiff=0.5967,  sumDiff=27.4475
            [junit] (0:0 to 14, f=0.2459 d=0.1464)
            [junit] (17:17 to 38, f=0.3607 d=0.6057)
            [junit] (48:51 to 59, f=0.1475 d=1.3251)

       for alt sorting which includes diff angles:
         [junit] frac=0.7049, avgDiff=0.4138,  sumDiff=17.7920
         [junit] (0:0 to 14, f=0.2459 d=0.1464)
         [junit] (17:17 to 38, f=0.3607 d=0.6057)
         [junit] (40:40 to 45, f=0.0984 d=0.3784)
        */

        assertNotNull(sequences);
        assertTrue(sequences.fractionOfWhole > 0.6);

        // assert the correspondence range

    }

    public void testDescriptors3() throws Exception {

        // rotate points p so that start points are 
        // different and assert that wrap around is
        // handled correctly
        
        PairIntArray p = getScissors1();
        p.rotateLeft(16);
        //plot(p, 200);
        
        PairIntArray q = getScissors2();
        //plot(q, 201);
        
        System.out.println("p.n=" + p.getN() 
            + " q.n=" + q.getN());
        
        //q.rotateLeft(q.getN() - 3);
                
        PartialShapeMatcher shapeMatcher = new PartialShapeMatcher();
        shapeMatcher.overrideSamplingDistance(1);
        
        // articulated:
        Sequences sequences = shapeMatcher.match(p, q);        
    
        /*
        unless improve the first image blade positions:
        expecting roughly for case without
        rotation:
            frac=0.7541, avgDiff=0.5967,  sumDiff=27.4475
            (0:0 to 14, f=0.2459 d=0.1464)
            (17:17 to 38, f=0.3607 d=0.6057)
            (48:51 to 59, f=0.1475 d=1.3251)
        rotation left by 16 becomes:
            (45:0 to 14, f=0.2459 d=0.1464)
            (1:17 to 38, f=0.3607 d=0.6057)
            (32:51 to 59, f=0.1475 d=1.3251)

        found:
            frac=0.6721, avgDiff=0.9792,  sumDiff=40.1460
            (0:16 to 38,  f=0.3770 d=1.4907)
            (37:57 to 60, f=0.0656 d=0.2804)
            (45:0 to 13,  f=0.2295 d=0.3385)
        */

        assertNotNull(sequences);
        assertTrue(sequences.fractionOfWhole > 0.6);

        // assert the correspondence range
        List<Sequence> list = sequences.sequences;
        //System.out.println("SEQ0=" + list.get(0).toString());
        //System.out.println("p=" + p.toString());
        //System.out.println("q=" + q.toString());

    }
    
    protected PairIntArray getScissors1() {
        
        PairIntArray p = new PairIntArray();
        p.add(95,55); p.add(102,52); p.add(108,42);
        p.add(110,35); p.add(118,28); p.add(128,25);
        p.add(135,25); p.add(142,27); p.add(152,32);
        p.add(150,40); p.add(142,48); p.add(135,52);
        p.add(128,58); p.add(119,59); p.add(110,60);
        p.add(102,62); p.add(95,65); p.add(100,72);
        p.add(108,80); p.add(115,87); p.add(122,95);      // 18,19,20
        p.add(128,104); p.add(130,111); p.add(130,120);   // 21,22, 23
        p.add(125,120); p.add(116,122); p.add(108,122);
        p.add(100,112); p.add(95,105); p.add(94,98);      // 27
        p.add(96,88); p.add(95,81); p.add(89,72);
        p.add(80,70); p.add(72,70); p.add(63,70);    //33, 34, 35
        p.add(57,70); p.add(49,70); p.add(39,70);
        p.add(32,70); p.add(23,67); p.add(20,64);    // 39, 40, 41
        p.add(28,62); p.add(37,63); p.add(45,62);    // 42
        p.add(53,61); p.add(60,60); p.add(70,60);    // 45
        p.add(62,55); p.add(53,51); p.add(45,45);    // 48
        p.add(38,40); p.add(29,37); p.add(30,33);    // 51
        p.add(38,34); p.add(47,36); p.add(54,40);    // 54
        p.add(62,44); p.add(70,49); p.add(78,52);
        p.add(87,58);
        
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
        p.add(50,147); p.add(55,155);  // 6,7
        p.add(60,160); p.add(67,156);
        p.add(74,149); p.add(78,141);
        p.add(82,133); p.add(85,125);
        p.add(85,116); p.add(89,109);
        p.add(92,100); p.add(100,94);   // 16,17
        p.add(109,91); p.add(117,87);
        p.add(125,85); p.add(132,82);   // 20,21
        p.add(141,78); p.add(149,71);
        p.add(155,63); p.add(150,59);   // 24,25
        p.add(142,52); p.add(135,52);
        p.add(127,53); p.add(118,58);
        p.add(112,66); p.add(110,74);   // 30,31
        p.add(105,82); p.add(97,87);
        p.add(97,77); p.add(97,69);
        p.add(97,61); p.add(97,53);
        p.add(97,45); p.add(95,36);    // 38, 39
        p.add(92,28); p.add(90,30);
        p.add(89,38); p.add(89,47);
        p.add(88,55); p.add(87,62);    // 44, 45
        p.add(87,72); p.add(85,79);    // 46, 47
        p.add(84,88); p.add(78,90);        
        p.add(70,92); p.add(62,93);    // 50, 51
        p.add(52,94); p.add(45,96); 
        p.add(38,97); p.add(30,98); 
        p.add(20,100); p.add(29,104); 
        p.add(38,107); p.add(45,107);  // 58,59
        p.add(52,106); p.add(62,105); 
        p.add(79,103); p.add(77,102);
        
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

    private void plot(PairIntArray p, int fn) throws Exception {

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
        
        plot.writeFile(fn);
    }

}
