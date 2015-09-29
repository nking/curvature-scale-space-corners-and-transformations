package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.ArrayPair;

import java.util.logging.Logger;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 *
 * @author nichole
 */
public class TwoPointVoidStats2Test extends BaseTwoPointTest {

    boolean debug = true;

    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    /**
     *
     * @throws Exception
     */
    public void testTopCentroid() throws Exception {

        /* (0)                     | (1)                     | (2)                   | (3)
         *8          peak          |           peak          |         peak          |            peak
         *7
         *6                        |     pt_i+1    pt_j-1    |                       |
         *5    pt_i        pt_j    |   pt_i           pt_j   |                       |
         *4   ------------------   |  ---------------------  |  -------------------  |  ----------------------
         *3                        |                         |                       |
         *2  pt_0             pt_k | pt_0              pt_k  |   pt_0      pt_k      |  <no point>      <no point>
         *                         |                         |                       |
         *   2  3  4  5  6  7  8     2  3  4  5  6  7  8  9      2  3  4  5  6  7       0  2  3  4  5  6  8
         *
         *                                                                           |8(4)     peak
         *                                                                           |7
         *                                                                           |6     pt_i       pt_j
         *                                                                           |5
         *                                                                           |4 ----------------------
         *                                                                           | <no point>      <no point>
         */

        // case 0
        float[] x = new float[]{2.5f, 3.0f, 5.0f, 7.0f, 8.0f};
        float[] y = new float[]{2.0f, 5.0f, 8.0f, 5.1f, 2.0f};
        // assert results internal to calculateCentroidOfTop first
        ArrayPair resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, 
            null, null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 6);
        assertTrue(resultXY.getY().length == 6);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[4] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 5.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 5.1f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[0] - 2.83f) < 0.1f);
        assertTrue(Math.abs(resultXY.getX()[1] - 3.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 5.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 7.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[4] - 7.35f) < 0.1f);
        assertTrue(Math.abs(resultXY.getX()[5] - resultXY.getX()[0]) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[5] - resultXY.getY()[0]) < 0.01f);

        //area, xc, yc
        float[] a = TwoPointVoidStats.calculateCentroidOfTop(x, y, 0.5f);
        assertNotNull(a);
        assertTrue(Math.abs(a[1] - 5.0f) < 0.1f);
        assertTrue(Math.abs(a[2] - 6.0f) < 1.0f);


        // case 1
        x = new float[]{2.2f,   3.0f, 4.0f,   5.0f,  6.8f, 8.0f,   8.5f};
        y = new float[]{2.0f,   5.0f, 6.0f,   8.0f,  6.0f, 5.1f,   2.0f};
        resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, 
            null, null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 8);
        assertTrue(resultXY.getY().length == 8);
        assertTrue(Math.abs(resultXY.getX()[0] - 2.73f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[1] - 3.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 5.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[4] - 6.8f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[5] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[6] - 8.18f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[7] - 2.73f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 5.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 6.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[4] - 6.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[5] - 5.1f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[6] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[7] - 4.0f) < 0.01f);

        // case 2
        x = new float[]{2.5f, 4.5f, 6.0f};
        y = new float[]{2.0f, 8.0f, 2.0f};
        resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, null,
            null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 4);
        assertTrue(resultXY.getY().length == 4);
        assertTrue(Math.abs(resultXY.getX()[0] - 3.17f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[1] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 5.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 3.17f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 4.0f) < 0.01f);

        // case 4
        x = new float[]{4.0f, 5.5f, 7.5f};
        y = new float[]{6.0f, 8.0f, 6.0f};
        resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, null,
            null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 6);
        assertTrue(resultXY.getY().length == 6);
        assertTrue(Math.abs(resultXY.getX()[0] - 2.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[1] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 5.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 7.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[4] - 9.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[5] - 2.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 6.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 6.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[4] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[5] - 4.0f) < 0.01f);

        // case 3 left
        x = new float[]{4.5f, 5.5f, 7.8f};
        y = new float[]{8.0f, 6.0f, 2.0f};
        resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, null,
            null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 5);
        assertTrue(resultXY.getY().length == 5);
        assertTrue(Math.abs(resultXY.getX()[0] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[1] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 5.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 6.65f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[4] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 6.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[4] - 4.0f) < 0.01f);

        // case 3 right
        x = new float[]{3.5f, 4.5f};
        y = new float[]{6.5f, 8.0f};
        resultXY = LinesAndAngles.createPolygonOfTopFWFractionMax(x, y, null,
            null, 0.5f);
        assertNotNull(resultXY);
        assertTrue(resultXY.getX().length == 5);
        assertTrue(resultXY.getY().length == 5);
        assertTrue(Math.abs(resultXY.getX()[0] - 1.83f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[1] - 3.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[2] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[3] - 4.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getX()[4] - 1.83f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[0] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[1] - 6.5f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[2] - 8.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[3] - 4.0f) < 0.01f);
        assertTrue(Math.abs(resultXY.getY()[4] - 4.0f) < 0.01f);

    }
}
