package algorithms.compGeometry;

import algorithms.util.PairFloatArray;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class EllipseHelperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EllipseHelperTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testFitEllipseToPoints() throws Exception {
                
        PairFloatArray xyPoints = getEllipse();
        
        EllipseHelper curveHelper = new EllipseHelper();
        double[] params = curveHelper.fitEllipseToPointsWithAlgLSQ(xyPoints);
        
        double a = params[2];
        double b = params[3];
        
        assertTrue(Math.abs(a - 2) < 0.01);
        assertTrue(Math.abs(b - 1) < 0.01);
    }
    
    public PairFloatArray getEllipse0() {
        /*
        from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, & Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        
        points on ellipse
            a/b=2
        */
        PairFloatArray xy = new PairFloatArray();
        xy.add(3.8760f,  9.8104f);
        xy.add(19.1326f,  2.9130f);
        xy.add(-8.6445f, -9.0177f);
        xy.add(-8.5955f, -9.0294f);
        xy.add(18.3397f, -3.9892f);
        xy.add(-14.8771f,  6.6834f);
        xy.add(-19.8514f, -1.2169f);
        xy.add(9.7412f, -8.7337f);
        
        return xy;
    }
    
    public PairFloatArray getEllipse() {
        // generated from 0 to 360 degrees w/ 10 intervals 
        // and a = 2, b = 1
        PairFloatArray xy = new PairFloatArray();
        xy.add(2.0f, 0.0f);
        xy.add(1.96961550602f, 0.173648177667f);
        xy.add(1.87938524157f, 0.342020143326f);
        xy.add(1.73205080757f, 0.5f);
        xy.add(1.53208888624f, 0.642787609687f);
        xy.add(1.28557521937f, 0.766044443119f);
        xy.add(1.0f, 0.866025403784f);
        xy.add(0.684040286651f, 0.939692620786f);
        xy.add(0.347296355334f, 0.984807753012f);
        xy.add(1.22464679915e-16f, 1.0f);
        xy.add(-0.347296355334f, 0.984807753012f);
        xy.add(-0.684040286651f, 0.939692620786f);
        xy.add(-1.0f, 0.866025403784f);
        xy.add(-1.28557521937f, 0.766044443119f);
        xy.add(-1.53208888624f, 0.642787609687f);
        xy.add(-1.73205080757f, 0.5f);
        xy.add(-1.87938524157f, 0.342020143326f);
        xy.add(-1.96961550602f, 0.173648177667f);
        xy.add(-2.0f, 1.22464679915e-16f);
        xy.add(-1.96961550602f, -0.173648177667f);
        xy.add(-1.87938524157f, -0.342020143326f);
        xy.add(-1.73205080757f, -0.5f);
        xy.add(-1.53208888624f, -0.642787609687f);
        xy.add(-1.28557521937f, -0.766044443119f);
        xy.add(-1.0f, -0.866025403784f);
        xy.add(-0.684040286651f, -0.939692620786f);
        xy.add(-0.347296355334f, -0.984807753012f);
        xy.add(-3.67394039744e-16f, -1.0f);
        xy.add(0.347296355334f, -0.984807753012f);
        xy.add(0.684040286651f, -0.939692620786f);
        xy.add(1.0f, -0.866025403784f);
        xy.add(1.28557521937f, -0.766044443119f);
        xy.add(1.53208888624f, -0.642787609687f);
        xy.add(1.73205080757f, -0.5f);
        xy.add(1.87938524157f, -0.342020143326f);
        xy.add(1.96961550602f, -0.173648177667f);
        return xy;
    }    
   
}
