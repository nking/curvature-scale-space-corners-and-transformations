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
    
    public void testCalculateEllipseResiduals() throws Exception {
        
        EllipseHelper ellipseHelper = new EllipseHelper();
       
        PairFloatArray points = getEllipse();
        float[] xP = new float[points.getN()];
        float[] yP = new float[points.getN()];
        for (int i = 0; i < xP.length; i++) {
            xP[i] = points.getX(i);
            yP[i] = points.getY(i);
        }
        
        double[] params = ellipseHelper.fitEllipseToPointsWithAlgLSQ(points);
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        double[] stats = ellipseHelper.calculateEllipseResidualStats(xP, yP,
            xc, yc, a, b, alpha);
        
        assertTrue(Math.abs(stats[0] - 0) < 1e-6);
        
        assertTrue(Math.abs(stats[1] - 0) < 1e-6);        
    }
    
    public void testCalculateEllipseResiduals1() throws Exception {
        
        EllipseHelper ellipseHelper = new EllipseHelper();
       
        PairFloatArray points = getEllipse0();
        float[] xP = new float[points.getN()];
        float[] yP = new float[points.getN()];
        for (int i = 0; i < xP.length; i++) {
            xP[i] = points.getX(i);
            yP[i] = points.getY(i);
        }
        
        double[] params = ellipseHelper.fitEllipseToPoints(points);
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        // a=20.0 ,b=10.0 , xc=0 , yc=0 , alpha=0
        double[] stats = ellipseHelper.calculateEllipseResidualStats(xP, yP,
            xc, yc, a, b, alpha);
        
        assertTrue(Math.abs(stats[0] - 0) < 1e-5);
        
        assertTrue(Math.abs(stats[1] - 0) < 1e-5);        
    }
    
    public void testCalculateEllipseResiduals4() throws Exception {
        
        EllipseHelper ellipseHelper = new EllipseHelper();
       
        PairFloatArray points = getEllipse4();
        float[] xP = new float[points.getN()];
        float[] yP = new float[points.getN()];
        for (int i = 0; i < xP.length; i++) {
            xP[i] = points.getX(i);
            yP[i] = points.getY(i);
        }
        
        double[] params = ellipseHelper.fitEllipseToPoints(points);
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        double[] stats = ellipseHelper.calculateEllipseResidualStats(xP, yP,
            xc, yc, a, b, alpha);
        
        assertTrue(Math.abs(stats[0] - 0) < 1e-4);
        assertTrue(Math.abs(stats[1] - 0) < 1e-3);
    }
    
    public void testCalculateEllipseResiduals2() throws Exception {
        
        EllipseHelper ellipseHelper = new EllipseHelper();
       
        PairFloatArray points = getEllipse2();
        float[] xP = new float[points.getN()];
        float[] yP = new float[points.getN()];
        for (int i = 0; i < xP.length; i++) {
            xP[i] = points.getX(i);
            yP[i] = points.getY(i);
        }
        
        double[] params = ellipseHelper.fitEllipseToPoints(points);
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        // a=25.208, b=8.4729, (xc, yc)=(-4.468, 2.7686), alpha=-0.27766
        double[] stats = ellipseHelper.calculateEllipseResidualStats(xP, yP,
            xc, yc, a, b, alpha);
        
        assertTrue(Math.abs(stats[0] - 0) < 0.1);
        
        /*
        to check this visually:
        
        from math import *
        from pylab import *

        def ellipse():
            a = 25.208
            b = 8.4729
            xc = -4.468
            yc = 2.7686
            alpha = -0.27766
            ca = math.cos(alpha)
            sa = math.sin(alpha)
            x = []
            y = []
            for angle in range(0,360, 10):
                t = angle*math.pi/180.
                ct = math.cos(t)
                st = math.sin(t)
                g = xc + (a * ca * ct) - (b * sa * st)
                x.append(g)
                h = yc + (a * sa * ct) + (b * ca * st)
                y.append(h)

            xx = [-0.7783, 14.6672, -8.3475, -6.8840, 13.4167, -16.0429, -24.1829, 8.9161]
            yy = [11.6781, 3.8028, -4.7133, -5.5677, -3.7199, 2.6031, 0.3222, -9.5737]
            plot(x,y)
            scatter(xx,yy, color = 'red')
            show()

        ellipse()
        */  
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
    
    public PairFloatArray getEllipse2() {
        /*
        from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, & Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        
        perturbed ellipse points
        */
        PairFloatArray xy = new PairFloatArray();
        xy.add(-0.7783f, 11.6781f);
        xy.add(14.6672f, 3.8028f);
        xy.add(-8.3475f, -4.7133f);
        xy.add(-6.8840f, -5.5677f);
        xy.add(13.4167f, -3.7199f);
        xy.add(-16.0429f, 2.6031f);
        xy.add(-24.1829f, 0.3222f);
        xy.add(8.9161f, -9.5737f);
        
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
   
    public PairFloatArray getEllipse4() {
        // generated from 0 to 360 degrees w/ 40 intervals 
        // and a = 2, b = 1, angle=30 degrees
        
        float[] x = new float[]{401.7320508075689f, 401.0054340914946f, 
            399.80836358985476f, 398.70096189432337f, 398.20139456563845f, 
            398.54341470896406f, 399.5669872981078f, 400.79317134286697f, 
            401.64822170118117f};
        float[] y= new float[]{126.0f, 126.32271484234539f, 126.02651670961937f, 
            125.25f, 124.35650551194011f, 123.76410924648808f, 123.75f, 
            124.3207796457145f, 125.20937404389255f};
        
        PairFloatArray xy = new PairFloatArray();
        
        for (int i = 0; i < x.length; i++) {
            xy.add(x[i], y[i]);
        }
        
        return xy;
    }    
}
