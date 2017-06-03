package algorithms.random;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class QuasiMonteCarloTest extends TestCase {

    public QuasiMonteCarloTest() {
    }

    public void testMC_0() {

        /*
        from chap 1 of 
         "Monte Carlo theory, methods and examples"
            by Owen, 2013,
        http://statweb.stanford.edu/~owen/mc/
        
        function is for the average distance between 
            randomly chosen points in a region.
         */
        // rectangle is 4 X 5
        int a = 4;
        int b = 5;
        int nPoints = 100000;
        Dist distFunc = new Dist(a, b);

        MonteCarlo mc = new MonteCarlo();
        double avg = mc.f(distFunc, nPoints);

        // expected:
        double g = expectedAvg(a, b);

        System.out.println("mc avg=" + avg + " expected=" + g);

        assertTrue(Math.abs(g - avg) < 0.05 * g);
    }

    public void testQMC_0() {

        /*
        from chap 1 of 
         "Monte Carlo theory, methods and examples"
            by Owen, 2013,
        http://statweb.stanford.edu/~owen/mc/
        
        function is for the average distance between 
            randomly chosen points in a region.
         */
        // rectangle is 4 X 5
        int a = 4;
        int b = 5;
        int nPoints = 100;
        Dist distFunc = new Dist(a, b);

        QuasiMonteCarlo qmc = new QuasiMonteCarlo();
        double avg0 = qmc.haltonAdvanced(distFunc, nPoints);
        double avg = qmc.lds(distFunc, nPoints);

        // expected:
        double g = expectedAvg(a, b);

        System.out.println("qmc \navg0=" + avg0 + "\navg=" + avg
            + "\nexpected=" + g);

        assertTrue(Math.abs(g - avg0) < 0.05 * g);
    }

    private double expectedAvg(int a, int b) {
        double a2 = a * a;
        double b2 = b * b;
        double a3 = a * a2;
        double b3 = b * b2;
        double g
            = (1. / 15.)
            * ((a3 / b2) + (b3 / a2)
            + (Math.sqrt(a2 + b2) * (3 - (a2 / b2) - (b2 / a2))))
            + (1. / 6.)
            * ((b2 / a) * arccosh(Math.sqrt(a2 + b2) / b)
            + (a2 / b) * arccosh(Math.sqrt(a2 + b2) / a));
        return g;
    }

    private double arccosh(double t) {
        return Math.log(t + Math.sqrt(t * t - 1));
    }

    private class Dist extends AFunction {

        private int a, b;

        public Dist(int rectDim1, int rectDim2) {
            super(4, 0);
            this.a = rectDim1;
            this.b = rectDim2;
        }

        @Override
        public double f(double[] coords) {

            double dx = a * coords[0] - a * coords[2];
            double dy = b * coords[1] - b * coords[3];

            double dist = Math.sqrt(dx * dx + dy * dy);

            return dist;
        }

        @Override
        public double[] der(double[] doubles) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

    }
}
