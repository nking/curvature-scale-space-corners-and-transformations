package thirdparty.edu.princeton.cs.algs4;

import junit.framework.TestCase;


/**
 * these tests are from the main method of StdRandom.java.
 * Please see the license information there.
 * 
 * @author nichole
 */
public class StdRandomTest extends TestCase {
    
    public StdRandomTest() {
    }

    public void test0() {
        
        /*
        %  java StdRandom 5
 *  seed = 1316600602069
 *  59 16.81826  true 8.83954  0 
 *  32 91.32098  true 9.11026  0 
 *  35 10.11874  true 8.95396  3 
 *  92 32.88401  true 8.87089  0 
 *  72 92.55791  true 9.46241  0 
 *
 *  % java StdRandom 5
 *  seed = 1316600616575
 *  96 60.17070  true 8.72821  0 
 *  79 32.01607  true 8.58159  0 
 *  81 59.49065  true 9.10423  1 
 *  96 51.65818  true 9.02102  0 
 *  99 17.55771  true 8.99762  0 
 *
 *  % java StdRandom 5 1316600616575
 *  seed = 1316600616575
 *  96 60.17070  true 8.72821  0 
 *  79 32.01607  true 8.58159  0 
 *  81 59.49065  true 9.10423  1 
 *  96 51.65818  true 9.02102  0 
 *  99 17.55771  true 8.99762  0 
 *
 *
        */
        
        for (int t = 0; t < 3; ++t) {

            String[] args = null;
            if (t == 0) {
                args = new String[]{"5"};
            } else if (t == 1) {
                args = new String[]{"5"};
            } else if (t == 2) {
                args = new String[]{"5", "1316600616575"};
            }
            
            int n = Integer.parseInt(args[0]);
            if (args.length == 2) StdRandom.setSeed(Long.parseLong(args[1]));
            double[] probabilities = { 0.5, 0.3, 0.1, 0.1 };
            int[] frequencies = { 5, 3, 1, 1 };
            String[] a = "A B C D E F G".split(" ");

            System.out.println("seed = " + StdRandom.getSeed());
            for (int i = 0; i < n; i++) {
                System.out.printf("%2d ",   StdRandom.uniform(100));
                System.out.printf("%8.5f ", StdRandom.uniform(10.0, 99.0));
                System.out.printf("%5b ",   StdRandom.bernoulli(0.5));
                System.out.printf("%7.5f ", StdRandom.gaussian(9.0, 0.2));
                System.out.printf("%1d ",   StdRandom.discrete(probabilities));
                System.out.printf("%1d ",   StdRandom.discrete(frequencies));
                StdRandom.shuffle(a);
                for (String s : a)
                    System.out.print(s);
                System.out.println();
            }

        }
        
    }
    
}
