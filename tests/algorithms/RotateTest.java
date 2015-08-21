package algorithms;

import java.util.Arrays;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class RotateTest extends TestCase {
    
    public RotateTest() {
    }
   
    public void testRotate0() {
       
        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g', (int)'h'
        };
        
        Rotate r = new Rotate();
        r.rotate(a, 3);
        
        int[] expected = new int[]{
            (int)'d', (int)'e', (int)'f', (int)'g', (int)'h', (int)'a', 
            (int)'b', (int)'c'
        };
        
        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testRotate1() {

        int shift = -3;

        System.out.println("SHIFT=" + shift);

        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g'};

        /*
        efgabcd
        */        
        
        Rotate r = new Rotate();
        r.rotate(a, shift);

        int[] expected = new int[]{
            (int)'e', (int)'f', (int)'g',
            (int)'a', (int)'b', (int)'c', (int)'d'
        };

        assertTrue(Arrays.equals(expected, a));
        
        //========================
        
        shift = -1 * (a.length + 3);

        System.out.println("SHIFT=" + shift);

        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g'};

        /*
        efgabcd
        */        
        
        r.rotate(a, shift);
        
        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testRotate2() {
       
        int[] shifts = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
        
        for (int shift : shifts) {
            
            System.out.println("SHIFT=" + shift);
            
            int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
                (int)'f', (int)'g', (int)'h'
            };
            
            Rotate r = new Rotate();
            r.rotate2(a, shift);

            int[] expected;
            
            switch(shift) {
                case 1:
                    expected = new int[]{
                        (int)'b', (int)'c', (int)'d', (int)'e', 
                        (int)'f', (int)'g', (int)'h', (int)'a'
                    };
                    break;
                case 2:
                    expected = new int[]{
                        (int)'c', (int)'d', (int)'e', (int)'f', (int)'g', 
                        (int)'h', (int)'a', (int)'b'
                    };
                    break;
                case 3:
                    expected = new int[]{
                        (int)'d', (int)'e', (int)'f', (int)'g', (int)'h', (int)'a', 
                        (int)'b', (int)'c'
                    };
                    break;
                case 4:
                    expected = new int[]{
                        (int)'e', (int)'f', (int)'g', (int)'h', (int)'a', 
                        (int)'b', (int)'c', (int)'d'
                    };
                    break;
                case 5:
                    expected = new int[]{
                        (int)'f', (int)'g', (int)'h', (int)'a', (int)'b', 
                        (int)'c', (int)'d', (int)'e'
                    };
                    break;
                case 6:
                    expected = new int[]{
                        (int)'g', (int)'h', (int)'a', (int)'b', (int)'c', 
                        (int)'d', (int)'e', (int)'f'
                    };
                    break;
                case 7:
                    expected = new int[]{
                        (int)'h', (int)'a', (int)'b', (int)'c', (int)'d', 
                        (int)'e', (int)'f', (int)'g', 
                    };
                    break;
               default:
                    expected = new int[]{
                        (int)'a', (int)'b', (int)'c', (int)'d', 
                        (int)'e', (int)'f', (int)'g', (int)'h'
                    };
                    break;
            }

            assertTrue(Arrays.equals(expected, a));
        }
    }
    
    public void testRotate4() {
               
        int shift = 6;
            
        System.out.println("SHIFT=" + shift);

        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        Rotate r = new Rotate();
        r.rotate2(a, shift);

        int[] expected = new int[]{
            (int)'c', (int)'d', (int)'a', (int)'b'
        };

        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testReverse() {
               
        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        Rotate r = new Rotate();
        r.reverse(a, 0, a.length - 1);

        int[] expected = new int[]{
            (int)'d', (int)'c', (int)'b', (int)'a'
        };

        assertTrue(Arrays.equals(expected, a));
        
        //------------
        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        r.reverse(a, 1, 3);

        expected = new int[]{
            (int)'a', (int)'d', (int)'c', (int)'b'
        };

        assertTrue(Arrays.equals(expected, a));
        
        //------------
        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        r.reverse(a, 1, 2);

        expected = new int[]{
            (int)'a', (int)'c', (int)'b', (int)'d'
        };

        assertTrue(Arrays.equals(expected, a));
        
        assertTrue(Arrays.equals(expected, a));
        
        //------------
        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        r.reverse(a, 0, 2);

        expected = new int[]{
            (int)'c', (int)'b', (int)'a', (int)'d'
        };

        assertTrue(Arrays.equals(expected, a));
        
        //------------
        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        r.reverse(a, 0, 1);

        expected = new int[]{
            (int)'b', (int)'a', (int)'c', (int)'d'
        };

        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testRotate5() {
               
        int shift = -5;
            
        System.out.println("SHIFT=" + shift);

        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d'};

        /*
        a b c d   0
        d a b c  -1
        c d a b  -2
        b c d a  -3
        a b c d  -4
        d a b c  -5
        */        
        
        Rotate r = new Rotate();
        r.rotate2(a, shift);

        int[] expected = new int[]{
            (int)'d', (int)'a', (int)'b', (int)'c'
        };

        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testRotate6() {
       
        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g', (int)'h'
        };
        
        Rotate r = new Rotate();
        r.rotate2(a, 3);
        
        int[] expected = new int[]{
            (int)'d', (int)'e', (int)'f', (int)'g', (int)'h', (int)'a', 
            (int)'b', (int)'c'
        };
        
        assertTrue(Arrays.equals(expected, a));
    }
    
    public void testRotate7() {

        int shift = -3;

        System.out.println("SHIFT=" + shift);

        int[] a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g'};

        /*
        efgabcd
        */        
        
        Rotate r = new Rotate();
        r.rotate2(a, shift);

        int[] expected = new int[]{
            (int)'e', (int)'f', (int)'g',
            (int)'a', (int)'b', (int)'c', (int)'d'
        };

        assertTrue(Arrays.equals(expected, a));
        
        //========================
        
        shift = -1 * (a.length + 3);

        System.out.println("SHIFT=" + shift);

        a = new int[]{(int)'a', (int)'b', (int)'c', (int)'d', (int)'e',
            (int)'f', (int)'g'};

        /*
        efgabcd
        */        
        
        r.rotate2(a, shift);
        
        assertTrue(Arrays.equals(expected, a));
    }
}
