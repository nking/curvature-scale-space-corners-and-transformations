package thirdparty.libquantum;

import algorithms.misc.MiscMath;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GroverTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public GroverTest() {
    }

    public void testRun() {

        log.info("testRun");

        int number = 3;
        int nbits = 3;

        int nTests = 1;

        for (int i = 0; i < nTests; ++i) {

            Grover grover = new Grover();

            grover.run(number, nbits);
        }
    }

    public void testRun_list() {

        log.info("testRun_list2");

        int[] list;
        int r, number, nbits;
        //list = new int[]{0,1,2,3,4,5,6,7,8};
        //list = new int[]{1,2,3,0,5,4,6,7,8};
        list = new int[]{1, 2, 3, 0, 1, 4, 2, 7, 8};

        nbits = 4;

        Grover grover = new Grover();

        number = 3;

        r = grover.run(number, nbits, list);
        System.out.println("search=" + number + " found=" + r);
        assertEquals(number, r);

        number = 8;
        grover = new Grover();
        r = grover.run(number, nbits, list);
        System.out.println("search=" + number + " found=" + r);
        assertEquals(number, r);

        number = 5;
        grover = new Grover();
        r = grover.run(number, nbits, list);
        System.out.println("search=" + number + " found=" + r);
        assertEquals(-1, r);
    }

    public void testSetBits() {

        System.out.println("testSetBits");

        // 1 1 0 1
        //int setBits = (1 << 0) + (1 << 2) + (1 << 3);
        // 0 1 1 0
        //int setBits = (1 << 1) + (1 << 2);
        // 1 1 1 1
        int setBits = (1 << 0) + (1 << 1) + (1 << 2) + (1 << 3);

        /*
        qubits
        0000
        0001
        0100
        0101
        1000
        1001
        1100
        1101
         */
        System.out.println("setBits=" + setBits);
        TIntSet set = new TIntHashSet();
        int nBits = MiscMath.numberOfBits(setBits);
        for (int i = 0; i < nBits; ++i) {
            if ((setBits & (1 << i)) != 0) {
                set.add(1 << i);
            }
        }
        //int[] list = set.toArray(new int[set.size()]);

        int v = Integer.MAX_VALUE;
        System.out.println("v=" + v);
        v ^= (1 << 32);
        System.out.println("hb=" + v);

        int width = 4;
        for (int i = 0; i < 10; ++i) {
            int number = i;

            System.out.println("searching for " + number);

            Grover grover = new Grover();

            int r = grover.run(number, width, setBits);
            //int r = grover.run(number, width, list);

            System.out.println("number=" + number + " set="
                + (set.contains(number)) + " r=" + r);
            if (set.contains(number)) {
                assertEquals(number, r);
            } else {
                assertFalse(r == number);
            }
        }

    }
}
