package algorithms.compGeometry;

import algorithms.imageProcessing.Image;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import junit.framework.TestCase;
import org.junit.Test;

import java.util.*;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class PerimeterFinder3Test extends TestCase {

    @Test
    public void test100() {
        int w = 640;
        int h = 279;
        int[] idxs = new int[]{103457,104097,104737,103456,104096,104736,104739,104098,104738,104741,104740,102815,
                103455,104095,104735,100894,101534,102174,102814,103454,104094,104734,};
        Set<Integer> points = Arrays.stream(idxs)
                .boxed() // Convert IntStream to Stream<Integer>
                .collect(Collectors.toSet());
        //(414,157)->(414,163), 415,163)->(421,163)  ->(419,163)  418,162,
        PairIntArray boundary = PerimeterFinder3.mooreTracingWithJacob(points, w, h);
        //System.out.printf("boundary=%s\n", boundary.toString());

        int[][] expected = new int[][]{{414, 157},{414, 158},{414, 159},{414, 160},{414, 161},{414, 162}, {414, 163},
            {415, 163},{416, 163},{417, 163},{418, 163},{419, 163},
            {420, 163},{421, 163},{420, 163},{419, 163},
                {418, 162},{417, 161}, {416, 161},{415, 160},{414, 159},{414, 158}};
        assertEquals(expected.length, boundary.getN());
        for (int i = 0; i < boundary.getN(); ++i) {
            assertEquals(expected[i][0], boundary.getX(i));
            assertEquals(expected[i][1], boundary.getY(i));
        }
    }

    @Test
    public void testGetLabeledPointsMap() {

        int[] labels = new int[]{2,0,1,0};
        Map<Integer, Set<Integer>> labelPointsMap = PerimeterFinder3.getLabeledPointsMap(labels);

        assertEquals(3, labelPointsMap.size());
        assertTrue(labelPointsMap.containsKey(0));
        assertEquals(2, labelPointsMap.get(0).size());
        assertTrue(labelPointsMap.get(0).contains(1));
        assertTrue(labelPointsMap.get(0).contains(3));

        assertTrue(labelPointsMap.containsKey(1));
        assertEquals(1, labelPointsMap.get(1).size());
        assertTrue(labelPointsMap.get(1).contains(2));

        assertTrue(labelPointsMap.containsKey(2));
        assertEquals(1, labelPointsMap.get(2).size());
        assertTrue(labelPointsMap.get(2).contains(0));
    }

    @Test
    public void testGetLabelAdjMap() {
        /*
        0 2 2 3
        0 0 1 3

        using 8-neighbor adjacency
        */
        int w = 4;
        int h = 2;
        int[] labels = new int[]{0, 0, 1, 3, 0, 2, 2, 3};
        Map<Integer, Set<Integer>> adjMap = PerimeterFinder3.getLabelAdjMap(labels, w, h);

        assertEquals(4, adjMap.size());
        assertTrue(adjMap.containsKey(0));
        assertTrue(adjMap.containsKey(1));
        assertTrue(adjMap.containsKey(2));
        assertTrue(adjMap.containsKey(3));

        assertEquals(2, adjMap.get(0).size());
        assertTrue(adjMap.get(0).contains(1));
        assertTrue(adjMap.get(0).contains(2));

        assertEquals(3, adjMap.get(1).size());
        assertTrue(adjMap.get(1).contains(0));
        assertTrue(adjMap.get(1).contains(2));
        assertTrue(adjMap.get(1).contains(3));

        assertEquals(3, adjMap.get(2).size());
        assertTrue(adjMap.get(2).contains(0));
        assertTrue(adjMap.get(2).contains(1));
        assertTrue(adjMap.get(2).contains(3));

        assertEquals(2, adjMap.get(3).size());
        assertTrue(adjMap.get(3).contains(1));
        assertTrue(adjMap.get(3).contains(2));
    }

    @Test
    public void test0() {
        int w = 10;
        int h = 9;
        int n = w*h;
        int[] labels = new int[] {
            1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
            1, 3, 3, 3, 3, 3, 3, 3, 3, 0,
            1, 3, 3, 4, 4, 6, 6, 6, 3, 0,
            1, 3, 3, 5, 5, 6, 7, 6, 3, 0,
            1, 3, 3, 5, 6, 6, 6, 6, 3, 0,
            1, 3, 3, 3, 3, 3, 3, 3, 3, 0,
            1, 2, 2, 2, 2, 2, 2, 2, 2, 0,
            1, 2, 2, 2, 2, 2, 2, 2, 2, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        /*
        0 | 1  1  1  0  0  0  0  0  0  0
        1 | 1  3  3  3  3  3  3  3  3  0
        2 | 1  3  3  4  4  6  6  6  3  0
        3 | 1  3  3  5  5  6  7  6  3  0
        4 | 1  3  3  5  6  6  6  6  3  0
        5 | 1  3  3  3  3  3  3  3  3  0
        6 | 1  2  2  2  2  2  2  2  2  0
        7 | 1  2  2  2  2  2  2  2  2  0
        8 | 0  0  0  0  0  0  0  0  0  0
          | -  -  -  -  -  -  -  -  -  -
            0  1  2  3  4  5  6  7  8  9
        */
        int r, c;
        List<Integer> expected0 = new ArrayList<>();
        r = 8;
        for (c = 0; c < w; ++c) {
            expected0.add(r*w + c);
        }
        c = 9;
        for (r = 0; r <= 7; ++r) {
            expected0.add(r*w + c);
        }
        r = 0;
        for (c = 3; c <= 8; ++c) {
            expected0.add(r*w + c);
        }
        Collections.sort(expected0);

        List<Integer> expected1 = new ArrayList<>();
        c = 0;
        for (r = 0; r <= 7; ++r) {
            expected1.add(r*w + c);
        }
        r = 0;
        for (c = 1; c <= 2; ++c) {
            expected1.add(r*w + c);
        }
        Collections.sort(expected1);

        List<Integer> expected2 = new ArrayList<>();
        for (r = 6; r <= 7; ++r) {
            for (c = 1; c <= 8; ++c) {
                expected2.add(r * w + c);
            }
        }
        Collections.sort(expected2);

        List<Integer> expected3 = new ArrayList<>();
        for (r = 1; r <= 5; ++r) {
            for (c = 1; c <= 8; ++c) {
                expected3.add(r*w + c);
            }
        }
        Collections.sort(expected3);

        List<Integer> expected4 = new ArrayList<>();
        for (r = 2; r <= 2; ++r) {
            for (c = 3; c <= 4; ++c) {
                expected4.add(r*w + c);
            }
        }
        Collections.sort(expected4);

        List<Integer> expected5 = new ArrayList<>();
        for (r = 3; r <= 4; ++r) {
            for (c = 3; c <= 3; ++c) {
                expected5.add(r*w + c);
            }
        }
        expected5.add(3*w + 4);
        Collections.sort(expected5);

        List<Integer> expected6 = new ArrayList<>();
        for (r = 2; r <= 4; ++r) {
            for (c = 5; c <= 7; ++c) {
                expected6.add(r*w + c);
            }
        }
        expected6.add(4*w + 4);
        Collections.sort(expected6);

        List<Integer> expected7 = new ArrayList<>();
        expected7.add(3*w + 6);

        List<List<Integer>>  expectedList = new ArrayList<>();
        expectedList.add(expected0);
        expectedList.add(expected1);
        expectedList.add(expected2);
        expectedList.add(expected3);
        expectedList.add(expected4);
        expectedList.add(expected5);
        expectedList.add(expected6);
        expectedList.add(expected7);

        assertEquals(n, labels.length);

        Map<Integer, Set<Integer>> filledRegionsMap = PerimeterFinder3.regionFill(labels, w, h);
        assertEquals(8, filledRegionsMap.size());

        BiPredicate<List<Integer>, List<Integer>> listsAreEqual = List::equals;
        for (Map.Entry<Integer, Set<Integer>> entry : filledRegionsMap.entrySet()) {
            List<Integer> region = new ArrayList<>(entry.getValue());
            Collections.sort(region);
            int label = entry.getKey();
            List<Integer> expected = expectedList.get(label);
            assertEquals(expected.size(), region.size());
            assertTrue(listsAreEqual.test(expected, region));
        }

        // ======   test moore tracing ========
        /*
        8 | 0  0  0  0  0  0  0  0  0  0
        7 | 1  2  2  2  2  2  2  2  2  0
        6 | 1  2  2  2  2  2  2  2  2  0
        5 | 1  3  3  3  3  3  3  3  3  0
        4 | 1  3  3  5  6  6  6  6  3  0
        3 | 1  3  3  5  5  6  7  6  3  0
        2 | 1  3  3  4  4  6  6  6  3  0
        1 | 1  3  3  3  3  3  3  3  3  0
        0 | 1  1  1  0  0  0  0  0  0  0
          | -  -  -  -  -  -  -  -  -  -
            0  1  2  3  4  5  6  7  8  9
        */
        List<PairIntArray> expectedBounds = new ArrayList<>();
        for (int i = 0; i < 8; ++i) {
            expectedBounds.add(new PairIntArray());
        }
        int i;
        // region0 =====
        i = 0;
        r = 8;
        for (c = 0; c <= 9; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        c = 9;
        for (r = 7; r >= 0; --r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 0;
        for (c = 8; c >= 3; --c) {
            expectedBounds.get(i).add(c, r);
        }
        for (c = 4; c <= 8; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        c = 9;
        for (r = 1; r <= 7; ++r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 8;
        for (c = 8; c >= 1; --c) {
            expectedBounds.get(i).add(c, r);
        }
        assertEquals(44, expectedBounds.get(i).getN());

        // region1 =======
        i = 1;
        c = 0;
        for (r = 0; r <= 7; ++r) {
            expectedBounds.get(i).add(c, r);
        }
        for (r = 6; r >= 1; --r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 0;
        for (c = 1; c <= 2; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        expectedBounds.get(i).add(1, 0);

        // region2 =======
        i = 2;
        c = 1;
        for (r = 6; r <= 7; ++r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 7;
        for (c = 2; c <= 8; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        r = 6;
        for (c = 8; c >= 2; --c) {
            expectedBounds.get(i).add(c, r);
        }

        // region3 =======
        i = 3;
        c = 1;
        for (r = 1; r <= 5; ++r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 5;
        for (c = 2; c <= 8; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        c = 8;
        for (r = 4; r >= 1; --r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 1;
        for (c = 7; c >= 2; --c) {
            expectedBounds.get(i).add(c, r);
        }

        // region4 =======
        i = 4;
        r = 2;
        for (c = 3; c <= 4; ++c) {
            expectedBounds.get(i).add(c, r);
        }

        // region5 =======
        i = 5;
        c = 3;
        for (r = 3; r <= 4; ++r) {
            expectedBounds.get(i).add(c, r);
        }
        c = 4;
        r = 3;
        expectedBounds.get(i).add(c, r);

        // region6 =======
        i = 6;
        r = 4;
        for (c = 4; c <= 7; ++c) {
            expectedBounds.get(i).add(c, r);
        }
        c = 7;
        for (r = 3; r >= 2; --r) {
            expectedBounds.get(i).add(c, r);
        }
        r = 2;
        for (c = 6; c >= 5; --c) {
            expectedBounds.get(i).add(c, r);
        }
        c = 5;
        r = 3;
        expectedBounds.get(i).add(c, r);

        // region7 =======
        i = 7;
        expectedBounds.get(i).add(6, 3);

        /*
        8 | 0  0  0  0  0  0  0  0  0  0
        7 | 1  2  2  2  2  2  2  2  2  0
        6 | 1  2  2  2  2  2  2  2  2  0
        5 | 1  3  3  3  3  3  3  3  3  0
        4 | 1  3  3  5  6  6  6  6  3  0
        3 | 1  3  3  5  5  6  7  6  3  0
        2 | 1  3  3  4  4  6  6  6  3  0
        1 | 1  3  3  3  3  3  3  3  3  0
        0 | 1  1  1  0  0  0  0  0  0  0
          | -  -  -  -  -  -  -  -  -  -
            0  1  2  3  4  5  6  7  8  9
        */

        for (Map.Entry<Integer, Set<Integer>> entry : filledRegionsMap.entrySet()) {
            //System.out.printf("label=%d\n", entry.getKey());
            PairIntArray boundary = PerimeterFinder3.mooreTracingWithJacob(entry.getValue(),
                    w, h);
            PairIntArray expected = expectedBounds.get(entry.getKey());
            assertEquals(expected.getN(), boundary.getN());
            for (i = 0; i < expected.getN(); ++i) {
                assertEquals(expected.getX(i), boundary.getX(i));
                assertEquals(expected.getY(i), boundary.getY(i));
            }
        }

    }


}