package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

public class MergeLabelsTest extends TestCase {

    public void testCreateLabels() {

        /*
        make a simple square of 4 colors with an edge outlining
        each square
        0 1 2 3  4
        A A A A  B B B B  0
        A A A A  B B B B  1
        A A A A  B B B B  2
        A A A A  B B B B  3
        C C C C  D D D D  4
        C C C C  D D D D
        C C C C  D D D D
        C C C C  D D D D  7

        E E E E  F F F F  8
        E E E E  F F F F
        E E E E  F F F F
        E E E E  F F F F  11
         */

        int w = 8;
        int h = 12;

        Image img = new Image(w, h);
        int[] labels = new int[img.getNPixels()];
        int r, c;
        for (int i = 0; i < img.getNPixels(); ++i) {
            r = img.getRow(i);
            c = img.getCol(i);
            if (r < 4) {
                int[] clr = ImageIOHelper.getNextRGB(0);
                if (c < 4) {
                    img.setRGB(c, r, clr[0], clr[1], clr[2]);
                    labels[i] = 0;
                } else {
                    clr[0] -= 1;
                    clr[1] -= 2;
                    img.setRGB(c, r, clr[0], clr[1], clr[2]);
                    labels[i] = 1;
                }
            } else if (r < 8) {
                if (c < 4) {
                    img.setRGB(c, r, ImageIOHelper.getNextColorRGB(2));
                    labels[i] = 2;
                } else {
                    img.setRGB(c, r, ImageIOHelper.getNextColorRGB(3));
                    labels[i] = 3;
                }
            } else {
                if (c < 4) {
                    img.setRGB(c, r, ImageIOHelper.getNextColorRGB(4));
                    labels[i] = 4;
                } else {
                    img.setRGB(c, r, ImageIOHelper.getNextColorRGB(5));
                    labels[i] = 5;
                }
            }
        }

        int[][] bounds = new int[][] {
                {3, 3}, {3, 7}, {7, 3}, {7,7}, {11, 3}, {11,7}
        };

        //
        double thresh = 3; //2.5

        int nLabels2 = MergeLabels.mergeUsingDeltaE2000(img, labels, thresh);

        // expect A and B to be merged
        int blockNumber = 0;
        for (int[] bound : bounds) {

            int r1 = bound[0];
            int r0 = bound[0] - 3;
            int c1 = bound[1];
            int c0 = bound[1] - 3;
            int label = labels[img.getInternalIndex(c0, r0)];
            for (r = r0; r <= r1; ++r) {
                for (c = c0; c <= c1; ++c) {
                    if (blockNumber == 1) {
                        assertEquals(label, labels[0]);
                    } else {
                        assertEquals(label, labels[img.getInternalIndex(c, r)]);
                    }
                }
            }
            ++blockNumber;
        }

        Image im2 = img.copyImage();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im2, labels);
        MiscDebug.writeImage(im2, "_merged_");
    }
}