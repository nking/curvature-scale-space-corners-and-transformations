package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.util.GroupAverageColors;
import algorithms.util.PairInt;
import java.util.Collection;
import java.util.Set;

/**
 * descriptions black, white, or other used in describing the template
     * shape color.  the extremes black and white can be used to limit
     * the regions created.
 * @author nichole
 */
public enum CMODE {
    WHITE, BLACK, OTHER;
    
    public static CMODE determineColorMode(ImageExt img, Set<PairInt> set) {

        GroupAverageColors clrs = new GroupAverageColors(img, set);
        
        int limit1 = 150;
        int limit2 = 55;
        if (clrs.getR() >= limit1 && clrs.getG() >= limit1 &&
            clrs.getB() >= limit1) {
            return CMODE.WHITE;
        } else if (clrs.getR() <= limit2 && clrs.getG() <= limit2 &&
            clrs.getB() <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    public static CMODE determineColorMode(ImageExt img) {

        GroupAverageColors clrs = new GroupAverageColors(img);

        int limit1 = 150;
        int limit2 = 55;
        if (clrs.getR() >= limit1 && clrs.getG() >= limit1 &&
                clrs.getB() >= limit1) {
            return CMODE.WHITE;
        } else if (clrs.getR() <= limit2 && clrs.getG() <= limit2 &&
                clrs.getB() <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    public static CMODE determinePolarThetaMode(GreyscaleImage luvTheta, 
        Collection<PairInt> points) {
    
        double avg = 0;
        for (PairInt p : points) {
            avg += luvTheta.getValue(p);
        }
        avg /= (double)points.size();
        
        int limit1 = 220;
        int limit2 = 25;
        if (avg >= limit1) {
            return CMODE.WHITE;
        } else if (avg <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    public static CMODE determinePolarThetaMode(GreyscaleImage luvTheta) {

        double avg = 0;
        for (int p = 0; p < luvTheta.getNPixels(); ++p) {
            avg += luvTheta.getValue(p);
        }
        avg /= (double)luvTheta.getNPixels();

        int limit1 = 220;
        int limit2 = 25;
        if (avg >= limit1) {
            return CMODE.WHITE;
        } else if (avg <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    static CMODE determineColorMode(GreyscaleImage rImg, 
        GreyscaleImage gImg, GreyscaleImage bImg, 
        Collection<PairInt> points) {
    
        GroupAverageColors clrs = new GroupAverageColors(
            rImg, gImg, bImg, points);
        
        int limit1 = 150;
        int limit2 = 55;
        if (clrs.getR() >= limit1 && clrs.getG() >= limit1 &&
            clrs.getB() >= limit1) {
            return CMODE.WHITE;
        } else if (clrs.getR() <= limit2 && clrs.getG() <= limit2 &&
            clrs.getB() <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

}
