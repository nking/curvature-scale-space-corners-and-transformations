package algorithms.imageProcessing.features;

import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.compGeometry.convexHull.GrahamScanPairInt;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.features.CornerRegion.CornerRegionDegneracyException;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

public class FeatureMatcher {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public FeatureMatcher() {
    }
      
    /**
     * same as ditherAndRotateForBestLocation2 except have added a filter
     * to remove matches with a cosine similarity larger than 0.95.
     * 
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param img1
     * @param img2
     * @return 
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation3(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2) {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];

        return ditherAndRotateForBestLocation3(features1, features2,
            x1, y1, x2, y2, dither, img1, img2);
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation2(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2) {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];

        return ditherAndRotateForBestLocation2(features1, features2,
            x1, y1, x2, y2, dither, img1, img2);
    }
    
    protected FeatureComparisonStat ditherAndRotateForBestLocation2(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int x2, final int y2,      
        int dither, GreyscaleImage img1, GreyscaleImage img2) {
        
        FeatureComparisonStat best = null;
        
        int rot2 = features2.calculateOrientation(x2, y2);
        
        IntensityDescriptor desc2 = features2.extractIntensity(img2, x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        int[] rotations = new int[3];
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(img1, x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(img1, y1d)) {
                    continue;
                }
                
                int rot1 = features1.calculateOrientation(x1d, y1d);
                
                // fetch rotation for this point (x1d, y1d) and try this
                // rotation and -20, -10, +10 and +20
                rotations[0] = rot1;
                rotations[1] = rot1 - 10;
                rotations[2] = rot1 + 10;
                //rotations[3] = rot1 - 20;
                //rotations[4] = rot1 + 20;
        
                for (int rotD1 : rotations) {
                    
                    if (rotD1 > 359) {
                        rotD1 = rotD1 - 360;
                    } else if (rotD1 < 0) {
                        rotD1 += 360;
                    }
                    IntensityDescriptor desc1 = features1.extractIntensity(
                        img1, x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        continue;
                    }
                    
                    FeatureComparisonStat stat = IntensityFeatures.calculateStats(
                        desc1, x1d, y1d, desc2, x2, y2);
                   
                    if (stat.getSumIntensitySqDiff() < stat.getImg2PointIntensityErr()) {
                       
                        if (best == null) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumIntensitySqDiff() > stat.getSumIntensitySqDiff()) {
                                best = stat;
                                best.setImg1PointRotInDegrees(rotD1);
                                best.setImg2PointRotInDegrees(rot2);
                            }
                        }
                    }
                }
            }
        }
        
        return best;
    }
    
    /**
     * 
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @return 
     */
    public FeatureComparisonStat matchDescriptors(
        IntensityClrFeatures features1, IntensityClrFeatures features2, 
        CornerRegion region1, CornerRegion region2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2
        ) {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        
        return matchDescriptors(features1, features2, x1, y1, x2, y2, 
            redImg1, greenImg1, blueImg1, redImg2, greenImg2, blueImg2);
    }
    
    /**
     * 
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @return 
     */
    public FeatureComparisonStat findBestMatch(
        IntensityClrFeatures features1, IntensityClrFeatures features2, 
        CornerRegion region1, CornerRegion region2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int dither
        ) {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        
        return findBestMatch(features1, features2, x1, y1, x2, y2, 
            redImg1, greenImg1, blueImg1, redImg2, greenImg2, blueImg2,
            dither);
    }
  
    /**
     * uses delta E 1994 based upon CIE LAB color space
     * @param features1
     * @param features2
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @return 
     */
    public FeatureComparisonStat matchDescriptors(
        IntensityClrFeatures features1, IntensityClrFeatures features2,
        int x1, int y1, int x2, int y2,
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1,
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2
        ) {

        int rot2;
        try {
            rot2 = features2.calculateOrientation(x2, y2);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
        
        IntensityDescriptor desc2_l = features2.extractIntensityLOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_l == null) {
            return null;
        }
        IntensityDescriptor desc2_a = features2.extractIntensityAOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_a == null) {
            return null;
        }
        IntensityDescriptor desc2_b = features2.extractIntensityBOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_b == null) {
            return null;
        }
                
        int rot1;
        try {
            rot1 = features1.calculateOrientation(x1, y1);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
                
        IntensityDescriptor desc1_l = features1.extractIntensityLOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_l == null) {
            return null;
        }
        IntensityDescriptor desc1_a = features1.extractIntensityAOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_a == null) {
            return null;
        }
        IntensityDescriptor desc1_b = features1.extractIntensityBOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_b == null) {
            return null;
        }

        FeatureComparisonStat stat_deltaE = IntensityClrFeatures.calculateStats(
            desc1_l, desc1_a, desc1_b, x1, y1, desc2_l, desc2_a, desc2_b, x2, y2);
        
        if (Float.isNaN(stat_deltaE.getSumIntensitySqDiff()) || 
            Float.isNaN(stat_deltaE.getImg2PointIntensityErr())
            ) {
            return null;
        }
        
        return stat_deltaE;
    }
    
    /**
     * uses delta E 1994 based upon CIE LAB color space
     * @param features1
     * @param features2
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @return 
     */
    public FeatureComparisonStat findBestMatch(
        IntensityClrFeatures features1, IntensityClrFeatures features2,
        int x1, int y1, int x2, int y2,
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1,
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int dither
        ) {

        FeatureComparisonStat best = null;
        
        int rot2;
        try {
            rot2 = features2.calculateOrientation(x2, y2);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
        
        IntensityDescriptor desc2_l = features2.extractIntensityLOfCIELAB(
            redImg2, greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_l == null) {
            return null;
        }
        IntensityDescriptor desc2_a = features2.extractIntensityAOfCIELAB(
            redImg2, greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_a == null) {
            return null;
        }
        IntensityDescriptor desc2_b = features2.extractIntensityBOfCIELAB(
            redImg2, greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_b == null) {
            return null;
        }
                
        int[] rotations = new int[3];
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(redImg1, x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(redImg1, y1d)) {
                    continue;
                }
                
                int rot1;
                try {
                    rot1 = features1.calculateOrientation(x1d, y1d);
                } catch (CornerRegionDegneracyException e) {
                    continue;
                }
                // fetch rotation for this point (x1d, y1d) and try this
                // rotation and -20, -10, +10 and +20
                rotations[0] = rot1;
                rotations[1] = rot1 - 10;
                rotations[2] = rot1 + 10;
                //rotations[3] = rot1 - 20;
                //rotations[4] = rot1 + 20;
        
                for (int rotD1 : rotations) {
                    
                    if (rotD1 > 359) {
                        rotD1 = rotD1 - 360;
                    } else if (rotD1 < 0) {
                        rotD1 += 360;
                    }

                    IntensityDescriptor desc1_l = features2.extractIntensityLOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1, rotD1);
                    if (desc1_l == null) {
                        continue;
                    }
                    IntensityDescriptor desc1_a = features2.extractIntensityAOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1, rotD1);
                    if (desc1_a == null) {
                        continue;
                    }
                    IntensityDescriptor desc1_b = features2.extractIntensityBOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1, rotD1);
                    if (desc1_b == null) {
                        continue;
                    }

                    FeatureComparisonStat stat = IntensityClrFeatures.calculateStats(
                        desc1_l, desc1_a, desc1_b, x1d, y1d, desc2_l, desc2_a, desc2_b, x2, y2);

                    if (Float.isNaN(stat.getSumIntensitySqDiff())
                        || Float.isNaN(stat.getImg2PointIntensityErr())) {
                        continue;
                    }
                    
                    if (stat.getSumIntensitySqDiff() < stat.getImg2PointIntensityErr()) {                       
                        if (best == null) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumIntensitySqDiff() > stat.getSumIntensitySqDiff()) {
                                best = stat;
                                best.setImg1PointRotInDegrees(rotD1);
                                best.setImg2PointRotInDegrees(rot2);
                            }
                        }
                    }
                }
            }
        }
        
        return best;
    }
    
    public FeatureComparisonStat findBestMatchO123(
        IntensityClrFeatures features1, IntensityClrFeatures features2,
        int x1, int y1, int x2, int y2,
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1,
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2
        ) {

        int rot2;
        try {
            rot2 = features2.calculateOrientation(x2, y2);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
        
        /*
        extractIntensityO3(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation)
        */
        
        IntensityDescriptor desc2_o1 = features2.extractIntensityO1(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_o1 == null) {
            return null;
        }
        IntensityDescriptor desc2_o2 = features2.extractIntensityO2(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_o2 == null) {
            return null;
        }
        IntensityDescriptor desc2_o3 = features2.extractIntensityO3(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_o3 == null) {
            return null;
        }
                
        int rot1;
        try {
            rot1 = features1.calculateOrientation(x1, y1);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
                
        IntensityDescriptor desc1_o1 = features1.extractIntensityO1(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_o1 == null) {
            return null;
        }
        IntensityDescriptor desc1_o2 = features1.extractIntensityO2(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_o2 == null) {
            return null;
        }
        IntensityDescriptor desc1_o3 = features1.extractIntensityO1(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_o3 == null) {
            return null;
        }

        FeatureComparisonStat stat_o1 = IntensityFeatures.calculateStats(
            desc1_o1, x1, y1, desc2_o1, x2, y2);
        
        if (stat_o1.getSumIntensitySqDiff() >= stat_o1.getImg2PointIntensityErr()) {
            return null;
        }
        
        FeatureComparisonStat stat_o2 = IntensityFeatures.calculateStats(
            desc1_o2, x1, y1, desc2_o2, x2, y2);
        
        if (stat_o2.getSumIntensitySqDiff() >= stat_o2.getImg2PointIntensityErr()) {
            return null;
        }
        
        FeatureComparisonStat stat_o3 = IntensityFeatures.calculateStats(
            desc1_o3, x1, y1, desc2_o3, x2, y2);
        
        if (stat_o3.getSumIntensitySqDiff() >= stat_o3.getImg2PointIntensityErr()) {
            return null;
        }
        
        /*
        double ssdIntensity = 
            Math.sqrt(stat_o1.getSumIntensitySqDiff()) +
            Math.sqrt(stat_o2.getSumIntensitySqDiff()) +
            Math.sqrt(stat_o3.getSumIntensitySqDiff());
        
        double err2SqIntensity = 
            Math.sqrt(stat_o1.getImg2PointIntensityErr()) +
            Math.sqrt(stat_o2.getImg2PointIntensityErr()) +
            Math.sqrt(stat_o3.getImg2PointIntensityErr());
        */
        double ssdIntensity = 
            Math.sqrt(stat_o1.getSumIntensitySqDiff()) +
            Math.sqrt(stat_o2.getSumIntensitySqDiff());
        
        double err2SqIntensity = 
            Math.sqrt(stat_o1.getImg2PointIntensityErr()) +
            Math.sqrt(stat_o2.getImg2PointIntensityErr());
        
        log.info(String.format("(%d,%d), (%d,%d)  %.1f  %.1f  %.1f => %.1f (%.1f)",
            x1, y1, x2, y2, 
            stat_o1.getSumIntensitySqDiff(), stat_o2.getSumIntensitySqDiff(),
            stat_o3.getSumIntensitySqDiff(), ssdIntensity, err2SqIntensity));
        
        FeatureComparisonStat best = new FeatureComparisonStat();
        best.setImg1Point(new PairInt(x1, y1));
        best.setImg2Point(new PairInt(x2, y2));
        best.setSumIntensitySqDiff((float) ssdIntensity);
        best.setImg2PointIntensityErr((float) err2SqIntensity);

        return best;
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation3(
        IntensityClrFeatures features1,
        IntensityClrFeatures features2,
        final int x1, final int y1, final int x2, final int y2,      
        int dither, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1,
        GreyscaleImage blueImg1,
        GreyscaleImage redImg2, GreyscaleImage greenImg2,
        GreyscaleImage blueImg2
    ) {
        
        FeatureComparisonStat best = null;
        
        int rot2;
        try {
            rot2 = features2.calculateOrientation(x2, y2);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
        
        // make L, A, B descriptors
        IntensityDescriptor desc2_l = features2.extractIntensityLOfCIELAB(redImg2,
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_l == null) {
            return null;
        }
        IntensityDescriptor desc2_a = features2.extractIntensityAOfCIELAB(redImg2,
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_a == null) {
            return null;
        }
        IntensityDescriptor desc2_b = features2.extractIntensityBOfCIELAB(redImg2,
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_b == null) {
            return null;
        }
        
        int[] rotations = new int[3];
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(redImg1, x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(redImg1, y1d)) {
                    continue;
                }
                
                int rot1;
                try {
                    rot1 = features1.calculateOrientation(x1d, y1d);
                } catch (CornerRegionDegneracyException e) {
                    continue;
                }
                
                // fetch rotation for this point (x1d, y1d) and try this
                // rotation and -20, -10, +10 and +20
                rotations[0] = rot1;
                rotations[1] = rot1 - 10;
                rotations[2] = rot1 + 10;
                //rotations[3] = rot1 - 20;
                //rotations[4] = rot1 + 20;
        
                for (int rotD1 : rotations) {
                    
                    if (rotD1 > 359) {
                        rotD1 = rotD1 - 360;
                    } else if (rotD1 < 0) {
                        rotD1 += 360;
                    }
                    
                    IntensityDescriptor desc1_l = features1.extractIntensityLOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1d, rot1);
                    if (desc1_l == null) {
                        return null;
                    }
                    IntensityDescriptor desc1_a = features1.extractIntensityAOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1d, rot1);
                    if (desc1_a == null) {
                        return null;
                    }
                    IntensityDescriptor desc1_b = features1.extractIntensityBOfCIELAB(redImg1,
                        greenImg1, blueImg1, x1d, y1d, rot1);
                    if (desc1_b == null) {
                        return null;
                    }
                
                   FeatureComparisonStat stat_deltaE = 
                       IntensityClrFeatures.calculateStats(
                       desc1_l, desc1_a, desc1_b, x1, y1, 
                       desc2_l, desc2_a, desc2_b, x2, y2);

                    if (Float.isNaN(stat_deltaE.getSumIntensitySqDiff())
                        || Float.isNaN(stat_deltaE.getImg2PointIntensityErr())) {
                        return null;
                    }
                                        
                    if (stat_deltaE.getSumIntensitySqDiff() 
                        < stat_deltaE.getImg2PointIntensityErr()) {
                       
                        if (best == null) {
                            best = stat_deltaE;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumIntensitySqDiff() 
                                > stat_deltaE.getSumIntensitySqDiff()) {
                                best = stat_deltaE;
                                best.setImg1PointRotInDegrees(rotD1);
                                best.setImg2PointRotInDegrees(rot2);
                            }
                        }
                    }
                }
            }
        }
        
        return best;
    }
    
    protected FeatureComparisonStat ditherAndRotateForBestLocation3(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int x2, final int y2,      
        int dither, GreyscaleImage img1, GreyscaleImage img2) {
        
        FeatureComparisonStat best = null;
        
        int rot2 = features2.calculateOrientation(x2, y2);
       
        IntensityDescriptor desc2 = features2.extractIntensity(img2, x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        int[] rotations = new int[3];
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(img1, x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(img1, y1d)) {
                    continue;
                }
                
                int rot1 = features1.calculateOrientation(x1d, y1d);
                
                // fetch rotation for this point (x1d, y1d) and try this
                // rotation and -20, -10, +10 and +20
                rotations[0] = rot1;
                rotations[1] = rot1 - 10;
                rotations[2] = rot1 + 10;
                //rotations[3] = rot1 - 20;
                //rotations[4] = rot1 + 20;
        
                for (int rotD1 : rotations) {
                    
                    if (rotD1 > 359) {
                        rotD1 = rotD1 - 360;
                    } else if (rotD1 < 0) {
                        rotD1 += 360;
                    }
                    IntensityDescriptor desc1 = features1.extractIntensity(
                        img1, x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        continue;
                    }
                    
                    FeatureComparisonStat stat = IntensityFeatures.calculateStats(
                        desc1, x1d, y1d, desc2, x2, y2);
                    
                    float cSim = desc1.calculateCosineSimilarity(desc2);
                    
                    if (cSim < 0.95) {
                        continue;
                    }
                   
                    if (stat.getSumIntensitySqDiff() < stat.getImg2PointIntensityErr()) {
                       
                        if (best == null) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumIntensitySqDiff() > stat.getSumIntensitySqDiff()) {
                                best = stat;
                                best.setImg1PointRotInDegrees(rotD1);
                                best.setImg2PointRotInDegrees(rot2);
                            }
                        }
                    }
                }
            }
        }
        
        return best;
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation2(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        final int expectedRotationInDegrees, final int rotationTol,
        GreyscaleImage img1, GreyscaleImage img2) {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        
        return ditherAndRotateForBestLocation2(features1, features2, 
            x1, y1, x2, y2, dither,
            expectedRotationInDegrees, rotationTol, img1, img2);
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation2(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int x2, final int y2,      
        int dither, int expectedRotationInDegrees,
        final int rotationTol, GreyscaleImage img1, GreyscaleImage img2) {
        
        FeatureComparisonStat best = null;
        
        int rot2 = features2.calculateOrientation(x2, y2);
        
        IntensityDescriptor desc2 = features2.extractIntensity(img2, x2, y2, rot2);
     
        if (desc2 == null) {
            return null;
        }
        
        // because anglediff compares closest angles:
        if (expectedRotationInDegrees > 180.) {
            expectedRotationInDegrees = 360 - expectedRotationInDegrees;
        }
        
        int[] rotations = new int[3];
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(img1, x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(img1, y1d)) {
                    continue;
                }
                
                int rot1 = features1.calculateOrientation(x1d, y1d);
                
                // fetch rotation for this point (x1d, y1d) and try this
                // rotation and -20, -10, +10 and +20
                rotations[0] = rot1;
                rotations[1] = rot1 - 10;
                rotations[2] = rot1 + 10;
                //rotations[3] = rot1 - 20;
                //rotations[4] = rot1 + 20;
                        
                for (int rotD1 : rotations) {
                    if (rotD1 > 359) {
                        rotD1 = rotD1 - 360;
                    } else if (rotD1 < 0) {
                        rotD1 += 360;
                    }
            
                    // only try rotations within expected rotation limits
                    float rotDiffs = AngleUtil.getAngleDifference(rotD1, rot2);
                    if (rotDiffs < 0) {
                        rotDiffs *= -1;
                    }
 //NOTE: change to rotDiffs to use -1* instead of 360+=             
                    float rotDiff = AngleUtil.getAngleDifference(
                        expectedRotationInDegrees, rotDiffs);
                    if (Math.abs(rotDiff) > rotationTol) {
                        continue;
                    }
            
                    IntensityDescriptor desc1 = features1.extractIntensity(img1, 
                        x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        continue;
                    }
                    
                    FeatureComparisonStat stat = IntensityFeatures.calculateStats(
                        desc1, x1d, y1d, desc2, x2, y2);
                   
                    if (stat.getSumIntensitySqDiff() > stat.getImg2PointIntensityErr()) {
                        continue;
                    }
                    if ((best == null) || (best.getSumIntensitySqDiff() > stat.getSumIntensitySqDiff())) {
                        best = stat;
                        best.setImg1PointRotInDegrees(rotD1);
                        best.setImg2PointRotInDegrees(rot2);
                    }
                }
            }
        }
        
        return best;
    }
    
    public CorrespondenceList findSimilarFeatures(GreyscaleImage gsImg1, 
        CornerRegion[] cr1, GreyscaleImage gsImg2, 
        CornerRegion[] cr2, TransformationParameters params, float scaleTol, 
        float rotationInRadiansTol, int transXYTol, int dither,
        RotatedOffsets rotatedOffsets) {
      
        List<FeatureComparisonStat> stats = 
            findSimilarFeaturesAsStats(gsImg1, cr1, gsImg2, cr2, params, 
                scaleTol, rotationInRadiansTol, transXYTol, dither,
                rotatedOffsets);
        
        if (stats == null) {
            return null;
        }
            
        List<PairInt> matched1 = new ArrayList<PairInt>();
        List<PairInt> matched2 = new ArrayList<PairInt>();

        for (FeatureComparisonStat stat : stats) {
            matched1.add(stat.getImg1Point());
            matched2.add(stat.getImg2Point());
        }
        
        int rangeRotation = Math.round(params.getStandardDeviations()[1]);
        int rangeTranslationX = Math.round(params.getStandardDeviations()[2]);
        int rangeTranslationY = Math.round(params.getStandardDeviations()[3]);

        CorrespondenceList cl = new CorrespondenceList(params.getScale(),
            Math.round(params.getRotationInDegrees()),
            Math.round(params.getTranslationX()), Math.round(params.getTranslationY()),
            rangeRotation, rangeTranslationX, rangeTranslationY,
            matched1, matched2);

        return cl;
    }    
        
    public List<FeatureComparisonStat> findSimilarFeaturesAsStats(GreyscaleImage gsImg1, 
        CornerRegion[] cr1, GreyscaleImage gsImg2, 
        CornerRegion[] cr2, TransformationParameters params, float scaleTol, 
        float rotationInRadiansTol, int transXYTol, int dither,
        RotatedOffsets rotatedOffsets) {
        
        final int blockHalfWidth = 5;
        final boolean useNormalizedIntensities = true;
        
        IntensityFeatures features1 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features1.calculateGradientWithGreyscale(gsImg1);
        IntensityFeatures features2 = new IntensityFeatures(blockHalfWidth,
            useNormalizedIntensities, rotatedOffsets);
        features2.calculateGradientWithGreyscale(gsImg2);
        
        List<FeatureComparisonStat> stats = findSimilarFeaturesAsStats(gsImg1, 
            cr1, gsImg2, cr2, features1, features2, params, scaleTol, 
            rotationInRadiansTol, transXYTol, dither, rotatedOffsets);
        
        return stats;
    }    
    
    public List<FeatureComparisonStat> findSimilarFeaturesAsStats(
        GreyscaleImage gsImg1, CornerRegion[] cr1s, CornerRegion[] cr1Trs, 
        GreyscaleImage gsImg2, CornerRegion[] cr2s, 
        IntensityFeatures features1, IntensityFeatures features2, 
        TransformationParameters parameters, int dither, int transXYTol,  
        RotatedOffsets rotatedOffsets) {
        
        //for each combination of cr1 and cr2, find best stat if any between filters
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        for (int idx1 = 0; idx1 < cr1s.length; ++idx1) {
            
            CornerRegion cr1 = cr1s[idx1];
            CornerRegion cr1Tr = cr1Trs[idx1];
            
            int x1 = cr1.getX()[cr1.getKMaxIdx()];
            int y1 = cr1.getY()[cr1.getKMaxIdx()];
            
            FeatureComparisonStat bestStat = null;
            
            for (int idx2 = 0; idx2 < cr2s.length; ++idx2) {
                
                CornerRegion cr2 = cr2s[idx2];
                
                int x2 = cr2.getX()[cr2.getKMaxIdx()];
                int y2 = cr2.getY()[cr2.getKMaxIdx()];
                
                double diffX = cr1Tr.getX()[cr1Tr.getKMaxIdx()] - x2;
                double diffY = cr1Tr.getY()[cr1Tr.getKMaxIdx()] - y2;
                                
                if ((diffX > transXYTol) || (diffY > transXYTol)) {
                    continue;
                }
                                
                FeatureComparisonStat stat = ditherAndRotateForBestLocation2(
                    features1, features2, x1, y1, x2, y2, dither, gsImg1, gsImg2);
                
                if ((stat == null) || 
                    (stat.getSumIntensitySqDiff() > stat.getImg2PointIntensityErr())) {
                    continue;
                }
                
                if ((bestStat == null) ||
                    stat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                    
                    bestStat = stat;
                }
            }
            
            if (bestStat != null) {
                stats.add(bestStat);
            }
        }
        
        return stats;
    }
    
    public List<FeatureComparisonStat> findSimilarFeaturesAsStats(GreyscaleImage gsImg1, 
        CornerRegion[] cr1, GreyscaleImage gsImg2, 
        CornerRegion[] cr2, IntensityFeatures features1, IntensityFeatures features2,
        TransformationParameters params, float scaleTol, 
        float rotationInRadiansTol, int transXYTol, int dither,
        RotatedOffsets rotatedOffsets) {
        
        List<CornerRegion> filteredTransformedC1 = new ArrayList<CornerRegion>();
        List<CornerRegion> filteredC1 = new ArrayList<CornerRegion>();
        List<CornerRegion> filteredC2 = new ArrayList<CornerRegion>();
        
        filterForIntersection4(params, transXYTol, 
            cr1, cr2, filteredTransformedC1, filteredC1, filteredC2,
            gsImg1.getWidth(), gsImg1.getHeight(), gsImg2.getWidth(), gsImg2.getHeight());

        if (true) {
            try {
                MiscDebug.writeImage(filteredC1, gsImg1.copyToColorGreyscale(),
                    "filtered_1_corners_");
                MiscDebug.writeImage(filteredC2, gsImg2.copyToColorGreyscale(), 
                    "filtered_2_corners_");
                MiscDebug.writeImage(filteredTransformedC1, gsImg2.copyToColorGreyscale(), 
                    "filtered_1_trans_corners_");
            } catch (IOException ex) {
                Logger.getLogger(FeatureMatcher.class.getName()).log(
                    Level.SEVERE, null, ex);
            }
        }
        
        int n1 = filteredC1.size();
        int n2 = filteredC2.size();
        
        /*
        when transformation params are known ahead of time:
        cr1 can be transformed into crTr1 (including the internal points).
        then the matching is faster than n1 * n2 because can discard some 
        possible matches immediately.
        bipartite matching when points are present within tolerance.
        
        bipartite is n^3 but the n is < n1.
        */
               
        Map<PairInt, FeatureComparisonStat> statMap = null;
        float[][] cost = null;
                
        final boolean useBipartite = false;//(nMaxMatchable < 251);
        
        Map<Integer, Integer> index1Map = null;
        Map<Integer, Set<Integer>> index2Map = null;
        Map<Integer, FeatureComparisonStat> index1StatMap = null;
        
        if (useBipartite) {
            cost = new float[n1][n2];
            statMap = new HashMap<PairInt, FeatureComparisonStat>();
        } else {
            index1Map = new HashMap<Integer, Integer>();
            index2Map = new HashMap<Integer, Set<Integer>>();
            index1StatMap = new HashMap<Integer, FeatureComparisonStat>();
        }
        
        int rotationInDegrees = Math.round(params.getRotationInDegrees());
        int rotationToleranceInDegrees = (int)Math.round(rotationInRadiansTol * 180/Math.PI);
                
        Transformer transformer = new Transformer();
        
        int count = 0;
        
        for (int i1 = 0; i1 < n1; ++i1) {
            
            if (useBipartite) {
                cost[i1] = new float[n2];
                Arrays.fill(cost[i1], Float.MAX_VALUE);
            }
            
            CornerRegion c1Tr = filteredTransformedC1.get(i1);
            CornerRegion c1 = filteredC1.get(i1);
            
            int x1Tr = c1Tr.getX()[c1Tr.getKMaxIdx()];
            int y1Tr = c1Tr.getY()[c1Tr.getKMaxIdx()];
            
            int x1 = c1.getX()[c1.getKMaxIdx()];
            int y1 = c1.getY()[c1.getKMaxIdx()];
            
            // additional tolerance due to rotation error of 10 degrees
            double[] errTr = transformer.applyTransformation(params.getScale(),
                params.getRotationInRadians() + 0.1745,
                params.getOriginX(), params.getOriginY(),
                params.getTranslationX(), params.getTranslationY(), x1, y1);
            
            double xTolAdd = Math.abs(errTr[0] - x1Tr);
            double yTolAdd = Math.abs(errTr[1] - y1Tr);
            
            double bestCost = Double.MAX_VALUE;
            int bestIdx2 = -1;
            FeatureComparisonStat bestStat = null;
            
//TODO: replace w/ a nearest neighbors structure ***** <=====
            
            for (int i2 = 0; i2 < n2; ++i2) {
                
                CornerRegion c2 = filteredC2.get(i2);
                
                int x2 = c2.getX()[c2.getKMaxIdx()];
                int y2 = c2.getY()[c2.getKMaxIdx()];
                                
                int diffX = Math.abs(x1Tr - x2);
                int diffY = Math.abs(y1Tr - y2);
                if (diffX > (transXYTol + xTolAdd) || diffY > (transXYTol + yTolAdd)) {
                    continue;
                }
                // use the untransformed cr1 to be able to filter by rotation
                FeatureComparisonStat stat = ditherAndRotateForBestLocation2(
                    features1, features2, c1, c2, dither,
                    rotationInDegrees, rotationToleranceInDegrees,
                    gsImg1, gsImg2);

                if (stat != null && 
                    (stat.getSumIntensitySqDiff() < stat.getImg2PointIntensityErr())) {

                    if (useBipartite) {
                        cost[i1][i2] = stat.getSumIntensitySqDiff();
                        PairInt p = new PairInt(i1, i2);
                        statMap.put(p, stat);
                    } else {
                        if ((bestIdx2 == -1) || (bestCost > stat.getSumIntensitySqDiff())) {
                            bestIdx2 = i2;
                            bestCost = stat.getSumIntensitySqDiff();
                            bestStat = stat;
                            bestStat.setIndex1(i1);
                            bestStat.setIndex2(i2);
                        }
                    }

                    count++;
                }
            }
            
            if (!useBipartite && (bestStat != null)) {
                Integer key1 = Integer.valueOf(i1);
                Integer key2 = Integer.valueOf(bestIdx2);
                index1Map.put(key1, key2);
                Set<Integer> set = index2Map.get(key2);
                if (set == null) {
                    set = new HashSet<Integer>();
                    index2Map.put(key2, set);
                }
                set.add(key1);
                
                index1StatMap.put(key1, bestStat);
            }
        }
        
        if (useBipartite) {
            return useBipartiteMatchingAsStats(cost, statMap, params);
        }
        
        // resolve any double matchings, but discard the higher cost matches
        //   from conflicted matches rather than re-trying a solution for them
        
        Set<Integer> resolved = new HashSet<Integer>();
        for (Entry<Integer, Set<Integer>> entry : index2Map.entrySet()) {
            if (resolved.contains(entry.getKey())) {
                continue;
            }
            Set<Integer> set = entry.getValue();
            if (set.size() > 1) {
                double bestCost = Double.MAX_VALUE;
                Integer bestIndex1 = -1;
                for (Integer index1 : set) {
                    FeatureComparisonStat fcs = index1StatMap.get(index1);
                    assert(fcs != null);
                    double cost2 = fcs.getSumIntensitySqDiff();
                    if (cost2 < bestCost) {
                        bestCost = cost2;
                        bestIndex1 = index1;
                    }
                    resolved.add(index1);
                }
                for (Integer index1 : set) {
                    if (index1.equals(bestIndex1)) {
                        continue;
                    }
                    index1Map.remove(index1);
                    index1StatMap.remove(index1);
                }
            }
        }
        
        int nc = index1StatMap.size();
        
        if (nc < 7) {
            return null;
        }
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        float[] weights = new float[nc];
        double sumW = 0;

        nc = 0;
        for (Entry<Integer, FeatureComparisonStat> entry : index1StatMap.entrySet()) {
            FeatureComparisonStat fcs = entry.getValue();
            stats.add(fcs);
            weights[nc] = fcs.getSumIntensitySqDiff();
            sumW += weights[nc];
            nc++;
        }
        if (sumW > 0) {
            double tot = 0;
            for (int i = 0; i < nc; ++i) {
                double div = (sumW - weights[i]) / ((nc - 1) * sumW);
                weights[i] = (float) div;
                tot += div;
            }
            assert (Math.abs(tot - 1.) < 0.03);
        } else {
            float a = 1.f / (float) nc;
            Arrays.fill(weights, a);
        }

        return stats;
    }
    
    public static float[] calculateThetaDiff(List<FeatureComparisonStat> compStats) {
        
        if (compStats == null || compStats.isEmpty()) {
            return null;
        }
        
        float[] values = new float[compStats.size()];
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            float diff = AngleUtil.getAngleDifference(
                stat.getImg1PointRotInDegrees(), stat.getImg2PointRotInDegrees());
            
            values[i] = diff;
        }
        
        return values;
    }

    public static void filterForIntersection3(TransformationParameters params, 
        int transXYTol, CornerRegion[] c1, CornerRegion[] c2, 
        List<CornerRegion> outFilteredTransformedC1, 
        List<CornerRegion> outFilteredC1, 
        List<CornerRegion> outFilteredC2,
        int img1Width, int img1Height, int img2Width, int img2Height) {
        
        /*
        transform corners1 to image2 reference frame and trim any points
           in it that are out of the image2 frame.
        then make a
        */
        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        Transformer transformer = new Transformer();
        
        TransformationParameters revParams = tc.swapReferenceFrames(params);
        
        for (int i = 0; i < c1.length; ++i) {
                        
            CornerRegion ctr = transformer.applyTransformation(params, c1[i]);
                
            int xTr = ctr.getX()[ctr.getKMaxIdx()];
            int yTr = ctr.getY()[ctr.getKMaxIdx()];

            if ((xTr < 0) || (xTr > (img2Width - 1))) {
                continue;
            }
            if ((yTr < 0) || (yTr > (img2Height - 1))) {
                continue;
            }

            outFilteredTransformedC1.add(ctr);
            outFilteredC1.add(c1[i]);
            
        }

        for (int i = 0; i < c2.length; ++i) {
                        
            int x = c2[i].getX()[c2[i].getKMaxIdx()];
            int y = c2[i].getY()[c2[i].getKMaxIdx()];
            
            double[] xyTr = transformer.applyTransformation(revParams, x, y);
                            
            if ((xyTr[0] < 0) || (xyTr[0] > (img1Width - 1))) {
                continue;
            }
            if ((xyTr[1] < 0) || (xyTr[1] > (img1Height - 1))) {
                continue;
            }

            outFilteredC2.add(c2[i]);
        }
    }

    public static List<Integer> removeDiscrepantThetaDiff(
        List<FeatureComparisonStat> compStats) {
        
        if (compStats == null || compStats.isEmpty()) {
            return null;
        }
        
        float[] values = calculateThetaDiff(compStats);
        for (int i = 0; i < values.length; ++i) {
            if (values[i] < 0) {
                values[i] += 360;
            }
        }
        
        // 20 degree wide bins
        HistogramHolder hist = Histogram.createSimpleHistogram(20.f, values, 
            Errors.populateYErrorsBySqrt(values));
        
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        if ((yMaxIdx > -1) && (hist.getYHist()[yMaxIdx] == 1)) {
            hist = Histogram.createSimpleHistogram(40.f, values, 
                Errors.populateYErrorsBySqrt(values));
            yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        }
        
        float thetaDiff;
        if (yMaxIdx == -1) {
            float[] thetaDiffMeanStDev = MiscMath.getAvgAndStDev(values);
            thetaDiff = thetaDiffMeanStDev[0];
        } else {
            thetaDiff = hist.getXHist()[yMaxIdx];
        }
        
        //TODO: consider a bin larger than 20 degrees... 25
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < values.length; ++i) {
            float diffRot = AngleUtil.getAngleDifference(values[i], thetaDiff);
            if (diffRot > 20) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        for (int i = remove.size() - 1; i > -1; --i) {
            int idx = remove.get(i);
            compStats.remove(idx);
        }
        
        return remove;
    }
    
    public static List<Integer> removeDiscrepantThetaDiff(
        List<FeatureComparisonStat> compStats, float rotationInDegrees) {
        
        if (compStats == null || compStats.isEmpty()) {
            return null;
        }
        
        List<Integer> remove = new ArrayList<Integer>();
        
        float[] values = calculateThetaDiff(compStats);
        for (int i = 0; i < values.length; ++i) {
            
            if (values[i] < 0) {
                values[i] += 360;
            }
            
            int diffRot = Math.round(AngleUtil.getAngleDifference(values[i], 
                rotationInDegrees));
            if (diffRot > 20) {
                if ((diffRot < 30) && compStats.get(i).getSumIntensitySqDiff() < 100) {
                    continue;
                }
                remove.add(Integer.valueOf(i));
            }
        }
        
        for (int i = remove.size() - 1; i > -1; --i) {
            int idx = remove.get(i);
            compStats.remove(idx);
        }
        
        return remove;
    }
    
    public static float[] calcIntensitySSDMeanAndStDev(List<FeatureComparisonStat> compStats) {
        
        int n = compStats.size();
        
        float[] ssds = new float[n];
        
        for (int i = 0; i < n; ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            ssds[i] = stat.getSumIntensitySqDiff();
        }
        
        float[] meanStDv = MiscMath.getAvgAndStDev(ssds);
        
        return meanStDv;
    }
    
    public static List<Integer> removeIntensityOutliers(List<FeatureComparisonStat> 
        compStats) {
        
        if (compStats.size() < 3) {
            return new ArrayList<Integer>();
        }
        
        float sigmaFactor = 1.25f;
        
        return removeIntensityOutliers(compStats, sigmaFactor);
    }

    public static List<Integer> removeIntensityOutliers(List<FeatureComparisonStat> 
        compStats, float sigmaFactor) {
        
        if (compStats.size() < 3) {
            return new ArrayList<Integer>();
        }
        
        int n = compStats.size();
        
        float[] meanStDv = calcIntensitySSDMeanAndStDev(compStats);
        
        List<Integer> rm = new ArrayList<Integer>();
        
        for (int i = 0; i < n; ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            float diff = stat.getSumIntensitySqDiff() - meanStDv[0];
            
            if (diff > (sigmaFactor * meanStDv[1])) {
                rm.add(Integer.valueOf(i));
            }
        }
        
        for (int i = rm.size() - 1; i > -1; --i) {
            
            int idx = rm.get(i).intValue();
            
            compStats.remove(idx);
        }
        
        return rm;
    }

    private List<FeatureComparisonStat> useBipartiteMatchingAsStats(float[][] cost,
        Map<PairInt, FeatureComparisonStat> statMap, 
        TransformationParameters params) {

        if (cost == null) {
            return null;
        }
        
        boolean transposed = false;
        if (cost.length > cost[0].length) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }

        // one pass thru to count for array sizes
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);

        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            PairInt pI = new PairInt(idx1, idx2);
            if (!statMap.containsKey(pI)) {
                continue;
            }
        }

        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        int nc = 0;

        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairInt pI = new PairInt(idx1, idx2);
            FeatureComparisonStat stat = statMap.get(pI);
            if (stat == null) {
                continue;
            }
            
            stats.add(stat);
            
            nc++;
        }

        if (stats.size() < 7) {
            return null;
        }
        
        return stats;
    }
    
    public static void filterForIntersection4(TransformationParameters params, 
        int transXYTol, CornerRegion[] c1, CornerRegion[] c2, 
        List<CornerRegion> outFilteredTransformedC1, 
        List<CornerRegion> outFilteredC1, 
        List<CornerRegion> outFilteredC2,
        int img1Width, int img1Height, int img2Width, int img2Height) {
        
        if (c1.length < 3 || c2.length < 3) {
            return;
        }
        
        /*
        make convex hull of corners1 and transform the hull to corners2
           where will only keep the corners2 that are within the transformed
           hull.
        then will do the same with corners2 to corners1
        */
        
        Transformer transformer = new Transformer();
        PointInPolygon pip = new PointInPolygon();
        
        try {
            
            PairInt[] corners1 = Misc.convert(c1);
            GrahamScanPairInt<PairInt> scan = new GrahamScanPairInt<PairInt>();
            scan.computeHull(corners1);
            List<PairInt> hull = scan.getHull();
            
            float[] xTr1 = new float[hull.size()];
            float[] yTr1 = new float[hull.size()];
            for (int i = 0; i < hull.size(); ++i) {
                PairInt xy = hull.get(i);
                double[] xyTr = transformer.applyTransformation(params, xy.getX(), xy.getY());
                xTr1[i] = (float)xyTr[0];
                yTr1[i] = (float)xyTr[1];
                if (xTr1[i] < 0) {
                    xTr1[i] = 0;
                } else if (xTr1[i] > (img2Width - 1)) {
                    xTr1[i] = img2Width - 1;
                }
                if (yTr1[i] < 0) {
                    yTr1[i] = 0;
                } else if (yTr1[i] > (img2Height - 1)) {
                    yTr1[i] = img2Height - 1;
                }
            }
            
            for (int i = 0; i < c2.length; ++i) {
                CornerRegion cr = c2[i];
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                
                if (pip.isInSimpleCurve(x, y, xTr1, yTr1, yTr1.length)) {
                    outFilteredC2.add(cr);
                }
            }
            
            MatchedPointsTransformationCalculator tc
                = new MatchedPointsTransformationCalculator();

            TransformationParameters revParams = tc.swapReferenceFrames(params);
        
            if (outFilteredC2.size() < 3) {
                return;
            }
            
            PairInt[] corners2 = Misc.convert(outFilteredC2);
            scan = new GrahamScanPairInt<PairInt>();
            scan.computeHull(corners2);
            hull = scan.getHull();
            
            float[] xTr2 = new float[hull.size()];
            float[] yTr2 = new float[hull.size()];
            for (int i = 0; i < hull.size(); ++i) {
                PairInt xy = hull.get(i);
                double[] xyTr = transformer.applyTransformation(revParams, 
                    xy.getX(), xy.getY());
                xTr2[i] = (float)xyTr[0];
                yTr2[i] = (float)xyTr[1];
                if (xTr2[i] < 0) {
                    xTr2[i] = 0;
                } else if (xTr2[i] > (img1Width - 1)) {
                    xTr2[i] = img1Width - 1;
                }
                if (yTr2[i] < 0) {
                    yTr2[i] = 0;
                } else if (yTr2[i] > (img1Height - 1)) {
                    yTr2[i] = img1Height - 1;
                }
            }
            
            for (int i = 0; i < c1.length; ++i) {
                CornerRegion cr = c1[i];
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                
                if (pip.isInSimpleCurve(x, y, xTr2, yTr2, yTr2.length)) {
                    CornerRegion ctr = transformer.applyTransformation(params, c1[i]);
                    outFilteredTransformedC1.add(ctr);
                    outFilteredC1.add(c1[i]);
                }
            }
            
        } catch (GrahamScanTooFewPointsException e) {
        }
    }

    public List<FeatureComparisonStat> reviseStatsForFullImages(
        GreyscaleImage img1, GreyscaleImage img2,
        FeatureMatcherSettings settings,
        TransformationParameters params,
        List<FeatureComparisonStat> stats, 
        int prevBinFactor1, int prevBinFactor2,
        RotatedOffsets rotatedOffsets) {

        log.info("refine stats for full image reference frames");

        List<FeatureComparisonStat> revised = new ArrayList<FeatureComparisonStat>();

        FeatureMatcher featureMatcher = new FeatureMatcher();

        IntensityFeatures features1 = new IntensityFeatures(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features1.calculateGradientWithGreyscale(img1);

        IntensityFeatures features2 = new IntensityFeatures(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features2.calculateGradientWithGreyscale(img2);

        int dither2 = 1 * (Math.max(prevBinFactor1, prevBinFactor2));
        if (params.getStandardDeviations()[2] > 25 || params.getStandardDeviations()[3] > 25) {
            if (dither2 < 3) {
                dither2 = 3;
            }
        }

        int rotD = Math.round(params.getRotationInDegrees());

        final int rotationTolerance = 20;

        for (int i = 0; i < stats.size(); ++i) {

            FeatureComparisonStat stat = stats.get(i);

            int x1 = stat.getImg1Point().getX() * stat.getBinFactor1();
            int y1 = stat.getImg1Point().getY() * stat.getBinFactor1();
            int x2 = stat.getImg2Point().getX() * stat.getBinFactor2();
            int y2 = stat.getImg2Point().getY() * stat.getBinFactor2();

            // have to discard the best angles found in stat and derive new
            // for these higher resolution images
            FeatureComparisonStat compStat =
                featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, x1, y1, x2, y2, dither2, rotD,
                    rotationTolerance, img1, img2);

            if (compStat == null ||
                (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())) {
                continue;
            }

            revised.add(compStat);
        }
        
        return revised;
    }
    
    public FeatureComparisonStat matchHalfDescriptors(
        IntensityClrFeatures features1, IntensityClrFeatures features2, 
        KeyPointsAndBounds keyPointsAndBounds1, int bmaIndex1,
        KeyPointsAndBounds keyPointsAndBounds2, int bmaIndex2,
        PairInt keyPoint1, PairInt keyPoint2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, 
        GreyscaleImage blueImg1, GreyscaleImage redImg2, 
        GreyscaleImage greenImg2, GreyscaleImage blueImg2) {
        
        int x1 = keyPoint1.getX();
        int y1 = keyPoint1.getY();
        int x2 = keyPoint2.getX();
        int y2 = keyPoint2.getY();
        
        return matchHalfDescriptors(features1, features2, 
            keyPointsAndBounds1, bmaIndex1, 
            keyPointsAndBounds2, bmaIndex2, x1, y1, x2, y2, 
            redImg1, greenImg1, blueImg1, redImg2, greenImg2, blueImg2);
    }

    public FeatureComparisonStat matchHalfDescriptors(
        IntensityClrFeatures features1, IntensityClrFeatures features2, 
        KeyPointsAndBounds keyPointsAndBounds1, int bmaIndex1,
        KeyPointsAndBounds keyPointsAndBounds2, int bmaIndex2,
        int x1, int y1, int x2, int y2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, 
        GreyscaleImage blueImg1, GreyscaleImage redImg2, 
        GreyscaleImage greenImg2, GreyscaleImage blueImg2) {
        
        int rot1, rot2;
        try {
            rot1 = features1.calculateOrientation(x1, y1);
            rot2 = features2.calculateOrientation(x2, y2);
        } catch (CornerRegionDegneracyException e) {
            return null;
        }
                
        IntensityDescriptor desc2_l = features2.extractIntensityLOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_l == null) {
            return null;
        }
        IntensityDescriptor desc2_a = features2.extractIntensityAOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_a == null) {
            return null;
        }
        IntensityDescriptor desc2_b = features2.extractIntensityBOfCIELAB(redImg2, 
            greenImg2, blueImg2, x2, y2, rot2);
        if (desc2_b == null) {
            return null;
        }
            
        IntensityDescriptor desc1_l = features1.extractIntensityLOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_l == null) {
            return null;
        }
        IntensityDescriptor desc1_a = features1.extractIntensityAOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_a == null) {
            return null;
        }
        IntensityDescriptor desc1_b = features1.extractIntensityBOfCIELAB(redImg1, 
            greenImg1, blueImg1, x1, y1, rot1);
        if (desc1_b == null) {
            return null;
        }
        
        /*
                  90
           135    |    45
                  |
        180 ---------------  0
                  |
           225    |    315
                 270                
        */
        // determine which half of the descriptor to use, the top half which
        // is in the direction of the orientation, or the bottom half which
        // is 180 from that direction.
        
        BlobMedialAxes bma1 = keyPointsAndBounds1.getBoundingRegions().getBlobMedialAxes();
        PairInt xySkel1 = bma1.findClosestPoint(bmaIndex1, x1, y1);
        PairInt xyCen1 = bma1.getOriginalBlobXYCentroid(bmaIndex1);

        // direction away from skeleton or centroid
        int thetaOut1;
        if ((x1 != xySkel1.getX()) || (y1 != xySkel1.getY())) {
            double theta = Math.atan2(y1 - xySkel1.getY(), x1 - xySkel1.getX());
            // transform to 0 to 2*pi radians
            if (theta < 0) {
                theta += 2. * Math.PI;
            } 
            thetaOut1 = (int)Math.round(theta * 180./Math.PI);
        } else {
            double theta = Math.atan2(y1 - xyCen1.getY(), x1 - xyCen1.getX());
            // transform to 0 to 2*pi radians
            if (theta < 0) {
                theta += 2. * Math.PI;
            } 
            thetaOut1 = (int)Math.round(theta * 180./Math.PI);
        }

        BlobMedialAxes bma2 = keyPointsAndBounds2.getBoundingRegions().getBlobMedialAxes();
        PairInt xySkel2 = bma2.findClosestPoint(bmaIndex2, x2, y2);
        PairInt xyCen2 = bma2.getOriginalBlobXYCentroid(bmaIndex2);

        // direction away from skeleton or centroid
        int thetaOut2;
        if ((x2 != xySkel2.getX()) || (y2 != xySkel2.getY())) {
            double theta = Math.atan2(y2 - xySkel2.getY(), x2 - xySkel2.getX());
            // transform to 0 to 2*pi radians
            if (theta < 0) {
                theta += 2. * Math.PI;
            } 
            thetaOut2 = (int)Math.round(theta * 180./Math.PI);
        } else {
            double theta = Math.atan2(y2 - xyCen2.getY(), x2 - xyCen2.getX());
            // transform to 0 to 2*pi radians
            if (theta < 0) {
                theta += 2. * Math.PI;
            } 
            thetaOut2 = (int)Math.round(theta * 180./Math.PI);
        }
        
        boolean useTop1 = Math.abs(AngleUtil.getAngleDifference(rot1, thetaOut1)) < 90;
        
        boolean useTop2 = Math.abs(AngleUtil.getAngleDifference(rot2, thetaOut2)) < 90;
                
        FeatureComparisonStat stat_deltaE = IntensityClrFeatures.calculateHalfStats(
            desc1_l, desc1_a, desc1_b, x1, y1, useTop1,
            desc2_l, desc2_a, desc2_b, x2, y2, useTop2);
        
        if (Float.isNaN(stat_deltaE.getSumIntensitySqDiff()) || 
            Float.isNaN(stat_deltaE.getImg2PointIntensityErr())
            ) {
            return null;
        }
        
        return stat_deltaE;
    }

}
