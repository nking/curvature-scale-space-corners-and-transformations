package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
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

    protected void storeInMap(Map<PairInt, Map<PairInt, 
        Set<FeatureComparisonStat>>> comparisonMap, FeatureComparisonStat stat) {
        
        if (stat == null) {
            return;
        }
        
        if (comparisonMap == null) {
            throw new IllegalArgumentException("comparisonMap cannot be null");
        }
        
        PairInt p1 = stat.getImg1Point();
        PairInt p2 = stat.getImg2Point();
        
        Map<PairInt, Set<FeatureComparisonStat>> p1Map = comparisonMap.get(p1);
        
        if (p1Map == null) {
            p1Map = new HashMap<PairInt, Set<FeatureComparisonStat>>();
            comparisonMap.put(p1, p1Map);
        }
        
        Set<FeatureComparisonStat> p2Map = p1Map.get(p2);
        
        if (p2Map == null) {
            p2Map = new HashSet<FeatureComparisonStat>();
            p1Map.put(p2, p2Map);
        }
        
        p2Map.add(stat);
    }

    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param cornerRegion1
     * @param cornerRegion2
     * @param dither
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    protected FeatureComparisonStat ditherAndRotateForCorner1Location(
        Features features1, Features features2, CornerRegion cornerRegion1, 
        CornerRegion cornerRegion2, int dither) throws 
        CornerRegion.CornerRegionDegneracyException {
                
        int kMaxIdx1 = cornerRegion1.getKMaxIdx();
        final int x1 = cornerRegion1.getX()[kMaxIdx1];
        final int y1 = cornerRegion1.getY()[kMaxIdx1];
        int rot1 = Math.round(
            cornerRegion1.getRelativeOrientationInDegrees());
        
        int kMaxIdx2 = cornerRegion2.getKMaxIdx();
        final int x2 = cornerRegion2.getX()[kMaxIdx2];
        final int y2 = cornerRegion2.getY()[kMaxIdx2];
        int rot2 = Math.round(
            cornerRegion2.getRelativeOrientationInDegrees());
        
        FeatureComparisonStat best = null;
                
        GradientDescriptor gDesc2 = features2.extractGradient(x2, y2, rot2);
        
        if (gDesc2 == null) {
            return null;
        }
        
        // TODO: could decide not to find best rotation here and just discard
        // false matches due to wrong orientation at end of comparisons
        int[] rotations = new int[10];
        int i = 0;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(y1d)) {
                        continue;
                    }
                    
                    GradientDescriptor gDesc1 = features1.extractGradient(x1d, 
                        y1d, rotD1);
                    
                    if (gDesc1 == null) {
                        continue;
                    }
                       
                    FeatureComparisonStat stat = Features.calculateGradientStats(
                        gDesc1, x1d, y1d, gDesc2, x2, y2);
                   
                    if (stat.getSumGradientSqDiff() < stat.getImg2PointGradientErr()) {
                       
                        if (best == null) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff()) {
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param gsImg1
     * @param gsImg2
     * @param features1
     * @param features2
     * @param cornerRegion1
     * @param cornerRegion2
     * @param rotationInRadians
     * @param tolRotationInRadians
     * @param dither
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    protected FeatureComparisonStat ditherAndRotateForCorner1Location(
        GreyscaleImage gsImg1, GreyscaleImage gsImg2, 
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion cornerRegion1, CornerRegion cornerRegion2, 
        float rotationInRadians, float tolRotationInRadians, 
        int dither) throws CornerRegion.CornerRegionDegneracyException {
                
        float expectedRotationInDegrees = rotationInRadians 
            * (float)(180./Math.PI);
        final float rotTol = tolRotationInRadians * (float)(180./Math.PI);
                
        int kMaxIdx1 = cornerRegion1.getKMaxIdx();
        final int x1 = cornerRegion1.getX()[kMaxIdx1];
        final int y1 = cornerRegion1.getY()[kMaxIdx1];
        int rot1 = Math.round(
            cornerRegion1.getRelativeOrientationInDegrees());
        
        int kMaxIdx2 = cornerRegion2.getKMaxIdx();
        final int x2 = cornerRegion2.getX()[kMaxIdx2];
        final int y2 = cornerRegion2.getY()[kMaxIdx2];
        int rot2 = Math.round(
            cornerRegion2.getRelativeOrientationInDegrees());
        
        FeatureComparisonStat best = null;
                
        IntensityDescriptor desc2 = features2.extractIntensity(gsImg2, x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        // because anglediff compares closest angles:
        if (expectedRotationInDegrees > 180.) {
            expectedRotationInDegrees = 360 - expectedRotationInDegrees;
        }
        
        // TODO: could decide not to find best rotation here and just discard
        // false matches due to wrong orientation at end of comparisons
        int[] rotations = new int[11];
        rotations[0] = angleForResultDiff(rot2, Math.round(expectedRotationInDegrees));
        int i = 1;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(gsImg1, x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(gsImg1, y1d)) {
                        continue;
                    }
                    
                    // only try rotations within expected rotation limits
                    float rotDescriptors = AngleUtil.getAngleDifference(rotD1, rot2);
                    if (rotDescriptors < 0) {
                        rotDescriptors += 360;
                    }
                    float rotDiff = AngleUtil.getAngleDifference(
                        expectedRotationInDegrees, rotDescriptors);
                    if (Math.abs(rotDiff) > rotTol) {
                        continue;
                    }
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(
                        gsImg1, x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        return null;
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param cornerRegion1
     * @param cornerRegion2
     * @param rotationInRadians
     * @param tolRotationInRadians
     * @param dither
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    protected FeatureComparisonStat ditherAndRotateForCorner1Location(
        Features features1, Features features2, CornerRegion cornerRegion1, 
        CornerRegion cornerRegion2, 
        float rotationInRadians, float tolRotationInRadians, int dither) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        float expectedRotationInDegrees = rotationInRadians 
            * (float)(180./Math.PI);
        final float rotTol = tolRotationInRadians * (float)(180./Math.PI);
                
        int kMaxIdx1 = cornerRegion1.getKMaxIdx();
        final int x1 = cornerRegion1.getX()[kMaxIdx1];
        final int y1 = cornerRegion1.getY()[kMaxIdx1];
        int rot1 = Math.round(
            cornerRegion1.getRelativeOrientationInDegrees());
        
        int kMaxIdx2 = cornerRegion2.getKMaxIdx();
        final int x2 = cornerRegion2.getX()[kMaxIdx2];
        final int y2 = cornerRegion2.getY()[kMaxIdx2];
        int rot2 = Math.round(
            cornerRegion2.getRelativeOrientationInDegrees());
        
        FeatureComparisonStat best = null;
                
        GradientDescriptor gDesc2 = features2.extractGradient(x2, y2, rot2);
        
        if (gDesc2 == null) {
            return null;
        }
        
        // because anglediff compares closest angles:
        if (expectedRotationInDegrees > 180.) {
            expectedRotationInDegrees = 360 - expectedRotationInDegrees;
        }
        
        // TODO: could decide not to find best rotation here and just discard
        // false matches due to wrong orientation at end of comparisons
        int[] rotations = new int[11];
        rotations[0] = angleForResultDiff(rot2, Math.round(expectedRotationInDegrees));
        int i = 1;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(y1d)) {
                        continue;
                    }
                    
                    // only try rotations within expected rotation limits
                    float rotDescriptors = AngleUtil.getAngleDifference(rotD1, rot2);
                    if (rotDescriptors < 0) {
                        rotDescriptors += 360;
                    }
                    float rotDiff = AngleUtil.getAngleDifference(
                        expectedRotationInDegrees, rotDescriptors);
                    if (Math.abs(rotDiff) > rotTol) {
                        continue;
                    }
                    
                    GradientDescriptor gDesc1 = features1.extractGradient(x1d, 
                        y1d, rotD1);
                    
                    if (gDesc1 == null) {
                        continue;
                    }
                       
                    FeatureComparisonStat stat = Features.calculateGradientStats(
                        gDesc1, x1d, y1d, gDesc2, x2, y2);
                   
                    if (stat.getSumGradientSqDiff() < stat.getImg2PointGradientErr()) {
                       
                        if (best == null) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rotD1);
                            best.setImg2PointRotInDegrees(rot2);
                        } else {
                            if (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff()) {
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither) 
        throws CornerRegion.CornerRegionDegneracyException {
     
        final int x1 = region1.getX()[1];
        final int y1 = region1.getY()[1];
        int rot1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[1];
        final int y2 = region2.getY()[1];
        int rot2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rot1, x2, y2, rot2, dither);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2) 
        throws CornerRegion.CornerRegionDegneracyException {
     
        final int x1 = region1.getX()[1];
        final int y1 = region1.getY()[1];
        int rot1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[1];
        final int y2 = region2.getY()[1];
        int rot2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rot1, x2, y2, rot2, dither, img1, img2);
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotD1, x2, y2, rotD2, dither);
    }
    
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotD1, x2, y2, rotD2, dither, img1, img2);
    }
    
    /**
     * note that a degreeIntervals > 20 is not recommended
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param degreeIntervals
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither, 
        int degreeIntervals) throws CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());

        int n = 360/degreeIntervals;
        int[] rotations = new int[n];
        int i = 0;
        for (int rot1 = 0; rot1 < 360; rot1 += degreeIntervals) {
            rotations[i] = rot1 + rotD1;
            if (rotations[i] > 359) {
                rotations[i] -= 360;
            }
            i++;
        }
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotations, x2, y2, rotD2, dither);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param degreeIntervals a value > 20 is not recommended
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither,
        int degreeIntervals) throws CornerRegion.CornerRegionDegneracyException {
     
        final int x1 = region1.getX()[1];
        final int y1 = region1.getY()[1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[1];
        final int y2 = region2.getY()[1];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        int n = 360/degreeIntervals;
        int[] rotations = new int[n];
        int i = 0;
        for (int rot1 = 0; rot1 < 360; rot1 += degreeIntervals) {
            rotations[i] = rot1 + rotD1;
            if (rotations[i] > 359) {
                rotations[i] -= 360;
            }
            i++;
        }
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotations, x2, y2, rotD2, dither);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @param degreeIntervals a value > 20 is not recommended
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2, int degreeIntervals) 
        throws CornerRegion.CornerRegionDegneracyException {
     
        final int x1 = region1.getX()[1];
        final int y1 = region1.getY()[1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[1];
        final int y2 = region2.getY()[1];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        int n = 360/degreeIntervals;
        int[] rotations = new int[n];
        int i = 0;
        for (int rot1 = 0; rot1 < 360; rot1 += degreeIntervals) {
            rotations[i] = rot1 + rotD1;
            if (rotations[i] > 359) {
                rotations[i] -= 360;
            }
            i++;
        }
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotations, x2, y2, rotD2, dither, img1, img2);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @param degreeIntervals a value > 20 is not recommended
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        GreyscaleImage img1, GreyscaleImage img2, int degreeIntervals) 
        throws CornerRegion.CornerRegionDegneracyException {
     
        final int x1 = region1.getX()[region1.getKMaxIdx()];
        final int y1 = region1.getY()[region1.getKMaxIdx()];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[region2.getKMaxIdx()];
        final int y2 = region2.getY()[region2.getKMaxIdx()];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        int n = 360/degreeIntervals;
        int[] rotations = new int[n];
        int i = 0;
        for (int rot1 = 0; rot1 < 360; rot1 += degreeIntervals) {
            rotations[i] = rot1 + rotD1;
            if (rotations[i] > 359) {
                rotations[i] -= 360;
            }
            i++;
        }
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotations, x2, y2, rotD2, dither, img1, img2);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rot1
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int rot1,
        final int x2, final int y2, final int rot2,        
        int dither) {
        
        
        int[] rotations = new int[10];
        int i = 0;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        return ditherAndRotateForBestLocation(features1, features2, x1, y1, 
            rotations, x2, y2, rot2, dither);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rot1
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int rot1,
        final int x2, final int y2, final int rot2,        
        int dither, GreyscaleImage img1, GreyscaleImage img2) {
        
        
        int[] rotations = new int[10];
        int i = 0;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        return ditherAndRotateForBestLocation(features1, features2, x1, y1, 
            rotations, x2, y2, rot2, dither, img1, img2);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rotations
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int[] rotations,
        final int x2, final int y2, final int rot2,        
        int dither) {
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(y1d)) {
                        continue;
                    }
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        return null;
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rotations
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int[] rotations,
        final int x2, final int y2, final int rot2,        
        int dither, GreyscaleImage img1, GreyscaleImage img2) {
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(img2, x2, y2, 
            rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(img1, x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(img1, y1d)) {
                        continue;
                    }
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(img1, 
                        x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        return null;
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
    
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        final int expectedRotationInDegrees, final int rotationTol) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotD1, x2, y2, rotD2, dither, expectedRotationInDegrees,
            rotationTol);
    }
    
    /**
     * 
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @param expectedRotationInDegrees
     * @param rotationTol
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither,
        final int expectedRotationInDegrees, final int rotationTol,
        GreyscaleImage img1, GreyscaleImage img2) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherAndRotateForBestLocation(features1, features2,
            x1, y1, rotD1, x2, y2, rotD2, dither, expectedRotationInDegrees,
            rotationTol, img1, img2);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rot1
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @param expectedRotationInDegrees
     * @param rotationTol
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int rot1,
        final int x2, final int y2, final int rot2,        
        int dither, int expectedRotationInDegrees, final int rotationTol) {
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        // because anglediff compares closest angles:
        if (expectedRotationInDegrees > 180.) {
            expectedRotationInDegrees = 360 - expectedRotationInDegrees;
        }
        
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
        
        int[] rotations = new int[11];
        rotations[0] = angleForResultDiff(rot2, expectedRotationInDegrees);
        int i = 1;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            
            // only try rotations within expected rotation limits
            float rotDescriptors = AngleUtil.getAngleDifference(rotD1, rot2);
            if (rotDescriptors < 0) {
                rotDescriptors += 360;
            }
            float rotDiff = AngleUtil.getAngleDifference(
                expectedRotationInDegrees, rotDescriptors);
            if (Math.abs(rotDiff) > rotationTol) {
                continue;
            }
            
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(y1d)) {
                        continue;
                    }
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        return null;
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation. 
     * @param features1
     * @param x1
     * @param y1
     * @param features2
     * @param rot1
     * @param x2
     * @param dither
     * @param rot2
     * @param y2
     * @param expectedRotationInDegrees
     * @param rotationTol
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    protected FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int rot1,
        final int x2, final int y2, final int rot2,        
        int dither, int expectedRotationInDegrees, final int rotationTol,
        GreyscaleImage img1, GreyscaleImage img2) {
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(img2, x2, y2, 
            rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        // because anglediff compares closest angles:
        if (expectedRotationInDegrees > 180.) {
            expectedRotationInDegrees = 360 - expectedRotationInDegrees;
        }
        
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
                
        int[] rotations = new int[11];
        rotations[0] = angleForResultDiff(rot2, expectedRotationInDegrees);
        int i = 1;
        for (int rotD1 = (rot1 - 30); rotD1 <= (rot1 + 30); rotD1 += 10) {
            rotations[i] = rotD1;
            i++;
        }
        rotations[i] = rot1 + 90;
        i++;
        rotations[i] = rot1 + 180;
        i++;
        rotations[i] = rot1 + 270;
        i++;
        
        for (int rotD1 : rotations) {
            if (rotD1 > 359) {
                rotD1 = rotD1 - 360;
            } else if (rotD1 < 0) {
                rotD1 += 360;
            }
            
            // only try rotations within expected rotation limits
            float rotDescriptors = AngleUtil.getAngleDifference(rotD1, rot2);
            if (rotDescriptors < 0) {
                rotDescriptors += 360;
            }
            float rotDiff = AngleUtil.getAngleDifference(
                expectedRotationInDegrees, rotDescriptors);
            if (Math.abs(rotDiff) > rotationTol) {
                continue;
            }
            
            for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
                if (!features1.isWithinXBounds(img1, x1d)) {
                    continue;
                }
                for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                    if (!features1.isWithinYBounds(img1, y1d)) {
                        continue;
                    }
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(img1, 
                        x1d, y1d, rotD1);
        
                    if (desc1 == null) {
                        return null;
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
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @return
     */
    public FeatureComparisonStat ditherForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither) 
        throws CornerRegion.CornerRegionDegneracyException {
        
        final int x1 = region1.getX()[1];
        final int y1 = region1.getY()[1];
        int rot1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX()[1];
        final int y2 = region2.getY()[1];
        int rot2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherForBestLocation(features1, features2,
            x1, y1, rot1, x2, y2, rot2, dither);
    }
    
    public FeatureComparisonStat ditherForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        CornerRegion region1, CornerRegion region2, int dither) throws 
        CornerRegion.CornerRegionDegneracyException {
        
        int kMaxIdx1 = region1.getKMaxIdx();
        int x1 = region1.getX()[kMaxIdx1];
        int y1 = region1.getY()[kMaxIdx1];
        int rotD1 = Math.round(region1.getRelativeOrientationInDegrees());

        int kMaxIdx2 = region2.getKMaxIdx();
        int x2 = region2.getX()[kMaxIdx2];
        int y2 = region2.getY()[kMaxIdx2];
        int rotD2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        return ditherForBestLocation(features1, features2,
            x1, y1, rotD1, x2, y2, rotD2, dither);
    }
    
    /**
     * comparison of descriptors to tune the center of cornerRegion1 and the
     * orientation.  This doesn't compare the other descriptors to get overall
     * best.
     * @param features1
     * @param features2
     * @param x1
     * @param y1
     * @param rot1
     * @param x2
     * @param y2
     * @param rot2
     * @param dither
     * @return
     */
    public FeatureComparisonStat ditherForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        final int x1, final int y1, final int rot1,
        final int x2, final int y2, final int rot2,        
        int dither) {
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
                
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(y1d)) {
                    continue;
                }
                    
                IntensityDescriptor desc1 = features1.extractIntensity(x1d, y1d, rot1);

                if (desc1 == null) {
                    return null;
                }

                FeatureComparisonStat stat = IntensityFeatures.calculateStats(
                    desc1, x1d, y1d, desc2, x2, y2);

                if (stat.getSumIntensitySqDiff() < stat.getImg2PointIntensityErr()) {

                    if (best == null) {
                        best = stat;
                        best.setImg1PointRotInDegrees(rot1);
                        best.setImg2PointRotInDegrees(rot2);
                    } else {
                        if (best.getSumIntensitySqDiff() > stat.getSumIntensitySqDiff()) {
                            best = stat;
                            best.setImg1PointRotInDegrees(rot1);
                            best.setImg2PointRotInDegrees(rot2);
                        }
                    }
                }
            }
        }
        
        return best;
    }
    
    protected FeatureComparisonStat[] findBestMatch(Features features1, 
        Features features2, CornerRegion cornerRegion1, 
        CornerRegion[] cornerRegions2, int dither) {
        
        FeatureComparisonStat best = null;
        
        /*
        best3 ranks by intensity ssd only.  sometimes matches areas where
        projection has cause foreground or background changes that appear
        much more strongly in gradient and theta, so intensity alone may
        match, but the other two might not.
        */
        FeatureComparisonStat best3 = null;
        
        for (int idx2 = 0; idx2 < cornerRegions2.length; ++idx2) {

            CornerRegion cornerRegion2 = cornerRegions2[idx2];

            if (cornerRegion2 == null) {
                continue;
            }

            try {

                /*
                 for the given corner region1,
                 dither and make small rotation changes to find the 
                 best centering and rotation to minimize the differences
                 in gradient descriptors between cornerRegion1 and cornerRegion2.
                 */
                FeatureComparisonStat stat = ditherAndRotateForCorner1Location(
                    features1, features2, cornerRegion1, cornerRegion2,
                    dither);
               
                if (stat != null) {

                    stat = populateOtherDescriptors(features1, features2, stat,
                        cornerRegion2);
                    
                    if (stat == null) {
                        continue;
                    }

                    int rotD1 = Math.round(stat.getImg1PointRotInDegrees());

                    int rotD2 = Math.round(
                        cornerRegion2.getRelativeOrientationInDegrees());

//float diffRot = AngleUtil.getAngleDifference(rotD1, rotD2);
//log.info("diffRot=" + diffRot + " stat=" + stat.toString());

                    if (fitIsBetter(best, stat)) {
                        best = stat;
                        best.setImg1PointRotInDegrees(rotD1);
                        best.setImg2PointRotInDegrees(rotD2);
                    }
                    
                    if (fitIsBetter3(best3, stat)) {
                        best3 = stat;
                        best3.setImg1PointRotInDegrees(rotD1);
                        best3.setImg2PointRotInDegrees(rotD2);
                    }
                }
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
                //log.log(Level.SEVERE, null, ex);
            }
        }
        
        int n = 0;
        if (best != null) {
            n++;
        }
        if (best3 != null) {
            n++;
        }
        FeatureComparisonStat[] result = new FeatureComparisonStat[n];
        n = 0;
        if (best != null) {
            result[n] = best;
            n++;
        }
        if (best3 != null) {
            result[n] = best3;
            n++;
        }

        return result;
    }
    
    protected FeatureComparisonStat[] findBestMatch(Features features1, 
        Features features2, CornerRegion cornerRegion1, 
        CornerRegion[] cornerRegions2, 
        float rotationInRadians, float tolRotationInRadians, int dither) {
        
        FeatureComparisonStat best = null;
        
        /*
        best3 ranks by intensity ssd only.  sometimes matches areas where
        projection has cause foreground or background changes that appear
        much more strongly in gradient and theta, so intensity alone may
        match, but the other two might not.
        */
        FeatureComparisonStat best3 = null;
        
        for (int idx2 = 0; idx2 < cornerRegions2.length; ++idx2) {

            CornerRegion cornerRegion2 = cornerRegions2[idx2];

            if (cornerRegion2 == null) {
                continue;
            }

            try {

                /*
                 for the given corner region1,
                 dither and make small rotation changes to find the 
                 best centering and rotation to minimize the differences
                 in gradient descriptors between cornerRegion1 and cornerRegion2.
                */
                FeatureComparisonStat stat = ditherAndRotateForCorner1Location(
                    features1, features2, cornerRegion1, cornerRegion2,
                    rotationInRadians, tolRotationInRadians, dither);

                if (stat != null) {

                    //compare all descriptors to find best
                    int x1 = stat.getImg1Point().getX();
                    int y1 = stat.getImg1Point().getY();
                    int rotD1 = Math.round(stat.getImg1PointRotInDegrees());

                    int kMaxIdx2 = cornerRegion2.getKMaxIdx();
                    int x2 = cornerRegion2.getX()[kMaxIdx2];
                    int y2 = cornerRegion2.getY()[kMaxIdx2];
                    int rotD2 = Math.round(
                        cornerRegion2.getRelativeOrientationInDegrees());

                    IntensityDescriptor iDesc1 = features1.extractIntensity(
                        x1, y1, rotD1);
                    if (iDesc1 == null) {
                        continue;
                    }

                    ThetaDescriptor tDesc1 = features1.extractTheta(
                        x1, y1, rotD1);
                    if (tDesc1 == null) {
                        continue;
                    }

                    GradientDescriptor gDesc1 = features1.extractGradient(
                        x1, y1, rotD1);
                    if (gDesc1 == null) {
                        continue;
                    }

                    IntensityDescriptor iDesc2 = features2.extractIntensity(
                        x2, y2, rotD2);
                    if (iDesc2 == null) {
                        continue;
                    }

                    ThetaDescriptor tDesc2 = features2.extractTheta(
                        x2, y2, rotD2);
                    if (tDesc2 == null) {
                        continue;
                    }

                    GradientDescriptor gDesc2 = features2.extractGradient(
                        x2, y2, rotD2);
                    if (gDesc2 == null) {
                        continue;
                    }

                    stat = Features.calculateStats(iDesc1, gDesc1, tDesc1,
                        x1, y1, iDesc2, gDesc2, tDesc2, x2, y2);

//float diffRot = AngleUtil.getAngleDifference(rotD1, rotD2);
//log.info("diffRot=" + diffRot + " stat=" + stat.toString());

                    if (fitIsBetter(best, stat)) {
                        best = stat;
                        best.setImg1PointRotInDegrees(rotD1);
                        best.setImg2PointRotInDegrees(rotD2);
                    }
                    
                    if (fitIsBetter3(best3, stat)) {
                        best3 = stat;
                        best3.setImg1PointRotInDegrees(rotD1);
                        best3.setImg2PointRotInDegrees(rotD2);
                    }
                }
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
                //log.log(Level.SEVERE, null, ex);
            }
        }
        
        int n = 0;
        if (best != null) {
            n++;
        }
        if (best3 != null) {
            n++;
        }
        FeatureComparisonStat[] result = new FeatureComparisonStat[n];
        n = 0;
        if (best != null) {
            result[n] = best;
            n++;
        }
        if (best3 != null) {
            result[n] = best3;
            n++;
        }

        return result;
    }

    /**
     * compares best to stat and returns true if stat is better.
     * comparison uses all 3 descriptors, that is intensity, gradient, and
     * theta.
     * @param best
     * @param stat
     * @return 
     */
    protected boolean fitIsBetter(FeatureComparisonStat best, 
        FeatureComparisonStat stat) {
        
        if (stat == null) {
            return false;
        }
        
        if ((stat.getSumIntensitySqDiff() <= stat.getImg2PointIntensityErr()) &&
            (stat.getSumGradientSqDiff() <= stat.getImg2PointGradientErr())
            && (stat.getSumThetaSqDiff() <= stat.getImg2PointThetaErr())
            ) {

            if (best == null) {
                return true;
            } else {
                if (
                    (best.getSumIntensitySqDiff() >= stat.getSumIntensitySqDiff())
                    && 
                    (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff())
                    && (best.getSumThetaSqDiff() > stat.getSumThetaSqDiff())
                    ) {
                    return true;
                }
            }
        }
        
        return false;
    }

    /**
     * compares best to stat and returns true if stat is better.
     * comparison checks that all 3 descriptors are smaller than their errors,
     * then uses intensity, gradient to find best.
     * @param best
     * @param stat
     * @return 
     */
    protected boolean fitIsBetter2(FeatureComparisonStat best, 
        FeatureComparisonStat stat) {
        
        if (stat == null) {
            return false;
        }
        
        if ((stat.getSumIntensitySqDiff() <= stat.getImg2PointIntensityErr()) &&
            (stat.getSumGradientSqDiff() <= stat.getImg2PointGradientErr())
            && (stat.getSumThetaSqDiff() <= stat.getImg2PointThetaErr())
            ) {

            if (best == null) {
                return true;
            } else {
                if (
                    (best.getSumIntensitySqDiff() >= stat.getSumIntensitySqDiff())
                    && 
                    (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff())
                    ) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    /**
     * compares best to stat and returns true if stat is better.
     * comparison checks that all 3 descriptors are smaller than their errors,
     * then uses intensity, gradient to find best.
     * @param best
     * @param stat
     * @return 
     */
    protected boolean fitIsBetter3(FeatureComparisonStat best, 
        FeatureComparisonStat stat) {
        
        if (stat == null) {
            return false;
        }
        
        if ((stat.getSumIntensitySqDiff() <= stat.getImg2PointIntensityErr()) &&
            (stat.getSumGradientSqDiff() <= stat.getImg2PointGradientErr())
            && (stat.getSumThetaSqDiff() <= stat.getImg2PointThetaErr())
            ) {

            if (best == null) {
                return true;
            } else {
                if (
                    (best.getSumIntensitySqDiff() >= stat.getSumIntensitySqDiff())
                    ) {
                    return true;
                }
            }
        }
        
        return false;
    }

    protected Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>> 
        findMatchingFeatures(Features features1, Features features2, 
        CornerRegion[] cornerRegions1, CornerRegion[] cornerRegions2, 
        int blockHalfWidth, boolean useNormalizedIntensities) {
        
        Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>>
            comparisonMap = new HashMap<PairInt, Map<PairInt, 
            Set<FeatureComparisonStat>>>();
                    
        final int dither = 1;
                
        for (int idx1 = 0; idx1 < cornerRegions1.length; ++idx1) {
            
            CornerRegion cornerRegion1 = cornerRegions1[idx1];
            
            if (cornerRegion1 == null) {
                continue;
            }
            
            FeatureComparisonStat[] best = findBestMatch(features1, features2,
                cornerRegion1, cornerRegions2, dither);
            
            if ((best != null) && (best.length > 0)) {
                for (FeatureComparisonStat stat : best) {
                    storeInMap(comparisonMap, stat);
                }
            }
        }
        
        return comparisonMap;
    }
        
    protected Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>> 
        findMatchingFeatures(Features features1, Features features2, 
        CornerRegion[] cornerRegions1, CornerRegion[] cornerRegions2, 
        float rotationInRadians, float tolRotationInRadians,
        int blockHalfWidth, boolean useNormalizedIntensities) {
        
        Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>>
            comparisonMap = new HashMap<PairInt, Map<PairInt, 
            Set<FeatureComparisonStat>>>();
                    
        final int dither = 1;
                
        for (int idx1 = 0; idx1 < cornerRegions1.length; ++idx1) {
            
            CornerRegion cornerRegion1 = cornerRegions1[idx1];
            
            if (cornerRegion1 == null) {
                continue;
            }

            FeatureComparisonStat[] best = findBestMatch(features1, features2,
                cornerRegion1, cornerRegions2, rotationInRadians, 
                tolRotationInRadians, dither);
            
            if ((best != null) && (best.length > 0)) {
                for (FeatureComparisonStat stat : best) {
                    storeInMap(comparisonMap, stat);
                }
            }
        }
        
        return comparisonMap;
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

    public static float calculateDiffThetaMean(List<FeatureComparisonStat> 
        comparisonStats) {
        
        float[] values = calculateThetaDiff(comparisonStats);
        
        if (values == null || values.length == 0) {
            return Float.POSITIVE_INFINITY;
        }
        
        double sum = 0;
        
        for (int i = 0; i < values.length; ++i) {
            sum += values[i];
        }
        
        return (float) (sum / ((float) values.length));
    }
    
    public CorrespondenceList findSimilarFeatures(GreyscaleImage gsImg1, 
        CornerRegion[] cr1, GreyscaleImage gsImg2, CornerRegion[] cr2, 
        TransformationParameters params, float scaleTol, 
        float rotationInRadiansTol, int transXYTol) { 
        
        int dither = 1;
        
        return findSimilarFeatures(gsImg1, cr1, gsImg2, cr2, params, scaleTol, 
            rotationInRadiansTol, transXYTol, dither);
    }

    public CorrespondenceList findSimilarFeatures(GreyscaleImage gsImg1, 
        CornerRegion[] cr1, GreyscaleImage gsImg2, 
        CornerRegion[] cr2, TransformationParameters params, float scaleTol, 
        float rotationInRadiansTol, int transXYTol, int dither) {
        
        List<CornerRegion> filteredTransformedC1 = new ArrayList<CornerRegion>();
        List<CornerRegion> filteredC1 = new ArrayList<CornerRegion>();
        List<CornerRegion> filteredC2 = new ArrayList<CornerRegion>();
        
        filterForIntersection3(params, transXYTol, 
            cr1, cr2, filteredTransformedC1, filteredC1, filteredC2);
        
        if (true) {
            try {
                MiscDebug.writeImage(filteredC1, gsImg1.copyToColorGreyscale(),
                    "filtered_1_corners_");
                MiscDebug.writeImage(filteredC2, gsImg2.copyToColorGreyscale(), 
                    "filtered_2_corners_");
            } catch (IOException ex) {
                Logger.getLogger(FeatureMatcherWrapper.class.getName()).log(
                    Level.SEVERE, null, ex);
            }
        }
        
        /*
        when transformation params are known ahead of time:
        cr1 can be transformed into crTr1 (including the internal points).
        then the matching is faster than n1 * n2 because can discard some 
        possible matches immediately.
        bipartite matching when points are present within tolerance.
        
        bipartite is n^3 but the n is < n1.
        */
       
        final int blockHalfWidth = 5;
        final boolean useNormalizedIntensities = true;
        
        IntensityFeatures features1 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities);
        
        IntensityFeatures features2 = new IntensityFeatures(blockHalfWidth,
            useNormalizedIntensities);
        
        //Features features1 = new Features(gsImg1, gXY1, theta1, blockHalfWidth, 
        //    useNormalizedIntensities);
        //
        //Features features2 = new Features(gsImg2, gXY2, theta2, blockHalfWidth, 
        //    useNormalizedIntensities);
        
        int n1 = filteredC1.size();
        int n2 = filteredC2.size();
        int nMaxMatchable = Math.min(n1, n2);
        
        Map<PairInt, FeatureComparisonStat> statMap = null;
        float[][] cost = null;
                
        final boolean useBipartite = (nMaxMatchable < 251);
        
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
            
            double bestCost = Double.MAX_VALUE;
            int bestIdx2 = -1;
            FeatureComparisonStat bestStat = null;
            
            for (int i2 = 0; i2 < n2; ++i2) {
                
                CornerRegion c2 = filteredC2.get(i2);
                
                int x2 = c2.getX()[c2.getKMaxIdx()];
                int y2 = c2.getY()[c2.getKMaxIdx()];
                                
                int diffX = Math.abs(x1Tr - x2);
                int diffY = Math.abs(y1Tr - y2);
                if (diffX > transXYTol || diffY > transXYTol) {
                    continue;
                }
                
                try {
                  
                    // use the untransformed cr1 to be able to filter by rotation
                    FeatureComparisonStat stat = ditherAndRotateForCorner1Location(
                        gsImg1, gsImg2, features1, features2, c1, c2,
                        params.getRotationInRadians(), rotationInRadiansTol, dither);
                    
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
                    
                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                    log.fine(ex.getMessage());
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
            return useBipartiteMatching(cost, statMap, params);
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
        
        //PairIntArray matchedXY1 = new PairIntArray(nc);
        //PairIntArray matchedXY2 = new PairIntArray(nc);
        
        List<PairInt> matched1 = new ArrayList<PairInt>();
        List<PairInt> matched2 = new ArrayList<PairInt>();

        float[] weights = new float[nc];
        double sumW = 0;

        nc = 0;
        for (Entry<Integer, FeatureComparisonStat> entry : index1StatMap.entrySet()) {
            FeatureComparisonStat fcs = entry.getValue();
            int idx1 = fcs.getIndex1();
            int idx2 = fcs.getIndex2();
            PairInt p1 = fcs.getImg1Point();
            PairInt p2 = fcs.getImg2Point();
            matched1.add(p1.copy());
            matched2.add(p2.copy());
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

    public CorrespondenceList findSimilarFeatures(GreyscaleImage gsImg1, 
        GreyscaleImage gXY1, GreyscaleImage theta1, CornerRegion[] cr1, 
        GreyscaleImage gsImg2, GreyscaleImage gXY2, GreyscaleImage theta2, 
        CornerRegion[] cr2, float scale, float rotationInRadians, 
        float scaleTol, float rotationInRadiansTol) {
        
        /*
        for this one, the rotation can be used in the feature matcher
        to discard some pairs
        */
        
        if (gsImg1 == null) {
            throw new IllegalArgumentException("gsImg1 cannot be null");
        }
        if (gsImg2 == null) {
            throw new IllegalArgumentException("gsImg2 cannot be null");
        }
        if (gXY1 == null) {
            throw new IllegalArgumentException("gXY1 cannot be null");
        }
        if (gXY2 == null) {
            throw new IllegalArgumentException("gXY2 cannot be null");
        }
        if (theta1 == null) {
            throw new IllegalArgumentException("theta1 cannot be null");
        }
        if (theta2 == null) {
            throw new IllegalArgumentException("theta2 cannot be null");
        }
        if (cr1 == null) {
            throw new IllegalArgumentException("cr1 cannot be null");
        }
        if (cr2 == null) {
            throw new IllegalArgumentException("cr2 cannot be null");
        }
        
        int blockHalfWidth = 5;
        boolean useNormalizedIntensities = true;
        
        //NOTE: scale was not used in the feature descriptors, but
        // should be in the future.  For large scale factor, would want to
        // scale the grid sizes in the descriptor
        
        Features features1 = new Features(gsImg1, gXY1, theta1, blockHalfWidth, 
            useNormalizedIntensities);
        
        Features features2 = new Features(gsImg2, gXY2, theta2, blockHalfWidth, 
            useNormalizedIntensities);
        
        Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>> matchingMap =
            findMatchingFeatures(features1, features2, cr1, cr2, 
            rotationInRadians, rotationInRadiansTol, blockHalfWidth, 
            useNormalizedIntensities);
        
        // a quick rough look at the most frequent parameters of euclidean
        // transformations.
        
        //TODO: this may need alot of adjustment, espec for histogram bin sizes.
        
        List<PairInt> points1 = new ArrayList<PairInt>();
        List<PairInt> points2 = new ArrayList<PairInt>();
        List<Integer> rotations = new ArrayList<Integer>();
        List<Integer> translationsX = new ArrayList<Integer>();
        List<Integer> translationsY = new ArrayList<Integer>();
        
        for (Entry<PairInt, Map<PairInt, Set<FeatureComparisonStat>>> entry 
            : matchingMap.entrySet()) {
            
            PairInt p1 = entry.getKey();
            
            Map<PairInt, Set<FeatureComparisonStat>> p2Map = entry.getValue();
            
            for (Entry<PairInt, Set<FeatureComparisonStat>> entry2 : p2Map.entrySet()) {
                
                PairInt p2 = entry2.getKey();
                
                //TODO: consider only writing all best solutions if they
                // are different here:
                
                for (FeatureComparisonStat stat : entry2.getValue()) {
                   
                    float rotation = AngleUtil.getAngleDifference(
                        stat.getImg2PointRotInDegrees(), 
                        stat.getImg1PointRotInDegrees());
                    
                    if (rotation < 0) {
                        rotation += 360;
                    }
                    double rotationRadians = rotation * Math.PI/180.;
                    
                    double cosine = Math.cos(rotationRadians);
                    double sine = Math.sin(rotationRadians);
                    
                    //given scale and rotation, can calculate implied rotation
                    double xr = (p1.getX() * scale * cosine)
                        + (p1.getY() * scale * sine);

                    double yr = -(p1.getX() * scale * sine)
                        + (p1.getY() * scale * cosine);

                    int xt = Math.round(p2.getX() - (float)xr);
                    int yt = Math.round(p2.getY() - (float)yr);
                    
                    rotations.add(Integer.valueOf(Math.round(rotation)));
                    translationsX.add(Integer.valueOf(xt));
                    translationsY.add(Integer.valueOf(yt));
                    points1.add(p1);
                    points2.add(p2);
                }
            }
        }
        
        List<PairInt> matched1 = new ArrayList<PairInt>();
        List<PairInt> matched2 = new ArrayList<PairInt>();
        
        HistogramHolder rHist = Histogram.createSimpleHistogram(20, rotations);
        List<Integer> subsetRotation = new ArrayList<Integer>();
        List<Integer> subsetTransX = new ArrayList<Integer>();
        List<Integer> subsetTransY = new ArrayList<Integer>();
        List<PairInt> subsetPoints1 = new ArrayList<PairInt>();
        List<PairInt> subsetPoints2 = new ArrayList<PairInt>();
        int yMaxIdx = MiscMath.findYMaxIndex(rHist.getYHist());
        if (yMaxIdx != -1) {
            float rotMax = rHist.getXHist()[yMaxIdx];
            for (int i = 0; i < rotations.size(); ++i) {
                float rot = rotations.get(i);
                if (Math.abs(rot - rotMax) <= 20) {
                    subsetRotation.add(Math.round(rot));
                    subsetTransX.add(translationsX.get(i));
                    subsetTransY.add(translationsY.get(i));
                    subsetPoints1.add(points1.get(i));
                    subsetPoints2.add(points2.get(i));
                }
            }
            // A general euclidean transformation (possibly varying across
            // image due to projection effects) is needed to determine if
            // points are matched.
            
            // quick look at results to decide if need to use statistics of
            // fit here (as weights) or whether frequency of values is enough.
            
            HistogramHolder txHist = Histogram.createSimpleHistogram(50, subsetTransX);
            try {
                txHist.plotHistogram("translationX", MiscDebug.getCurrentTimeFormatted());
            } catch (IOException ex) {
                Logger.getLogger(FeatureMatcher.class.getName()).log(Level.SEVERE, null, ex);
            }
            int yMaxIdx2 = MiscMath.findYMaxIndex(txHist.getYHist());
            if (yMaxIdx2 != -1) {
                float txMax = txHist.getXHist()[yMaxIdx2];
                float txFWHM = Histogram.measureFWHM(txHist, yMaxIdx2);
                if (txFWHM < 50) {
                    txFWHM = 50;
                }
                List<Integer> subset2Rotation = new ArrayList<Integer>();
                List<Integer> subset2TransX = new ArrayList<Integer>();
                List<Integer> subset2TransY = new ArrayList<Integer>();
                List<PairInt> subset2Points1 = new ArrayList<PairInt>();
                List<PairInt> subset2Points2 = new ArrayList<PairInt>();
                for (int i = 0; i < subsetTransX.size(); ++i) {
                    int tx = subsetTransX.get(i);
                    if (Math.abs(tx - txMax) <= (txFWHM/2.f)) {
                        subset2Rotation.add(subsetRotation.get(i));
                        subset2TransX.add(tx);
                        subset2TransY.add(subsetTransY.get(i));
                        subset2Points1.add(subsetPoints1.get(i));
                        subset2Points2.add(subsetPoints2.get(i));
                    }
                }
                HistogramHolder tyHist = Histogram.createSimpleHistogram(50, subset2TransY);
                try {
                    tyHist.plotHistogram("translationY", MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(FeatureMatcher.class.getName()).log(Level.SEVERE, null, ex);
                }
                int yMaxIdx3 = MiscMath.findYMaxIndex(tyHist.getYHist());
                if (yMaxIdx3 != -1) {
                    float tyMax = tyHist.getXHist()[yMaxIdx3];
                    float tyFWHM = Histogram.measureFWHM(tyHist, yMaxIdx3);
                    if (tyFWHM < 50) {
                        tyFWHM = 50;
                    }
                    
                    log.info("solution rotation=" + rotMax + "+-20 "
                    + " translationX=" + txMax + "+-" + (txFWHM/2.)
                    + " translationY=" + tyMax + "+=" + (tyFWHM/2.));
                    for (int i = 0; i < subset2TransY.size(); ++i) {
                        int ty = subset2TransY.get(i);
                        if (Math.abs(ty - tyMax) <= (tyFWHM/2.f)) {
                            matched1.add(subset2Points1.get(i));
                            matched2.add(subset2Points2.get(i));
                        }
                    } 
                   
                    CorrespondenceList cl = new CorrespondenceList(
                        scale, 
                        Math.round(rotMax), Math.round(txMax), Math.round(tyMax),
                        20, Math.round(txFWHM), Math.round(tyFWHM), 
                        matched1, matched2);
                    
                    return cl;
                }
            }
        }
        
        return null;
    }

    private FeatureComparisonStat populateOtherDescriptors(Features 
        features1, Features features2, FeatureComparisonStat stat, 
        CornerRegion cornerRegion2) throws CornerRegion.CornerRegionDegneracyException {

        if (features1 == null) {
            throw new IllegalArgumentException("features1 cannot be null");
        }
        if (features2 == null) {
            throw new IllegalArgumentException("features2 cannot be null");
        }
        if (stat == null) {
            throw new IllegalArgumentException("stat cannot be null");
        }
        if (cornerRegion2 == null) {
            throw new IllegalArgumentException("cornerRegion2 cannot be null");
        }
        
        //compare all descriptors to find best
        int x1 = stat.getImg1Point().getX();
        int y1 = stat.getImg1Point().getY();
        int rotD1 = Math.round(stat.getImg1PointRotInDegrees());

        int kMaxIdx2 = cornerRegion2.getKMaxIdx();
        int x2 = cornerRegion2.getX()[kMaxIdx2];
        int y2 = cornerRegion2.getY()[kMaxIdx2];
        int rotD2 = Math.round(
            cornerRegion2.getRelativeOrientationInDegrees());

        IntensityDescriptor iDesc1 = features1.extractIntensity(
            x1, y1, rotD1);
        if (iDesc1 == null) {
            return null;
        }

        ThetaDescriptor tDesc1 = features1.extractTheta(
            x1, y1, rotD1);
        if (tDesc1 == null) {
            return null;
        }

        GradientDescriptor gDesc1 = features1.extractGradient(
            x1, y1, rotD1);
        if (gDesc1 == null) {
            return null;
        }

        IntensityDescriptor iDesc2 = features2.extractIntensity(
            x2, y2, rotD2);
        if (iDesc2 == null) {
            return null;
        }

        ThetaDescriptor tDesc2 = features2.extractTheta(
            x2, y2, rotD2);
        if (tDesc2 == null) {
            return null;
        }

        GradientDescriptor gDesc2 = features2.extractGradient(
            x2, y2, rotD2);
        if (gDesc2 == null) {
            return null;
        }

        return Features.calculateStats(iDesc1, gDesc1, tDesc1,
            x1, y1, iDesc2, gDesc2, tDesc2, x2, y2);
    }
    
    public static <T extends CornerRegion> void filterForIntersection(
        TransformationParameters params, float toleranceXY,
        List<List<T>> c1, List<List<T>> c2, 
        List<List<T>> outFilteredTransformedC1,
        List<List<T>> outFilteredC1, 
        List<List<T>> outFilteredC2) {
        
        float[] minXY2 = MiscMath.findMinXY(c2);
        float[] maxXY2 = MiscMath.findMaxXY(c2);
        
        Transformer transformer = new Transformer();
        
        float[] minXY1 = new float[]{Float.MAX_VALUE, Float.MAX_VALUE};
        float[] maxXY1 = new float[]{Float.MIN_VALUE, Float.MIN_VALUE};
   
        for (int i = 0; i < c1.size(); ++i) {
            List<T> outList = new ArrayList<T>();
            List<T> listTr = new ArrayList<T>();
            List<T> list = c1.get(i);
            for (int j = 0; j < list.size(); ++j) {
                
                T ctr = transformer.applyTransformation(params, 
                    list.get(j));
                
                int x = ctr.getX()[ctr.getKMaxIdx()];
                int y = ctr.getY()[ctr.getKMaxIdx()];
                
                if ((x < (minXY2[0] - toleranceXY)) || (x > (maxXY2[0] + toleranceXY))) {
                    continue;
                }
                if ((y < (minXY2[1] - toleranceXY)) || (y > (maxXY2[1] + toleranceXY))) {
                    continue;
                }
                listTr.add(ctr);
                outList.add(list.get(j));
                
                if (x < minXY1[0]) {
                    minXY1[0] = x;
                }
                if (y < minXY1[1]) {
                    minXY1[1] = y;
                }
                if (x > maxXY1[0]) {
                    maxXY1[0] = x;
                }
                if (y > maxXY1[1]) {
                    maxXY1[1] = y;
                }
            }
            outFilteredTransformedC1.add(listTr);
            outFilteredC1.add(outList);
        }

        for (int i = 0; i < c2.size(); ++i) {
            List<T> outList = new ArrayList<T>();
            List<T> list = c2.get(i);
            for (int j = 0; j < list.size(); ++j) {
                T c = list.get(j);
                int x = c.getX()[c.getKMaxIdx()];
                int y = c.getY()[c.getKMaxIdx()];
                
                if ((x < (minXY1[0] - toleranceXY)) || (x > (maxXY1[0] + toleranceXY))) {
                    continue;
                }
                if ((y < (minXY1[1] - toleranceXY)) || (y > (maxXY1[1] + toleranceXY))) {
                    continue;
                }
                outList.add(list.get(j));
            }
            outFilteredC2.add(outList);
        }        
    }
    
    public static void filterForIntersection2(TransformationParameters params, 
        int transXYTol, 
        List<List<CurvatureScaleSpaceContour>> c1, 
        List<List<CurvatureScaleSpaceContour>> c2, 
        List<List<CurvatureScaleSpaceContour>> outFilteredTransformedC1, 
        List<List<CurvatureScaleSpaceContour>> outFilteredC1, 
        List<List<CurvatureScaleSpaceContour>> outFilteredC2) {
        
        float[] minXY2 = MiscMath.findMinXY2(c2);
        float[] maxXY2 = MiscMath.findMaxXY2(c2);
        
        Transformer transformer = new Transformer();
        
        float[] minXY1 = new float[]{Float.MAX_VALUE, Float.MAX_VALUE};
        float[] maxXY1 = new float[]{Float.MIN_VALUE, Float.MIN_VALUE};
        
        for (int i = 0; i < c1.size(); ++i) {
            List<CurvatureScaleSpaceContour> outList = new ArrayList<CurvatureScaleSpaceContour>();
            List<CurvatureScaleSpaceContour> listTr = new ArrayList<CurvatureScaleSpaceContour>();
            List<CurvatureScaleSpaceContour> list = c1.get(i);
            for (int j = 0; j < list.size(); ++j) {
                
                CurvatureScaleSpaceContour ctr = transformer.applyTransformation(params, 
                    list.get(j));
                
                int x = ctr.getPeakDetails()[0].getXCoord();
                int y = ctr.getPeakDetails()[0].getYCoord();
                
                if ((x < (minXY2[0] - transXYTol)) || (x > (maxXY2[0] + transXYTol))) {
                    continue;
                }
                if ((y < (minXY2[1] - transXYTol)) || (y > (maxXY2[1] + transXYTol))) {
                    continue;
                }
                listTr.add(ctr);
                outList.add(list.get(j));
                
                if (x < minXY1[0]) {
                    minXY1[0] = x;
                }
                if (y < minXY1[1]) {
                    minXY1[1] = y;
                }
                if (x > maxXY1[0]) {
                    maxXY1[0] = x;
                }
                if (y > maxXY1[1]) {
                    maxXY1[1] = y;
                }
            }
            outFilteredTransformedC1.add(listTr);
            outFilteredC1.add(outList);
        }

        for (int i = 0; i < c2.size(); ++i) {
            List<CurvatureScaleSpaceContour> outList = new ArrayList<CurvatureScaleSpaceContour>();
            List<CurvatureScaleSpaceContour> list = c2.get(i);
            for (int j = 0; j < list.size(); ++j) {
                CurvatureScaleSpaceContour c = list.get(j);
                int x = c.getPeakDetails()[0].getXCoord();
                int y = c.getPeakDetails()[0].getYCoord();
                
                if ((x < (minXY1[0] - transXYTol)) || (x > (maxXY1[0] + transXYTol))) {
                    continue;
                }
                if ((y < (minXY1[1] - transXYTol)) || (y > (maxXY1[1] + transXYTol))) {
                    continue;
                }
                outList.add(list.get(j));
            }
            outFilteredC2.add(outList);
        }        
    }

    public static void filterForIntersection3(TransformationParameters params, 
        int transXYTol, CornerRegion[] c1, CornerRegion[] c2, 
        List<CornerRegion> outFilteredTransformedC1, 
        List<CornerRegion> outFilteredC1, 
        List<CornerRegion> outFilteredC2) {
        
        int[] minXY2 = MiscMath.findMinXY2(c2);
        int[] maxXY2 = MiscMath.findMaxXY2(c2);
        
        Transformer transformer = new Transformer();
        
        float[] minXY1 = new float[]{Float.MAX_VALUE, Float.MAX_VALUE};
        float[] maxXY1 = new float[]{Float.MIN_VALUE, Float.MIN_VALUE};
        
        for (int i = 0; i < c1.length; ++i) {
                
            CornerRegion ctr = transformer.applyTransformation(params, c1[i]);
                
            int x = ctr.getX()[ctr.getKMaxIdx()];
            int y = ctr.getY()[ctr.getKMaxIdx()];

            if ((x < (minXY2[0] - transXYTol)) || (x > (maxXY2[0] + transXYTol))) {
                continue;
            }
            if ((y < (minXY2[1] - transXYTol)) || (y > (maxXY2[1] + transXYTol))) {
                continue;
            }
            outFilteredTransformedC1.add(ctr);
            outFilteredC1.add(c1[i]);
                
            if (x < minXY1[0]) {
                minXY1[0] = x;
            }
            if (y < minXY1[1]) {
                minXY1[1] = y;
            }
            if (x > maxXY1[0]) {
                maxXY1[0] = x;
            }
            if (y > maxXY1[1]) {
                maxXY1[1] = y;
            }
        }

        for (int i = 0; i < c2.length; ++i) {
            CornerRegion c = c2[i];
            int x = c.getX()[c.getKMaxIdx()];
            int y = c.getY()[c.getKMaxIdx()];
                
            if ((x < (minXY1[0] - transXYTol)) || (x > (maxXY1[0] + transXYTol))) {
                continue;
            }
            if ((y < (minXY1[1] - transXYTol)) || (y > (maxXY1[1] + transXYTol))) {
                continue;
            }
            outFilteredC2.add(c);
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

    public static void removeIntensityOutliers(List<FeatureComparisonStat> compStats) {
        
        if (compStats.size() < 3) {
            return;
        }
        
        int n = compStats.size();
        
        float[] meanStDv = calcIntensitySSDMeanAndStDev(compStats);
        
        List<Integer> rm = new ArrayList<Integer>();
        
        for (int i = 0; i < n; ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            float diff = Math.abs(stat.getSumIntensitySqDiff() - meanStDv[0]);
            
            if (diff > (1.25 * meanStDv[1])) {
                rm.add(Integer.valueOf(i));
            }
        }
        
        for (int i = rm.size() - 1; i > -1; --i) {
            
            int idx = rm.get(i).intValue();
            
            compStats.remove(idx);
        }
    }

    private int angleForResultDiff(int rot2, int expectedRotationInDegrees) {
                
        int r0 = rot2 - expectedRotationInDegrees;
        if (r0 < 0) {
            r0 += 360;
        }
        
        int r = Math.abs(Math.round(AngleUtil.getAngleDifference(r0, rot2)));
        
        if (r != expectedRotationInDegrees) {
            r0 = expectedRotationInDegrees - rot2;
            if (r0 < 0) {
                r0 += 360;
            }
            r = Math.abs(Math.round(AngleUtil.getAngleDifference(r0, rot2)));
        }
        
        return r0;
    }

    private CorrespondenceList useBipartiteMatching(float[][] cost,
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

        List<PairInt> matched1 = new ArrayList<PairInt>();
        List<PairInt> matched2 = new ArrayList<PairInt>();
        
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
            
            int x1 = stat.getImg1Point().getX();
            int y1 = stat.getImg1Point().getY();
            matched1.add(new PairInt(x1, y1));
            int x2 = stat.getImg2Point().getX();
            int y2 = stat.getImg2Point().getY();
            matched2.add(new PairInt(x2, y2));

            nc++;
        }

        if (matched1.size() < 7) {
            return null;
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

}
