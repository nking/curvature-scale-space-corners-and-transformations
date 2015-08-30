package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

public class FeatureMatcher {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public FeatureMatcher() {
    }

    public CorrespondenceList findSimilarFeatures(
        GreyscaleImage img1, GreyscaleImage gXY1, GreyscaleImage theta1,
        CornerRegion[] cornerRegions1,
        GreyscaleImage img2, GreyscaleImage gXY2, GreyscaleImage theta2,
        CornerRegion[] cornerRegions2) {
        
        if (img1 == null) {
            throw new IllegalArgumentException("img1 cannot be null");
        }
        if (img2 == null) {
            throw new IllegalArgumentException("img2 cannot be null");
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
        if (cornerRegions1 == null) {
            throw new IllegalArgumentException("cornerRegions1 cannot be null");
        }
        if (cornerRegions2 == null) {
            throw new IllegalArgumentException("cornerRegions2 cannot be null");
        }
        
        //TODO: need to solve for scale using the scale space curves that went 
        // into making the corners. see doc/contours.pdf in this project.
        // the creation of scale space maps for inflection points takes longer
        // though, so might need to use what exists less precisely.
        // interestingly, to estimate scale well with scale space maps requires
        // finding closed curves ("contours") and then matching the peaks of the
        // scale space maps created from the inflection
        // points and that solves for scale, and translation and rotation
        // very quickly if there are contours in common with both frames and
        // then one would not need this feature based matching.
        // the combination of both is probably a very strong recognizer.
        
        // until then, making an assumption of '1' here for limited use
        float scale = 1.f;
        
        int blockHalfWidth = 5;
        boolean useNormalizedIntensities = true;
        
        Features features1 = new Features(img1, gXY1, theta1, blockHalfWidth, 
            useNormalizedIntensities);
        
        Features features2 = new Features(img2, gXY2, theta2, blockHalfWidth, 
            useNormalizedIntensities);
        
        Map<PairInt, Map<PairInt, Set<FeatureComparisonStat>>> matchingMap =
            findMatchingFeatures(features1, features2, cornerRegions1, 
                cornerRegions2, blockHalfWidth, useNormalizedIntensities);
        
        // a quick rough look at the most frequent parameters of euclidean
        // transformations.
        
        //TODO: this may need alot of adjustment, espec for histogram bin
        // sizes.
        
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
     * @param features1
     * @param features2
     * @param region1
     * @param region2
     * @param dither
     * @return
     */
    public FeatureComparisonStat ditherAndRotateForBestLocation(
        IntensityFeatures features1, IntensityFeatures features2, 
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither) {
        
        final int x1 = region1.getX();
        final int y1 = region1.getY();
        int rot1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX();
        final int y2 = region2.getY();
        int rot2 = Math.round(region2.getRelativeOrientationInDegrees());
        
        FeatureComparisonStat best = null;
        
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, rot2);
        
        if (desc2 == null) {
            return null;
        }
        
        /*
        NOTE: best distinguisher for improved rotation angle is usually 
        the gradient images because these are on edges.
        */
        
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
        BlobPerimeterRegion region1, BlobPerimeterRegion region2, int dither) {
        
        final int x1 = region1.getX();
        final int y1 = region1.getY();
        int rot1 = Math.round(region1.getRelativeOrientationInDegrees());
        
        final int x2 = region2.getX();
        final int y2 = region2.getY();
        int rot2 = Math.round(region2.getRelativeOrientationInDegrees());
        
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
}
