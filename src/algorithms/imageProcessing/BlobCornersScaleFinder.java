package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import algorithms.util.ScatterPointPlotterPNG;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 *
 * @author nichole
 */
public class BlobCornersScaleFinder extends AbstractBlobScaleFinder {

    private double ssdLimit = 1500;

    public MatchingSolution solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {

        List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
            type1, useBinned1);
        List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
            type2, useBinned2);
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(
            type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(
            type2, useBinned2);

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);

        assert(blobs1.size() == perimeters1.size());
        assert(blobs1.size() == corners1List.size());
        assert(blobs2.size() == perimeters2.size());
        assert(blobs2.size() == corners2List.size());

        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
  
        float dist = 2.5f;
        for (int i = 0; i < perimeters1.size(); ++i) {
            filterCorners(perimeters1.get(i), corners1List.get(i), dist);
        }
        for (int i = 0; i < perimeters2.size(); ++i) {
            filterCorners(perimeters2.get(i), corners2List.get(i), dist);
        }

        if (true) {
            for (List<CornerRegion> list : corners1List) {
                sortCornersToCCW(list);
            }
            for (List<CornerRegion> list : corners2List) {
                sortCornersToCCW(list);
            }
        }

        MatchingSolution soln = match(features1, features2, img1, img2, perimeters1,
            perimeters2, corners1List, corners2List, binFactor1, binFactor2);

        return soln;
    }

    private <T extends CornerRegion> MatchingSolution match(
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2,
        List<PairIntArray> perimeters1, List<PairIntArray> perimeters2,
        List<List<T>> corners1List, List<List<T>> corners2List,
        int binFactor1, int binFactor2) {

MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
float[] xPoints1 = new float[perimeters1.size()];
float[] yPoints1 = new float[perimeters1.size()];
double[][] xy1 = new double[perimeters1.size()][2];
for (int i = 0; i < perimeters1.size(); ++i) {
xy1[i] = curveHelper.calculateXYCentroids(perimeters1.get(i));
xPoints1[i] = (float)xy1[i][0];
yPoints1[i] = (float)xy1[i][1];
}
float[] xPoints2 = new float[perimeters2.size()];
float[] yPoints2 = new float[perimeters2.size()];
double[][] xy2 = new double[perimeters2.size()][2];
for (int i = 0; i < perimeters2.size(); ++i) {
xy2[i] = curveHelper.calculateXYCentroids(perimeters2.get(i));
xPoints2[i] = (float)xy2[i][0];
yPoints2[i] = (float)xy2[i][1];
}

ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
plotter.plotLabeledPoints(0, img1.getWidth(), 0, img1.getHeight(), xPoints1, yPoints1,
"img1", "X", "Y");
ScatterPointPlotterPNG plotter2 = new ScatterPointPlotterPNG();
plotter2.plotLabeledPoints(0, img2.getWidth(), 0, img2.getHeight(), xPoints2, yPoints2,
"img2", "X", "Y");
try {
    plotter.writeToFile("img1_labelled.png");
    plotter2.writeToFile("img2_labelled.png");
} catch (IOException ex) {
    Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
}

StringBuilder sb = new StringBuilder("xy1:\n");
for (int i = 0; i < xy1.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy1[i][0]), (int)Math.round(xy1[i][1])));
}
sb.append("xy2:\n");
for (int i = 0; i < xy2.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy2[i][0]), (int)Math.round(xy2[i][1])));
}
System.out.println(sb.toString());

/*
PairInt[] im1Chk = new PairInt[]{
    new PairInt(59, 178), new PairInt(42, 110), new PairInt(27, 105),
    new PairInt(68,  80), new PairInt(25,  55), new PairInt(80, 144)
};
PairInt[] im2Chk = new PairInt[]{
    new PairInt(189, 179), new PairInt(164, 109), new PairInt(164, 109),
    new PairInt(189, 179), new PairInt(154,  59), new PairInt(193, 146)
};
int[] im1ChkIdxs = new int[im1Chk.length];
int[] im2ChkIdxs = new int[im2Chk.length];
for (int i = 0; i < im1Chk.length; ++i) {
    PairInt p = im1Chk[i];
    for (int j = 0; j < xy1.length; ++j) {
        double diffX = p.getX() - xy1[j][0];
        double diffY = p.getY() - xy1[j][1];
        if (Math.abs(diffX) < 20 && Math.abs(diffY) < 20) {
            im1ChkIdxs[i] = j;
            break;
        }
    }
}
for (int i = 0; i < im2Chk.length; ++i) {
    PairInt p = im2Chk[i];
    for (int j = 0; j < xy2.length; ++j) {
        double diffX = p.getX() - xy2[j][0];
        double diffY = p.getY() - xy2[j][1];
        if (Math.abs(diffX) < 20 && Math.abs(diffY) < 20) {
            im2ChkIdxs[i] = j;
            break;
        }
    }
}
sb = new StringBuilder("expeected matches:\n");
for (int i = 0; i < im1ChkIdxs.length; ++i) {
    sb.append(String.format("[%d] to [%d]", im1ChkIdxs[i], im2ChkIdxs[i]));
    sb.append("\n");
}
System.out.println(sb.toString());
*/
        Map<PairInt, TransformationParameters> trMap
            = new HashMap<PairInt, TransformationParameters>();

        int n1 = corners1List.size();
        int n2 = corners2List.size();

        if (n1 == 0 || n2 == 0) {
            return null;
        }

        /*
        -- get best TransformationParameters for each idx1
           (this is at most 40 of them)
        -- consider combining similar parameters to reduce the eval step
        -- evaluate each param against all points.
            -- make 2 eval methods, hopefully the fastest is enough
               -- (1) eval by finding an existing point within tolerance of predicted position
               -- (2) eval by finding best SSD of points within tolerance of predicted and
                      use nEval, SSD and dist to return a normalized cost
        */

        for (int idx1 = 0; idx1 < n1; ++idx1) {

            if (corners1List.get(idx1).size() < 3) {
                continue;
            }

            List<T> corners1 = corners1List.get(idx1);

            /*
            first, see if nEval alone finds the true matches for curve to curve
            */
            int maxNEval = Integer.MIN_VALUE;
            TransformationParameters maxNEvalParams = null;
            Integer maxNEvalIndex2 = null;
            double minCost = Double.MAX_VALUE;

            for (int idx2 = 0; idx2 < n2; ++idx2) {

                if (corners2List.get(idx2).size() < 3) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                List<T> corners2 = corners2List.get(idx2);

                ClosedCurveCornerMatcher2<T> mapper =
                    new ClosedCurveCornerMatcher2<T>();

                boolean matched = mapper.matchCorners(features1, features2,
                    corners1, corners2, img1, img2, binFactor1, binFactor2);

                if (!matched) {
                    continue;
                }

                // NOTE: if solving for binned, the params are in binned reference frame
                TransformationParameters params = mapper.getSolution();
                int nEval = mapper.getNEval();
                double cost = mapper.getSolutionCost();

                if (cost < minCost) {
                    if ((nEval < 3) && (maxNEval > 4)) {
                        // do not accept if nEval is much lower than maxNEval
                        continue;
                    }
                    maxNEval = nEval;
                    maxNEvalParams = params;
                    maxNEvalIndex2 = index2;
                    minCost = cost;
                } else if ((maxNEval == 2) && (nEval > 3)) {
                    //TODO: may need to revise this
                    double avgCost = (cost + minCost)/2.;
                    if ((Math.abs(cost - avgCost)/(0.1*avgCost)) < 2) {
                        maxNEval = nEval;
                        maxNEvalParams = params;
                        maxNEvalIndex2 = index2;
                        minCost = cost;
                    }
                }
            }

            if (maxNEvalParams != null) {
                trMap.put(new PairInt(idx1, maxNEvalIndex2.intValue()),
                    maxNEvalParams);
            }
        }

        int n2c = 0;
        for (List<T> corners2 : corners2List) {
            n2c += corners2.size();
        }

        List<T> cr2List = new ArrayList<T>();
        int[] xC2 = new int[n2c];
        int[] yC2 = new int[n2c];
        n2c = 0;
        for (List<T> corners2 : corners2List) {
            for (T cr2: corners2) {
                xC2[n2c] = cr2.getX()[cr2.getKMaxIdx()];
                yC2[n2c] = cr2.getY()[cr2.getKMaxIdx()];
                cr2List.add(cr2);
                n2c++;
            }
        }

        NearestPoints np2 = new NearestPoints(xC2, yC2);

        MatchingSolution soln = evaluateForBestUsingFeatures(
            features1, features2, img1, img2,
            trMap, corners1List, cr2List, np2, binFactor1, binFactor2);

        return soln;
    }

    /**
     * given parameters map, evaluate the corner lists and return the solution.
     * Note that if binFactors are not equal to '1', it's assumed that the
     * paramsMap are in the binned reference frame and corrections to that
     * are made for the returned solution.
     * @param <T>
     * @param features1
     * @param features2
     * @param img1
     * @param img2
     * @param paramsMap
     * @param corners1List
     * @param corners2List
     * @param np2
     * @param binFactor1
     * @param binFactor2
     * @return
     */
    private <T extends CornerRegion> MatchingSolution evaluateForBestUsingFeatures(
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2,
        Map<PairInt, TransformationParameters> paramsMap,
        List<List<T>> corners1List, List<T> corners2List,
        NearestPoints np2, int binFactor1, int binFactor2) {

        int tolTransXY = 5;//10;

        double maxDistance = Math.sqrt(2) * tolTransXY;

        /*
        -- evaluate each param against all points.
           -- (2) eval by finding best SSD of points within tolerance of predicted and
              use nEval, SSD and dist to return a normalized cost
        */

        //List<TransformationParameters> parameterList =
        //    MiscStats.filterToSimilarParamSets(paramsMap, binFactor1, binFactor2);

        List<TransformationParameters> parameterList =
            MiscStats.filterToSimilarParamSets2(paramsMap, binFactor1, binFactor2);

        int n1 = 0;
        for (List<T> corners1 : corners1List) {
            n1 += corners1.size();
        }
        int n2 = corners2List.size();
        int nMaxMatchable = Math.min(n1, n2);

        FeatureMatcher featureMatcher = new FeatureMatcher();
        Transformer transformer = new Transformer();

        final int rotationTolerance = 20;

        final int dither = 4;

        TransformationParameters bestParams = null;
        float bestCost = Float.MAX_VALUE;
        float bestCost1Norm = Float.MAX_VALUE;
        List<FeatureComparisonStat> bestStats = null;

        boolean tolIsTooLarge = false;

        int nIter = 0;

        int deltaTol = 0;

        while (nIter < 2) {

            tolIsTooLarge = false;

            for (TransformationParameters params : parameterList) {

                double rotInRadians = params.getRotationInRadians();
                double cos = Math.cos(rotInRadians);
                double sin = Math.sin(rotInRadians);

                int rotD = Math.round(params.getRotationInDegrees());

                int tolTransXY2 = tolTransXY;
                if (params.getScale() < 1) {
                    tolTransXY2 = Math.round(tolTransXY * params.getScale());
                }
                if (tolTransXY2 == 0) {
                    tolTransXY2 = 1;
                }
                if (tolTransXY2 > 1) {
                    tolTransXY2 -= deltaTol;
                }
                int dither2 = Math.round(dither * params.getScale());
                if (dither2 == 0) {
                    dither2 = 1;
                } else if (dither2 > dither) {
                    // large dither makes runtime larger
                    dither2 = dither;
                }

                int nEval = 0;
                double sumSSD = 0;
                double sumDist = 0;
                List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();

                for (List<T> corners1 : corners1List) {

                    for (T cr : corners1) {

                        T crTr = transformer.applyTransformation(params, cr, cos, sin);

                        Set<Integer> indexes2 = np2.findNeighborIndexes(
                            crTr.getX()[crTr.getKMaxIdx()],
                            crTr.getY()[crTr.getKMaxIdx()], tolTransXY2);

                        for (Integer index : indexes2) {

                            int idx2 = index.intValue();

                            T corner2 = corners2List.get(idx2);

                            FeatureComparisonStat compStat = null;

                            try {
                                compStat = featureMatcher.ditherAndRotateForBestLocation(
                                    features1, features2, cr, corner2, dither2,
                                    rotD, rotationTolerance, img1, img2);
                            } catch (CornerRegion.CornerRegionDegneracyException ex) {
                                log.fine(ex.getMessage());
                            }
                            if (compStat == null || (compStat.getSumIntensitySqDiff() >= ssdLimit)
                                ||
                                (compStat.getSumIntensitySqDiff() >= compStat.getImg2PointIntensityErr())) {
                                continue;
                            }

                            double xTr = (compStat.getImg1Point().getX() *
                                params.getScale() * cos) +
                                (compStat.getImg1Point().getY() *
                                params.getScale() * sin);
                            xTr += params.getTranslationX();

                            double yTr = (-compStat.getImg1Point().getX() *
                                params.getScale() * sin) +
                                (compStat.getImg1Point().getY() *
                                params.getScale()* cos);
                            yTr += params.getTranslationY();

                            double dist = distance(xTr, yTr,
                                compStat.getImg2Point().getX(),
                                compStat.getImg2Point().getY());

                            stats.add(compStat);

                            sumDist += Math.abs(dist);
                            sumSSD += compStat.getSumIntensitySqDiff();
                            nEval++;
                        }
                    }
                }

                if (nEval == 0) {
                    continue;
                }

                // distance needs to be adjusted by scale, else the cost prefers
                // small scale solutions
                sumDist /= params.getScale();

                sumSSD /= (double)nEval;
                sumDist /= (double)nEval;

                float cost1Norm = 1.f/(float)nEval;
                float cost2Norm = (float)(sumSSD/ssdLimit);
                float cost3Norm = (float)sumDist;
                float normalizedCost = cost1Norm * cost2Norm * cost3Norm;
               
                //TODO: cost1Norm's proportion in normalizedCost should be higher
                
                boolean t1 = (normalizedCost < bestCost);
                
                if (t1 && (nEval > 2)) {

                    TransformationParameters combinedParams =
                        MiscStats.calculateTransformation(binFactor1, binFactor2,
                            stats, new float[4]);

                    if (combinedParams == null) {
                        continue;
                    }

                    if (MiscStats.standardDeviationsAreSmall(combinedParams)) {
                        tolIsTooLarge = false;
                        bestCost = normalizedCost;
                        bestParams = combinedParams;
                        bestStats = stats;
                        bestCost1Norm = cost1Norm;
                    } else if (MiscStats.standardDeviationsAreSmall(params)) {
                        //TODO: this suggests tolTransXY2 is too large
                        tolIsTooLarge = true;
                        bestCost = normalizedCost;
                        bestParams = params;
                        bestStats = stats;
                        bestCost1Norm = cost1Norm;
                    }
                }
            }
            if (!tolIsTooLarge) {
                break;
            }
            deltaTol++;

            nIter++;
        }

        if (bestParams != null) {

            if (binFactor1 != 1 || binFactor2 != 1) {
                for (int i = 0; i < bestStats.size(); ++i) {
                    FeatureComparisonStat stat = bestStats.get(i);
                    stat.setBinFactor1(binFactor1);
                    stat.setBinFactor2(binFactor2);
                }
            }

            MatchingSolution soln = new MatchingSolution(bestParams, bestStats);
            return soln;
        }

        return null;
    }

    protected void filterCorners(PairIntArray curve,
        List<CornerRegion> regions, float dist) {

        /*
        if there are more than 1 corner within dist of 2 or so of on another,
        remove all except strongest corner.
        */
        List<Set<Integer>> closeCornerIndexes = findCloseCorners(dist, regions);

        if (closeCornerIndexes.isEmpty()) {
            return;
        }

        List<Integer> remove = new ArrayList<Integer>();
        for (Set<Integer> set : closeCornerIndexes) {
            float maxK = Float.MIN_VALUE;
            Integer maxKIndex = null;
            for (Integer index : set) {
                CornerRegion cr = regions.get(index.intValue());
                float k = cr.getK()[cr.getKMaxIdx()];
                if (k > maxK) {
                    maxK = k;
                    maxKIndex = index;
                }
            }
            for (Integer index : set) {
                if (!index.equals(maxKIndex)) {
                    remove.add(index);
                }
            }
        }

        if (remove.size() > 1) {
            Collections.sort(remove);
        }
        for (int i = (remove.size() - 1); i > -1; --i) {
            regions.remove(remove.get(i).intValue());
        }
    }

    private static List<Set<Integer>> findCloseCorners(float tolD,
        List<CornerRegion> regions) {

        float[] x = new float[regions.size()];
        float[] y = new float[regions.size()];
        for (int i = 0; i < regions.size(); ++i) {
            CornerRegion cr = regions.get(i);
            x[i] = cr.getX()[cr.getKMaxIdx()];
            y[i] = cr.getY()[cr.getKMaxIdx()];
        }

        List<Set<Integer>> close = new ArrayList<Set<Integer>>();

        FixedDistanceGroupFinder groupFinder = new FixedDistanceGroupFinder(x, y);

        groupFinder.findGroupsOfPoints(tolD);

        int nGroups = groupFinder.getNumberOfGroups();

        for (int i = 0; i < nGroups; ++i) {
            Set<Integer> group = groupFinder.getGroupIndexes(i);
            if (group.size() > 1) {
                close.add(group);
            }
        }

        return close;
    }

    private double distance(double x1, double y1, double x2, double y2) {

        double diffX = x1 - x2;
        double diffY = y1 - y2;

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }

    private void sortCornersToCCW(List<CornerRegion> corners) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        PairIntArray cornerXY = new PairIntArray();
        for (int ii = 0; ii < corners.size(); ++ii) {
            CornerRegion cr = corners.get(ii);
            cornerXY.add(cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
        }
        boolean isCW = curveHelper.curveIsOrderedClockwise(cornerXY);
        if (isCW) {
            int n = corners.size();
            if (n > 1) {
                int end = n >> 1;
                for (int ii = 0; ii < end; ii++) {
                    int idx2 = n - ii - 1;
                    CornerRegion swap = corners.get(ii);
                    corners.set(ii, corners.get(idx2));
                    corners.set(idx2, swap);
                }
            }
        }
    }
}
