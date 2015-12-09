package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.ClosestPairBetweenSets;
import algorithms.compGeometry.ClosestPairBetweenSets.ClosestPairInt;
import algorithms.compGeometry.FurthestPair;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.LinearRegression;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class HoughTransformTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public HoughTransformTest(String testName) {
        super(testName);
    }

    public void estLines0() throws Exception {

        String fileName1, fileName2;

        for (int i = 5; i < 6; ++i) {
            //fileName1 = "valve_gaussian.png";
            //fileName2 = "valve_gaussian.png";
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    break;
                }
            }

            System.out.println("fileName1=" + fileName1);

            String bin = ResourceFinder.findDirectory("bin");
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            int idx = fileName1.lastIndexOf(".");
            String fileNameRoot = fileName1.substring(0, idx);

            ImageProcessor imageProcessor = new ImageProcessor();

            boolean useBinned = true;
            int binnedImageMaxDimension = 512;
            int binFactor1 = 1;
            int binFactor2 = 1;

            GreyscaleImage img1 = ImageIOHelper.readImage(filePath1).copyToGreyscale();
            GreyscaleImage img2 = ImageIOHelper.readImage(filePath2).copyToGreyscale();

            if (useBinned) {
                binFactor1 = (int) Math.ceil(Math.max((float) img1.getWidth() / binnedImageMaxDimension,
                (float) img1.getHeight() / binnedImageMaxDimension));

                binFactor2 = (int) Math.ceil(Math.max((float) img2.getWidth() / binnedImageMaxDimension,
                (float) img2.getHeight() / binnedImageMaxDimension));

                img1 = imageProcessor.binImage(img1, binFactor1);
                img2 = imageProcessor.binImage(img2, binFactor2);
            }

            //ImageSegmentation imageSegmentation = new ImageSegmentation();
            //GreyscaleImage segImg1 = imageSegmentation.createGreyscale5(img1);
            //GreyscaleImage segImg2 = imageSegmentation.createGreyscale5(img2);

            BlobPerimeterHelper bph = new BlobPerimeterHelper(
                ImageIOHelper.readImageExt(filePath1), fileNameRoot);
            bph.createBinnedGreyscaleImage(binnedImageMaxDimension);
            bph.applySegmentation(SegmentationType.GREYSCALE_HIST, useBinned);
            BlobCornerHelper bch = new BlobCornerHelper(bph, fileNameRoot);

            bch.turnOffCorrectionForLineArtifacts();

            List<List<CornerRegion>> cornerRegionLists =
                bch.generatePerimeterCorners(SegmentationType.GREYSCALE_HIST,
                useBinned);

            GreyscaleImage segImg1 = useBinned ?
                bph.getBinnedSegmentationImage(SegmentationType.GREYSCALE_HIST) :
                bph.getSegmentationImage(SegmentationType.GREYSCALE_HIST);

            ImageSegmentation imageSegmentation = new ImageSegmentation();
            GreyscaleImage segImg2 = imageSegmentation.createGreyscale5(img2);

            String outPath1 = bin + "/seg_1_" + fileNameRoot +".png";
            String outPath2 = bin + "/seg_2_" + fileNameRoot +".png";
            ImageIOHelper.writeOutputImage(outPath1, segImg1);
            ImageIOHelper.writeOutputImage(outPath2, segImg2);

            // a look at finding lines with Hough transform
            ZhangSuenLineThinner lineThinner = new ZhangSuenLineThinner();
            lineThinner.applyFilter(segImg1);
            String outPath1_0 = bin + "/seg_1__lt_" + fileNameRoot + ".png";
            ImageIOHelper.writeOutputImage(outPath1_0, segImg1);

            algorithms.compGeometry.HoughTransform ht0 =
                new algorithms.compGeometry.HoughTransform();
            Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap =
                ht0.calculateLineGivenEdges(segImg1);
            List<PairInt> outSortedKeys = ht0.sortByVotes(outputPolarCoordsPixMap);

            Image tmp1SegImg1 = segImg1.copyToColorGreyscale();

            for (int ii = 0; ii < outSortedKeys.size(); ++ii) {
                int[] c = ImageIOHelper.getNextRGB(ii);
                PairInt thetaRadius = outSortedKeys.get(ii);
                Set<PairInt> pixCoords = outputPolarCoordsPixMap.get(thetaRadius);

                for (PairInt p : pixCoords) {
                    tmp1SegImg1.setRGB(p.getX(), p.getY(), c[0], c[1], c[2]);
                }
            }
            ImageIOHelper.writeOutputImage(
                bin + "/seg_1_hough1_" + fileNameRoot + ".png", tmp1SegImg1);

            Map<PairInt, Set<PairInt>> polarCoordsPixMapOrig =
                new HashMap<PairInt, Set<PairInt>>(outputPolarCoordsPixMap);

            //TODO: the radiusTol probably has to be >= FWZI of sigma used in gradient
            int thetaTol = 2;
            int radiusTol = 20/binFactor1;
            int sizeLimit = 50/(binFactor1*binFactor1);

            //TODO: if there were lines found, remove them and try w/ a larger radius to
            //   capture the horizontal lines

            List<Set<PairInt>> outputSortedGroups = new ArrayList<Set<PairInt>>();
            Map<PairInt, PairInt> pixToTRMap = ht0.createPixTRMapsFromSorted(
                outSortedKeys, outputPolarCoordsPixMap, outputSortedGroups,
                thetaTol, radiusTol);
            Map<PairInt, Integer> trColorMap = new HashMap<PairInt, Integer>();
            Image tmp2SegImg1 = segImg1.copyToColorGreyscale();
            for (int col = 0; col < tmp2SegImg1.getWidth(); ++col) {
                for (int row = 0; row < tmp2SegImg1.getHeight(); ++row) {
                    int v = tmp2SegImg1.getR(col, row);
                    if (v > 0) {
                        PairInt p = new PairInt(col, row);
                        PairInt tr = pixToTRMap.get(p);
                        Integer cIndex = trColorMap.get(tr);
                        if (cIndex == null) {
                            cIndex = Integer.valueOf(trColorMap.size());
                            trColorMap.put(tr, cIndex);
                        }
                        int cIdx = cIndex.intValue();
                        int[] c = ImageIOHelper.getNextRGB(cIdx);
                        tmp2SegImg1.setRGB(col, row, c[0], c[1], c[2]);
                    }
                }
            }
            ImageIOHelper.writeOutputImage(
                bin + "/seg_1_hough2_" + fileNameRoot + ".png", tmp2SegImg1);

            Image tmp3SegImg1 = segImg1.copyToColorGreyscale();
            for (int ii = 0; ii < outputSortedGroups.size(); ++ii) {

                Set<PairInt> group = outputSortedGroups.get(ii);

                if (group.size() < sizeLimit) {
                    break;
                }

                int[] c = ImageIOHelper.getNextRGB(ii);

                for (PairInt p : group) {
                    tmp3SegImg1.setRGB(p.getX(), p.getY(), c[0], c[1], c[2]);
                }
            }

            for (List<CornerRegion> cornerRegions : cornerRegionLists) {
                for (CornerRegion cr : cornerRegions) {
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    ImageIOHelper.addPointToImage(x, y, tmp3SegImg1, 0, 0, 1,
                        255, 255, 255);
                }
            }
            for (List<CornerRegion> cornerRegions : cornerRegionLists) {
                for (CornerRegion cr : cornerRegions) {
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    ImageIOHelper.addPointToImage(x, y, tmp3SegImg1, 0, 0, 0,
                        255, 0, 0);
                }
            }
            ImageIOHelper.writeOutputImage(
                bin + "/seg_1_hough3_" + fileNameRoot + ".png", tmp3SegImg1);

            // make inverse map from polarCoordsPixMapOrig
            Map<PairInt, PairInt> origPixToTRMap = new HashMap<PairInt, PairInt>();
            for (Entry<PairInt, Set<PairInt>> entry : polarCoordsPixMapOrig.entrySet()) {
                for (PairInt p : entry.getValue()) {
                    origPixToTRMap.put(p, entry.getKey());
                }
            }

            Set<PairInt> group0 = outputSortedGroups.get(1);
            int n = group0.size();
            int count = 0;
            float[] t = new float[n];
            float[] r = new float[n];
            for (PairInt p : group0) {
                PairInt tr = origPixToTRMap.get(p);
                t[count] = tr.getX();
                r[count] = tr.getY();
                count++;
            }
            float[] tAvgStDev = MiscMath.getAvgAndStDev(t);
            float[] rAvgStDev = MiscMath.getAvgAndStDev(r);

            float minXT = (float)Math.floor(tAvgStDev[0] - 3*tAvgStDev[1]);
            float maxXT = (float)Math.ceil(tAvgStDev[0] + 3*tAvgStDev[1]);
            if (maxXT == minXT) {
                maxXT += 2;
            }

            float minXR = (float)Math.floor(rAvgStDev[0] - 3*rAvgStDev[1]);
            float maxXR = (float)Math.ceil(rAvgStDev[0] + 3*rAvgStDev[1]);
            if (maxXR == minXR) {
                maxXR += 2;
            }

            HistogramHolder histT = Histogram.createSimpleHistogram(minXT, maxXT,
                1.f, t, Errors.populateYErrorsBySqrt(t));
            HistogramHolder histR = Histogram.createSimpleHistogram(minXR, maxXR,
                1.f, r, Errors.populateYErrorsBySqrt(r));

            List<Integer> indexesT = MiscMath.findStrongPeakIndexes(histT, 0.25f);
            List<Integer> indexesR = MiscMath.findStrongPeakIndexes(histR, 0.25f);

            histT.plotHistogram("theta", "_theta");
            histR.plotHistogram("radius", "_radius");

            int z0 = 1;
        }
    }

    public void testLines1() throws Exception {

        String fileName1, fileName2;

        for (int i = 0; i < 1; ++i) {
            //fileName1 = "valve_gaussian.png";
            //fileName2 = "valve_gaussian.png";
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    break;
                }
            }

            System.out.println("fileName1=" + fileName1);

            String bin = ResourceFinder.findDirectory("bin");
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            int idx = fileName1.lastIndexOf(".");
            String fileNameRoot = fileName1.substring(0, idx);

            ImageProcessor imageProcessor = new ImageProcessor();

            boolean useBinned = true;
            int binnedImageMaxDimension = 512;
            int binFactor1 = 1;
            int binFactor2 = 1;
            SegmentationType type = SegmentationType.GREYSCALE_HIST;

            GreyscaleImage img1 = ImageIOHelper.readImage(filePath1).copyToGreyscale();
            GreyscaleImage img2 = ImageIOHelper.readImage(filePath2).copyToGreyscale();

            if (useBinned) {
                binFactor1 = (int) Math.ceil(Math.max((float) img1.getWidth() / binnedImageMaxDimension,
                (float) img1.getHeight() / binnedImageMaxDimension));

                binFactor2 = (int) Math.ceil(Math.max((float) img2.getWidth() / binnedImageMaxDimension,
                (float) img2.getHeight() / binnedImageMaxDimension));

                img1 = imageProcessor.binImage(img1, binFactor1);
                img2 = imageProcessor.binImage(img2, binFactor2);
            }

            BlobPerimeterHelper bph = new BlobPerimeterHelper(
                ImageIOHelper.readImageExt(filePath1), fileNameRoot);
            bph.createBinnedGreyscaleImage(binnedImageMaxDimension);
            bph.applySegmentation(type, useBinned);
            BlobCornerHelper bch = new BlobCornerHelper(bph, fileNameRoot);
            
            bch.turnOffCorrectionForLineArtifacts(); // <==========

            List<List<CornerRegion>> cornerRegionLists =
                bch.generatePerimeterCorners(type, useBinned);

            GreyscaleImage segImg1 = useBinned ?
                bph.getBinnedSegmentationImage(type) :
                bph.getSegmentationImage(type);
            
            //NOTE: the kind of line artifact seen here appears to occur for
            // near vertical and near horizontal lines, so will limit
            // corrections to those

            //TODO: the radiusTol probably has to be >= FWZI of sigma used in gradient
            int thetaTol = 2;
            int radiusTol = 4;//20/binFactor1;
            int sizeLimit = 15;//30/(binFactor1*binFactor1);
            
   //TODO:  For edges that are closed circles, need to use that information
   // to distinguish orientation of curves.
   // There is an error introduced when have a narrow blob which has two sides
   // that are nearly parallel.  Without knowledge of orientation along the curve,
   // both of the parallel sides might be found to be the same "hough transform found line" if the
   // radius tolerance is high enough, and then the valid corner in between
   // them gets removed here.

            List<PairIntArray> edgeLists = bph.getBlobPerimeters(type, useBinned);

            //CornerCorrector.removeCornersFromLineArtifacts(edgeLists,
            //    cornerRegionLists, segImg1.getWidth(), segImg1.getHeight());

            ImageSegmentation imageSegmentation = new ImageSegmentation();
            GreyscaleImage segImg2 = imageSegmentation.createGreyscale5(img2);

            String outPath1 = bin + "/seg_1_" + fileNameRoot +".png";
            String outPath2 = bin + "/seg_2_" + fileNameRoot +".png";
            ImageIOHelper.writeOutputImage(outPath1, segImg1);
            ImageIOHelper.writeOutputImage(outPath2, segImg2);

            algorithms.compGeometry.HoughTransform ht0 =
                new algorithms.compGeometry.HoughTransform();

            Image tmp1SegImg1 = segImg1.copyToColorGreyscale();
            tmp1SegImg1 = new Image(segImg1.getWidth(), segImg1.getWidth());

            for (int ii = 0; ii < edgeLists.size(); ++ii) {

                //NOTE: in testable method for this, should allow ability to
                // pass in junctions and not delete corners that are in
                // junctions.
                // For these blob perimeters, there are not junctions.
                
                PairIntArray edge = edgeLists.get(ii);

                List<CornerRegion> cornerRegions = cornerRegionLists.get(ii);

                if (cornerRegions.size() < 2) {
                    continue;
                }
                
                String edgeStr = String.format("edge %d (%3d,%3d):(%3d,%3d) %d pts",
                    ii, edge.getX(0), edge.getY(0), 
                    edge.getX(edge.getN() - 1), edge.getY(edge.getN() - 1),
                    edge.getN());
                System.out.println(edgeStr);

                ImageIOHelper.addCurveToImage(edge, tmp1SegImg1, 1, 255, 255, 255);

                Map<PairInt, Integer> pointEdgeIndexMap = Misc.makePointIndexMap(edge);
            
                Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap =
                   ht0.calculateLineGivenEdge(edge, segImg1.getWidth(),
                       segImg1.getHeight());

                List<PairInt> outSortedKeys = ht0.sortByVotes(outputPolarCoordsPixMap);

                // === find indiv lines within the edge ====
                
                //Map<PairInt, Set<PairInt>> polarCoordsPixMapOrig =
                //    new HashMap<PairInt, Set<PairInt>>(outputPolarCoordsPixMap);

                List<Set<PairInt>> outputSortedGroups = new ArrayList<Set<PairInt>>();
                Map<PairInt, PairInt> pixToTRMap = ht0.createPixTRMapsFromSorted(
                    outSortedKeys, outputPolarCoordsPixMap, outputSortedGroups,
                    thetaTol, radiusTol);
               
                for (int iii = 0; iii < outputSortedGroups.size(); ++iii) {

                    Set<PairInt> group = outputSortedGroups.get(iii);

                    if (group.size() < sizeLimit) {
                        break;
                    }
                    
                    // quick look at line's theta.
                    // the artifacts are sometimes present in lines near
                    // vertical or near horizontal, so skipping if not those                    
                    PairInt aPoint = group.iterator().next();
                    PairInt tr = pixToTRMap.get(aPoint);
                    int vhLimit = 3;
                    if ((tr.getX() > vhLimit) && (Math.abs(tr.getX() - 180) > vhLimit)
                        && (Math.abs(tr.getX() - 90) > vhLimit)) {
                        continue;
                    }
                                
                    int minGroupIdx = Integer.MAX_VALUE;
                    int maxGroupIdx = Integer.MIN_VALUE;
                    for (PairInt p : group) {
                        int eIdx = pointEdgeIndexMap.get(p).intValue();
                        if (eIdx < minGroupIdx) {
                            minGroupIdx = eIdx;
                        }
                        if (eIdx > maxGroupIdx) {
                            maxGroupIdx = eIdx;
                        }
                    }
                    boolean wrapAround = (minGroupIdx == 0) &&
                        ((maxGroupIdx - minGroupIdx) > (group.size() + 5));
        //  /*
                    if (wrapAround) {
                        // determine points of furthest pair and then
                        // check every corner to see if it is in between and 
                        // larger than 2 pixels from furthest points
                        FurthestPair fp = new FurthestPair();
                        PairInt[] furthest = fp.find(group);
                        if (furthest == null) {
                            continue;
                        }
                        
                        double distBetweenEndPoints = distance(furthest[0], furthest[1]);
                        
                        for (int j = (cornerRegions.size() - 1); j > -1; --j) {
                            CornerRegion cr = cornerRegions.get(j);
                            double dist = distance(cr, furthest[0]);
                            if (dist < 3) {
                                continue;
                            }
                            dist = distance(cr, furthest[1]);
                            if (dist < 3) {
                                continue;
                            }
                            // if it is in between the 2 points, remove it
                            if (isInBetween(furthest, distBetweenEndPoints, cr)) {
                                
                                // print theta, radius for debugging
                                System.out.println("      *** deleting (theta, radius)=(" 
                                    + tr.getX() + "," + tr.getY() + ")");
                                    
                                cornerRegions.remove(cr);
                            }
                        }
                        
                        continue;
                    }
                    
                    int tx0 = (minGroupIdx > -1) ? edge.getX(minGroupIdx) : -1;
                    int ty0 = (minGroupIdx > -1) ? edge.getY(minGroupIdx) : -1;
                    int tx1 = (maxGroupIdx > -1) ? edge.getX(maxGroupIdx) : -1;
                    int ty1 = (maxGroupIdx > -1) ? edge.getY(maxGroupIdx) : -1;
                    String grpStr = String.format(
                        "  line segment %d:  (%3d,%3d) to (%3d,%3d) eIdxes: %d to %d,  %d pts",
                        iii, tx0, ty0, tx1, ty1, minGroupIdx, maxGroupIdx, group.size());
                    System.out.println(grpStr);

                    //the cornerRegions are ordered counter clock wise
                    
                    // NOTE: if junctions were present, would want to skip
                    // deleting a cornerRegion that was in a junction

                    // delete corners that are more than 2 indexes from
                    // line bounds, starting from last corner
                    for (int j = (cornerRegions.size() - 1); j > -1; --j) {
                        CornerRegion cr = cornerRegions.get(j);
                        int crIdx = cr.getIndexWithinCurve();
                        if ((crIdx > (minGroupIdx + 2)) &&
                            (crIdx < (maxGroupIdx - 2))) {
                            
                            // print theta, radius for debugging
                            System.out.println("      *** deleting (theta, radius)=(" 
                                + tr.getX() + "," + tr.getY() + ")");
                            System.out.println("deleting: (" 
                                + cr.getX()[cr.getKMaxIdx()] + "," 
                                + cr.getY()[cr.getKMaxIdx()] + ")  eIdx=" 
                                + cr.getIndexWithinCurve() + " (theta, radius)=(" 
                                + tr.getX() + "," + tr.getY() + ")");
                            
                            cornerRegions.remove(cr);
                        }
                    }
         // */
                    int[] c = ImageIOHelper.getNextRGB(iii);
                    for (PairInt p : group) {
                        tmp1SegImg1.setRGB(p.getX(), p.getY(), c[0], c[1], c[2]);
                    }
                    
                }
                for (int j = 0; j < cornerRegions.size(); ++j) {
                    CornerRegion cr = cornerRegions.get(j);
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    System.out.println("plotting: (" + x + "," + y + ")");
                    ImageIOHelper.addPointToImage(x, y, tmp1SegImg1, 0, 0, 2,
                        255, 0, 0);
                }
            }

            ImageIOHelper.writeOutputImage(
                bin + "/seg_1_hough1_" + fileNameRoot + ".png", tmp1SegImg1);

            int z0 = 1;
        }
    }

    private double distance(CornerRegion cr, PairInt p) {
        int x1 = cr.getX()[cr.getKMaxIdx()];
        int y1 = cr.getY()[cr.getKMaxIdx()];

        int diffX = x1 - p.getX();
        int diffY = y1 - p.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }
    
    private double distance(PairInt p1, PairInt p2) {

        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }

    // assuming cr was tested as further from endPoints than 2 pixels
    private boolean isInBetween(PairInt[] endPoints, double distBetweenEndPoints,
        CornerRegion cr) {
        
        double d0 = distance(cr, endPoints[0]);
        
        double d1 = distance(cr, endPoints[1]);
        
        return (d0 + d1) < distBetweenEndPoints;
    }

}
