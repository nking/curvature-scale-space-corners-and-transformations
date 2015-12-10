package algorithms.imageProcessing;

import algorithms.compGeometry.FurthestPair;
import algorithms.compGeometry.HoughTransform;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class CornerCorrector {

    /**
     * given parallel lists of edges and their corner regions, use a hough
     * transform for lines to help remove corner regions that are due
     * to line artifacts.  NOTE that cornerRegionLists needs to already be
     * ordered counter clockwise (as is done in BlobsAndCorners,
     * see BlobsAndCorners.orderCornersCCW(curve, cornerRegions) if needed).
     * @param edgeLists
     * @param cornerRegionLists
     * @param thetaTolerance
     * @param radiusTolerance
     * @param imageWidth
     * @param imageHeight
     */
    public static void removeCornersFromLineArtifacts(List<PairIntArray> edgeLists,
        List<List<CornerRegion>> cornerRegionLists,
        int thetaTolerance, int radiusTolerance,
        int imageWidth, int imageHeight) {
        
        if (true) {
            throw new UnsupportedOperationException("not yet implemented");
        }

        //TODO: sizeLimit may need to scale by binFactor, so be passed in as arg,
        // but should probably have a minimum size that could include
        // 3 corners, one of which would be looked at for removal
        int sizeLimit = 15;

        HoughTransform ht = new HoughTransform();

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

            Map<PairInt, Integer> pointEdgeIndexMap = Misc.makePointIndexMap(edge);

            Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap =
               ht.calculateLineGivenEdge(edge, imageWidth, imageHeight);

            List<PairInt> outSortedKeys = ht.sortByVotes(outputPolarCoordsPixMap);

            // === find indiv lines within the edge ====

            //Map<PairInt, Set<PairInt>> polarCoordsPixMapOrig =
            //    new HashMap<PairInt, Set<PairInt>>(outputPolarCoordsPixMap);

            List<Set<PairInt>> outputSortedGroups = new ArrayList<Set<PairInt>>();
            Map<PairInt, PairInt> pixToTRMap = ht.createPixTRMapsFromSorted(
                outSortedKeys, outputPolarCoordsPixMap, outputSortedGroups,
                thetaTolerance, radiusTolerance);
               
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
                if ((tr.getX() > vhLimit) && (Math.abs(tr.getX() - 90) > vhLimit)
                    && (Math.abs(tr.getX() - 180) > vhLimit)
                    && (Math.abs(tr.getX() - 270) > vhLimit)
                    && (Math.abs(tr.getX() - 360) > vhLimit)
                    ) {
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
                            cornerRegions.remove(cr);
                        }
                    }

                    continue;
                }

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
                        cornerRegions.remove(cr);
                    }
                }                    
            }
        }
    }

    private static double distance(CornerRegion cr, PairInt p) {
        int x1 = cr.getX()[cr.getKMaxIdx()];
        int y1 = cr.getY()[cr.getKMaxIdx()];

        int diffX = x1 - p.getX();
        int diffY = y1 - p.getY();

        return Math.sqrt(diffX * diffX + diffY * diffY);
    }

    private static double distance(PairInt p1, PairInt p2) {

        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();

        return Math.sqrt(diffX * diffX + diffY * diffY);
    }

    // assuming cr was tested as further from endPoints than 2 pixels
    private static boolean isInBetween(PairInt[] endPoints, double distBetweenEndPoints,
        CornerRegion cr) {

        double d0 = distance(cr, endPoints[0]);

        double d1 = distance(cr, endPoints[1]);

        return (d0 + d1) < distBetweenEndPoints;
    }
}
