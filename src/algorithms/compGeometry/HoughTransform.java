package algorithms.compGeometry;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.ClosestPairBetweenSets.ClosestPairInt;
import algorithms.imageProcessing.ConnectedGroupsWithGapsFinder;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.DFSSimilarThetaRadiusGroupsFinder;
import algorithms.imageProcessing.DFSConnectedGroupsFinder2;
import algorithms.imageProcessing.DFSConnectedHoughTransformGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.util.MiscMath;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;

/**
 * a class for Hough transforms of simple geometric shapes.  currently a line
 * is implemented.
 * 
 * Helpful in starting this was a look at the code available from
 * http://vase.essex.ac.uk/software/HoughTransform/
 * and 
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/flatjavasrc/Hough.java
 * though this implementation is different.
 * 
 * @author nichole
 */
public class HoughTransform {

    public HoughTransform() {
    }

    /**
     * given an edge of points, computes the Hough
     * transform of lines and returns results as an associate array with 
     * key = pair with x = polar theta in degrees and y = distance from
     * the origin in pixels; value = number of transformation points having
     * the key.   Note that the angle is calculated for expectations of a
     * counter clockwise ordered curve and the vector of the angle is 
     * perpendicular to p1 (direction given by right hand rule).
     * The angles are 0 to 360.
     * 
     * Note that if the edge has less than 3 points, an empty map is returned.
     * 
     * runtime complexity is O(N_edge_pts), but includes transcendental operations.
     * 
     * @param edge a curve defined by the points within
     * @param imageWidth
     * @param imageHeight
     * @return thetaRadiusPixCoords mappings
     */
    public Map<PairInt, Set<PairInt>> calculateLineGivenEdge(PairIntArray edge,
        int imageWidth, int imageHeight) {
        
        Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = new HashMap<PairInt, Set<PairInt>>();
        
        if (edge.getN() < 3) {
            return outputPolarCoordsPixMap;
        }
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
        
        boolean curveIsClosed = (edge instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)edge).isClosedCurve());

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                    
        int n = edge.getN();
        
        for (int i = 0; i < n; ++i) {
            
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            int xp, yp, xn, yn;
            
            if (i == 0) {
                if (curveIsClosed) {
                    xp = edge.getX(n - 1);
                    yp = edge.getY(n - 1);
                } else {
                    // use replication for boundary
                    xp = x;
                    yp = y;
                }
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            } else if (i == (n - 1)) {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                if (curveIsClosed) {
                    xn = edge.getX(0);
                    yn = edge.getY(0);
                } else {
                    xn = x;
                    yn = y;
                }
            } else {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            }
            
            // note, this is not the angle along the edge, it's perpendicular
            // to it, but the calculation is consistent
            double t = curveHelper.calculateAngleTangentToMidpoint(xp, yp, x, y, 
                xn, yn); 
            
            double tDegrees = t * 180./Math.PI;
            
            int tInt = (int)Math.round(tDegrees);
            
            if (tInt > 359) {
                tInt = tInt - 360;
            }
            
            Integer theta = Integer.valueOf(tInt);

            double ct = cosineMap.get(theta).doubleValue();
            double st = sineMap.get(theta).doubleValue();

            double r = (x * ct) + (y * st);

            if (r < 0) {
                r *= -1;
            }

            PairInt p = new PairInt(tInt, (int)Math.round(r));

            Set<PairInt> set = outputPolarCoordsPixMap.get(p);
            if (set == null) {
                set = new HashSet<PairInt>();
                outputPolarCoordsPixMap.put(p, set);
            }
            set.add(new PairInt(x, y));
        }
        
        return outputPolarCoordsPixMap;
    }
    
    /**
     * given lists of corner regions, computes the Hough
     * transform of lines and returns results lists of polar theta in degrees
     * and radius as distance from
     * the origin in pixels.
     * Note that the angle is calculated for expectations of a
     * counter clockwise ordered curve and the vector of the angle is 
     * perpendicular to p1 (direction given by right hand rule).
     * The angles are 0 to 360.
     * 
     * Note that if the edge has less than 3 points, an empty map is returned.
     * 
     * @param cornerLists
     * @param edges
     * @return thetaRadiusPixCoords lists
     */
    public List<List<PairInt>> calculateRoughHoughTransforms(
        List<List<CornerRegion>> cornerLists, List<PairIntArray> edges) {
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
                        
        List<List<PairInt>> trLists = new ArrayList<List<PairInt>>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                
        for (int cIdx = 0; cIdx < cornerLists.size(); ++cIdx) {
        
            PairIntArray edge = edges.get(cIdx);
            int nEdge = edge.getN();
            
            boolean curveIsClosed = (edge instanceof PairIntArrayWithColor) &&
                (((PairIntArrayWithColor)edge).isClosedCurve());

            List<CornerRegion> cornerList = cornerLists.get(cIdx);
            
            List<PairInt> trList = new ArrayList<PairInt>();
            
            for (int i = 0; i < cornerList.size(); ++i) {
            
                CornerRegion cr = cornerList.get(i);
                
                int eIdx = cr.getIndexWithinCurve();
                if (eIdx == -1) {
                    continue;
                }
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                
                int xp, yp, xn, yn;

                if (eIdx == 0) {
                    if (curveIsClosed) {
                        xp = edge.getX(nEdge - 1);
                        yp = edge.getY(nEdge - 1);
                    } else {
                        // use replication for boundary
                        xp = x;
                        yp = y;
                    }
                    xn = edge.getX(eIdx + 1);
                    yn = edge.getY(eIdx + 1);
                } else if (eIdx == (nEdge - 1)) {
                    xp = edge.getX(eIdx - 1);
                    yp = edge.getY(eIdx - 1);
                    if (curveIsClosed) {
                        xn = edge.getX(0);
                        yn = edge.getY(0);
                    } else {
                        xn = x;
                        yn = y;
                    }
                } else {
                    xp = edge.getX(eIdx - 1);
                    yp = edge.getY(eIdx - 1);

                    xn = edge.getX(eIdx + 1);
                    yn = edge.getY(eIdx + 1);
                }

                // note, this is not the angle along the edge, it's perpendicular
                // to it, but the calculation is consistent
                double t = curveHelper.calculateAngleTangentToMidpoint(xp, yp, x, y,
                    xn, yn);

                double tDegrees = t * 180. / Math.PI;

                int tInt = (int) Math.round(tDegrees);

                if (tInt > 359) {
                    tInt = tInt - 360;
                }

                Integer theta = Integer.valueOf(tInt);

                double ct = cosineMap.get(theta).doubleValue();
                double st = sineMap.get(theta).doubleValue();

                double r = (x * ct) + (y * st);

                if (r < 0) {
                    r *= -1;
                }

                PairInt p = new PairInt(tInt, (int) Math.round(r));

                trList.add(p);
            }
            trLists.add(trList);
        }
        
        return trLists;
    }
    
    /**
     * given the theta values and a set of point coordinates, returns a map
     * of the angle and distance from the image origin.
     * 
     * @param points
     * @param theta360
     * @return 
     */
    public Map<PairInt, PairInt> calculatehHoughTransforms(
        Set<PairInt> points, GreyscaleImage theta360) {
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
                        
        Map<PairInt, PairInt> pointTRMap = new HashMap<PairInt, PairInt>();
                
        for (PairInt p : points) {
        
            int x = p.getX();
            int y = p.getY();
            int t = theta360.getValue(x, y);
            
            /*if ((t - 90) > 0) {
                t -= 90;
            } else if ((t + 90) < 360) {
                t += 90;
            }*/
           
            Integer theta = Integer.valueOf(t);
            
            double ct = cosineMap.get(theta).doubleValue();
            double st = sineMap.get(theta).doubleValue();

            double r = (x * ct) + (y * st);

            if (r < 0) {
                r *= -1;
            }

            PairInt pTR = new PairInt(t, (int) Math.round(r));

            pointTRMap.put(p, pTR);
        }
        
        return pointTRMap;
    }
    
    /**
     * given a set of points, search for contiguous lines within the points
     * and return a map of the groups of points with value being their
     * theta and radius.
     * 
     * not ready for use yet
     * 
     * @param points
     * @return map with key = group of connected points having same theta and
     * radius, value = pairint of the theta and radius for the group.
     */
    public Map<Set<PairInt>, PairInt> findContiguousLines(Set<PairInt> points,
        int minimumGroupSize) {
        
        // key=(x,y); value = derived thetas.  
        // 0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5 where the 0.5 are rounded down
        Map<PairInt, Set<Integer>> pointThetasMap = new HashMap<PairInt, Set<Integer>>();
        
        //key = derived thetas; value = set of points with that theta
        Map<Integer, Set<PairInt>> thetaPointMap = new HashMap<Integer, Set<PairInt>>();
        
        /*
        for each point:
             look at the 24 point neighborhood

               2  2  2  2  2
               2  1  1  1  2
               2  1  @  1  2    test whether adjacent points fit patterns for:
               2  1  1  1  2       0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5
               2  2  2  2  2
               *note, no radius calc is needed yet
               store in 2 maps:  map<pairint, set<integer>> with key=x,y coordinates, value=set of derived theta
                                 map<integer, set<pairint>> with key=theta value=set of points in that group that fit that.
        
        after visiting all points to build the 2 maps, analyze the maps.
        
        for each entry in the theta, set of points map
        
            make a stack out of the set of points and place into contiguous groups.

            visit each of the contiguous groups within the current theta, radius:
                store results in 2 ways:
                     2 parallel lists: one holding the contiguous group of points, 
                         the other holding the theta
                     a map with key being the point x,y and value being a set of 
                         indexes to the parallel list that the point belongs to.

         when finished visiting all of the theta, points map:
             sort the 2 parallel lists into descending number of points,
                 (need to make an array to find the location of a list
                 item by the original unsorted list index)
             put the first list into a "decided" list of groups.
                 - for each point in its group:
                   - use the map to get the indexes to the groups it belongs in
                   - remove the point from each group in the first list that the point is part of
                   - remove the point's entry from the map
              presumably, it is not necessary to resort the parallel lists.
              if did need to re-sort those lists, would have to be careful
              to update the lookup array from original list index to new indexes.
        
              continue through the list of groups until the size is smaller than
                  minimumGroupSize.
              check whether any more of the remaining groups are larger than
                  minimumGroupSize 
                  and if so:
                  re-sort the lists and make sure the original list index 
                      lookup array is updated.
              
           now have groups of contiguous points that have rough theta.
        
           the refined value for the theta could be done in many ways.
               will use a bisector pattern within the range of 22.5 degrees
               and evaluate the angle for the points at each bisection.
        */
        
        boolean[] isPresent = new boolean[24];
        
        for (PairInt p : points) {
            
            populateNeighborhoodPresence(p.getX(), p.getY(), points, isPresent);

      /*
          19 20 21 22 23
          14 15 16 17 18
          10 11  @ 12 13    
           5  6  7  8  9       
           0  1  2  3  4
     */ 
            Set<Integer> thetas = new HashSet<Integer>();
            if ((isPresent[12] && isPresent[13]) || (isPresent[10] && isPresent[11])) {
                Integer t = Integer.valueOf(0);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            }
            if (isPresent[18] || isPresent[5]) {
                Integer t = Integer.valueOf(22);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if ((isPresent[17] && isPresent[23]) || (isPresent[6] && isPresent[0])) {
                Integer t = Integer.valueOf(45);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if (isPresent[22] || isPresent[1]) {
                Integer t = Integer.valueOf(67);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if ((isPresent[21] && isPresent[16]) || (isPresent[7] && isPresent[2])) {
                Integer t = Integer.valueOf(90);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if (isPresent[20] || isPresent[3]) {
                Integer t = Integer.valueOf(112);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if ((isPresent[19] && isPresent[15]) || (isPresent[8] && isPresent[4])) {
                Integer t = Integer.valueOf(135);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            } 
            if (isPresent[14] || isPresent[9]) {
                Integer t = Integer.valueOf(157);
                thetas.add(t);
                Set<PairInt> tPoints = thetaPointMap.get(t);
                if (tPoints == null) {
                    tPoints = new HashSet<PairInt>();
                    thetaPointMap.put(t, tPoints);
                }
                tPoints.add(p);
            }
            pointThetasMap.put(p, thetas);            
        }
        
        // connected groups of points with same thetas
        List<Set<PairInt>> contigGroups = new ArrayList<Set<PairInt>>();
        
        // the thetas of the groups in contigGroups
        List<Integer> contigGroupThetas = new ArrayList<Integer>();
        
        // indexes to groups in contigGroups to which the key point belongs
        Map<PairInt, Set<Integer>> pointIndexesToContigGroups = new HashMap<PairInt, Set<Integer>>();       
        for (Entry<Integer, Set<PairInt>> entry : thetaPointMap.entrySet()) {
                        
            // find contiguous groups
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.setToUse8Neighbors();
            finder.setMinimumNumberInCluster(minimumGroupSize);
            finder.findConnectedPointGroups(entry.getValue());
            
            for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
                
                Set<PairInt> group = finder.getXY(i);
                
                Integer index = Integer.valueOf(contigGroupThetas.size());
                
                contigGroups.add(group);
                contigGroupThetas.add(index);
                
                for (PairInt p : group) {

                    Set<Integer> indexes = pointIndexesToContigGroups.get(p);
                    if (indexes == null) {
                        indexes = new HashSet<Integer>();
                        pointIndexesToContigGroups.put(p, indexes);
                    }
                    indexes.add(index);
                }
            }
        }
        
        /*
        TODO: 
            need to search the aggregated points for chains between them
            that render a different composite angle on a larger scale.
        
            for example, this line is composed of 45 degree segments, but is
               slightly larger than 45 degrees
                              @
                           @
                        @
                        @
                     @
                  @
        
            Or this composed of 0 degree segments
        
                  @  @  @  @  @
                                 @  @  @  @  @
        
        A quick look shows that some lines such as present in obtuse triangles
            have a gap of 1 pixel and some a gap of 2 pixels in their lines 
            composed of staircase segments.
        
        Also, might need a check that parallel connected lines do not get 
        aggregated into one line.
           This can be done in an evaluation stage after thiel sen parameters,
        */
        
        /*
        //NOTE: had to use a generous theta diff of 44 degrees in the merging
        //   here, so may exclude this section in the end
        ConnectedGroupsWithGapsFinder cgFinder = new ConnectedGroupsWithGapsFinder();
        cgFinder.setMinimumNumberInCluster(2);
        cgFinder.findConnectedGroups(contigGroups, contigGroupThetas, 
            Math.sqrt(2.)*2);
        for (int i = 0; i < cgFinder.getNumberOfGroups(); ++i) {
            Set<Integer> groupedIndexes = cgFinder.getGroupedIndexes(i);
            Set<PairInt> combined = new HashSet<PairInt>();
            for (Integer index : groupedIndexes) {
                int idx = index.intValue();
                combined.addAll(contigGroups.get(idx));
            }
            Integer theta = contigGroupThetas.get(groupedIndexes.iterator().next().intValue());
            Integer index = Integer.valueOf(contigGroups.size());
            contigGroups.add(combined);
            contigGroupThetas.add(theta);
            for (PairInt p : combined) {
                Set<Integer> indexes = pointIndexesToContigGroups.get(p);
                if (indexes == null) {
                    indexes = new HashSet<Integer>();
                    pointIndexesToContigGroups.put(p, indexes);
                }
                indexes.add(index);
            }
        }
        */
        
        int[] indexes = new int[contigGroups.size()];
        int[] sizes = new int[contigGroups.size()];
        for (int i = 0; i < indexes.length; ++i) {
            indexes[i] = i;
            sizes[i] = contigGroups.get(i).size();
        }
        MultiArrayMergeSort.sortByDecr(sizes, indexes);
        
        // ---- calculate a refined theta and radius and evaluate, then
        //      store the acceptable solutions and remove them from those
        //      further down the list
        
        Map<Set<PairInt>, PairInt> outputPointsTR = new HashMap<Set<PairInt>, PairInt>();
        
        float tol0 = 0.5f;
        float tol1 = 0.45f;
        
        for (int i = 0; i < indexes.length; ++i) {
            int idx = indexes[i];
            Set<PairInt> group = contigGroups.get(idx);
            if (group.size() >= minimumGroupSize) {
                
                Integer roughTheta = contigGroupThetas.get(idx);
                
                //float[] trMeanStDev = calculateLinePolarCoordsAndStats(group, 
                //    roughTheta);
                
                float[] trMeanStDev = calculateReducedLinePolarCoordsAndStats(group, 
                    roughTheta);
                
                float t = trMeanStDev[0];
                // if it's highly inclined line, tolerance has to allow for step width
                /*if ((Math.abs(AngleUtil.getAngleDifference(t, 22.5f)) < 12) || 
                    (Math.abs(AngleUtil.getAngleDifference(t, 67.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 112.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 157.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 202.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 247.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 292.5f)) < 12) ||
                    (Math.abs(AngleUtil.getAngleDifference(t, 337.5f)) < 12)
                    ) {
                    tol1 = 1.5f;
                }*/
                
                // consider attempting to split parallel connected lines.
                if (trMeanStDev[2] >= 1) {
                    // 3 parallel lines: t=353, mn=1.05, stdv=0.4
                    // corner of 2 intersecting lines: t=311, mn=1.8, st=0.5
                }
                
                if (trMeanStDev[2] > tol0 && trMeanStDev[3] > tol1) {
                    continue;
                }
                //TODO: consider another limit here
                if (group.size() < minimumGroupSize) {
                    continue;
                }
                
                System.out.println("t=" + trMeanStDev[0] + " r=" + trMeanStDev[1] + 
                    " mn=" + trMeanStDev[2] + " stdv=" + trMeanStDev[3]);
                
                PairInt tr = new PairInt((int)trMeanStDev[0], 
                    Math.round(trMeanStDev[1]));
                      
                outputPointsTR.put(group, tr);
                
                for (PairInt p : group) {
                    Set<Integer> index2Set = pointIndexesToContigGroups.get(p);
                    for (Integer index2 : index2Set) {
                        int idx2 = index2.intValue();
                        if (idx2 != idx) {
                            contigGroups.get(idx2).remove(p);
                        }
                    }
                    pointIndexesToContigGroups.remove(p);
                }
            }
        }        
        
        return outputPointsTR;
    }
    
    /**
     * given the theta values and a set of point coordinates, 
     * finds contiguous groups of points with the same theta and distance
     * from origin and returns them as lists of contiguous sets of such
     * points.
     * 
     * should be able to find where the lines connect too as unassigned 
     * adjacent points to the line sets that have values in between
     * two lines... that will be another method
     * 
     * unfortunately, this method is not returning the best of results because
     * the angles in the theta map are determined from a very local 
     * neighborhood when the gradient is made and so even a thinned line has
     * a line with a large range in theta (as large as 20 degrees in a 
     * straight line) and hence the radius calculations have a large range too.
     * Instead, the user should prefer to use
     * Map<Set<PairInt>, PairInt> findContiguousLines(Set<PairInt> points,
        int minimumGroupSize) 
     * which has a larger factor times semi-linear runtime complexity, but the 
     * results should be better.
     * 
     * @param points
     * @param theta360
     * 
     * @return map with key = group of connected points having same theta and
     * radius, value = pairint of the theta and radius for the group.
     */
    public Map<Set<PairInt>, PairInt> findContiguousLines(Set<PairInt> points, 
        GreyscaleImage theta360) {
        
        Map<PairInt, PairInt> pointTRMap = calculatehHoughTransforms(points, 
            theta360);
        
        /*
        String[] fileNames = new String[]{"mask_01.png", "mask_02.png", 
            "mask_03.png", "mask_04.png", "mask_05.png", "mask_06.png", 
            "mask_07.png", "mask_08.png", "mask_09.png", "mask_10.png",
            "mask_11.png", "mask_12.png"};
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            for (String fileName : fileNames) {
                 String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
                 GreyscaleImage mask = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
                 Set<PairInt> edgePoints = new HashSet<PairInt>();
                 for (int i = 0; i < mask.getNPixels(); ++i) {
                     if (mask.getValue(i) > 127) {
                         PairInt p2 = new PairInt(mask.getCol(i), mask.getRow(i));
                         if (pointTRMap.containsKey(p2)) {
                             edgePoints.add(p2);
                         }
                     }
                 }
                 float[] theta = new float[edgePoints.size()];
                 float[] radius = new float[theta.length];
                 int count = 0;
                 for (PairInt p : edgePoints) {
                     PairInt tr = pointTRMap.get(p);
                     //theta[count] = tr.getX();
                     //theta[count] = theta360.getValue(p);
                     //radius[count] = tr.getY();
                     theta[count] = p.getX();
                     radius[count] = p.getY();
                     count++;
                 }
                 float[] xPolygon = null; 
                 float[] yPolygon = null;
                 float minRadius = MiscMath.findMin(radius);
                 float maxRadius = MiscMath.findMax(radius);
                 float minX = MiscMath.findMin(theta);
                 float maxX = MiscMath.findMax(theta);
                 plotter.addPlot(minX - 1, maxX + 1, minRadius - 1, maxRadius + 1, theta,
                     radius, xPolygon, yPolygon, "theta vs radius");
            }
            String filePath = plotter.writeFile();
            int z = 1;
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
        }
        */  
        
        // edit theta in pointTRMap to only use angles between 0 and 180 to
        // remove ambiguous edits made to theta360 or artifacts from the making
        // of the gradient
        for (Entry<PairInt, PairInt> entry : pointTRMap.entrySet()) {
            PairInt tr = entry.getValue();
            int t = tr.getY();
            if (t > 179) {
                tr.setY(t - 180);
            }
        }
        
        DFSConnectedHoughTransformGroupsFinder finder = 
            new DFSConnectedHoughTransformGroupsFinder();
        
        finder.setMinimumNumberInCluster(1);
        
        int thetaTolerance = 10;
        int radiusTolerance = 2;
        int wrapAroundValue = 180;
        
        finder.findConnectedPointGroups(pointTRMap, thetaTolerance,
            radiusTolerance, wrapAroundValue);
        
        Map<Set<PairInt>, PairInt> output = new HashMap<Set<PairInt>, PairInt>();
        
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            
            Set<PairInt> group = finder.getXY(i);
            
            if (group.size() < 2) {
                continue;
            }
            
            int averageTheta = 0;
            int averageRadius = 0;
            for (PairInt p : group) {
                PairInt tr = pointTRMap.get(p);
                averageTheta += tr.getX();
                averageRadius += tr.getY();
            }
            averageTheta /= group.size();
            averageRadius /= group.size();
            
            output.put(group, new PairInt(averageTheta, averageRadius));
        }
        
        return output;
    }
    
    /**
     * runtime complexity is O(N * lg_2(N)).
     * @param thetaRadiusPixMap
     * @return 
     */
    public List<PairInt> sortByVotes(Map<PairInt, Set<PairInt>> thetaRadiusPixMap) {
        
        int[] votes = new int[thetaRadiusPixMap.size()];
        int[] indexes = new int[votes.length];
        PairInt[] keys = new PairInt[votes.length];
        
        int count = 0;
        for (Entry<PairInt, Set<PairInt>> entry : thetaRadiusPixMap.entrySet()) {
            votes[count] = entry.getValue().size();
            keys[count] = entry.getKey();
            indexes[count] = count;
            count++;
        }
        
        MultiArrayMergeSort.sortByDecr(votes, indexes);

        List<PairInt> outSortedKeys = new ArrayList<PairInt>();
        
        for (int i = 0; i < indexes.length; ++i) {
            int idx = indexes[i];            
            outSortedKeys.add(keys[idx]);
        }
        
        return outSortedKeys;
    }
    
    public HoughTransformLines createPixTRMapsFromSorted(List<PairInt> sortedTRKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, 
        List<Set<PairInt>> outputSortedGroups) {
        
        int thetaTol = 2;
        int radiusTol = 8;
        
        return createPixTRMapsFromSorted(sortedTRKeys, thetaRadiusPixMap, 
            thetaTol, radiusTol);
    }

    /**
     * sets true and false in present array for presence in points set.
       <pre>
          19 20 21 22 23
          14 15 16 17 18
          10 11  @ 12 13    
           5  6  7  8  9       
           0  1  2  3  4
       </pre>
     * @param x
     * @param y
     * @param points
     * @param present 
    */
    private void populateNeighborhoodPresence(int x, int y, Set<PairInt> points, 
        boolean[] present) {
        
        int count = 0;
        for (int dy = -2; dy <= 2; ++dy) {
            int y2 = y + dy;
            for (int dx = -2; dx <= 2; ++dx) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                int x2 = x + dx;
                PairInt p = new PairInt(x2, y2);
                if (points.contains(p)) {
                    present[count] = true;
                } else {
                    present[count] = false;
                }
                ++count;
            }
        }
    }

    /**
     *  uses thiel sen estimator to calculate the slope and y intercept,
     * then calculates the average distance of the points from the line
     * and then the standard deviation of that average.  then removes
     * outliers that are further from the line than meanTolerance
     * and then recalculates.  Note that the argument set points is modified
     * to remove outliers.
     * @param linePoints input and output modified set
     * @param roughTheta
     * @param meanTolerance
     * @param stDvTolerance
     * @return 
     */
    private float[] calculateReducedLinePolarCoordsAndStats(
        Set<PairInt> linePoints, Integer roughTheta) {
        
        PairFloatArray xy = new PairFloatArray();
        
        for (PairInt p : linePoints) {
            xy.add(p.getX(), p.getY());
        }
        
        float[] x = Arrays.copyOf(xy.getX(), xy.getN());
        float[] y = Arrays.copyOf(xy.getY(), xy.getN());
        
        LinearRegression lReg = new LinearRegression();
        //lReg.plotTheLinearRegression(x, y);
        
        float[] yInterceptAndSlope = 
            lReg.calculateTheilSenEstimatorParams(x, y);
        
        double thetaRadians;
        double rSum = 0;
        
        boolean isVertical = (yInterceptAndSlope[1] == Float.MAX_VALUE);
        boolean isHorizontal = (yInterceptAndSlope[1] == 0); 
        
        if (isVertical) {
            thetaRadians = Math.PI/2.;
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += p.getX();
                ++count;
            }
            rSum /= count;
        } else if (isHorizontal) {
            thetaRadians = 0;
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += p.getY();
                ++count;
            }
            rSum /= count;
        } else {
            thetaRadians = Math.atan(yInterceptAndSlope[1]);
            double ct = Math.cos(thetaRadians);
            double st = Math.sin(thetaRadians);
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += (p.getX() * ct) + (p.getY() * st);
                ++count;
            }
            rSum /= count;
        }
        
        int t = (int)Math.round(thetaRadians*180./Math.PI);
        if (t < 0) {
            t += 360;
        }
        
        // ----- evaluate ------
        
        /*
        have y = m * x + k  where m is slope and k is yIntercept

        point (x0, y0)

        d = sqrt( 
            (((x0 + m*y0 - m*k)/(m*m + 1)) - x0)^2 +
            (((m*(x0 + m*y0 - m*k))/(m*m + 1)) + k - y0)^2
        )
        */
        
        float m = yInterceptAndSlope[1];
        float k = yInterceptAndSlope[0];
        float mn, stdv;
            
        double sum = 0;
        int count = 0;
        double[] d = new double[linePoints.size()];
        List<PairInt> points = new ArrayList<PairInt>(linePoints);
        for (PairInt p : points) {
            int x0 = p.getX();
            int y0 = p.getY();
            if (isVertical) {
                d[count] = x0 - rSum;
            } else if (isHorizontal) {
                d[count] = y0 - rSum;
            } else {
                /*
                m*x - y + yIntercept = 0
                */
                d[count] = Math.sqrt(Math.abs(m*x0 - y0 + k)/(m*m + 1));  
            }
            sum += d[count];
            count++;
        }
        sum /= (double)linePoints.size();
        mn = (float)sum;
        sum = 0;
        for (int i = 0; i < d.length; ++i) {
            double diff = d[i] - mn;
            sum += (diff * diff);
            count++;
        }
        stdv = (float)Math.sqrt(sum/((double)linePoints.size() - 1.));        
        
        if (stdv == 0) {
            return new float[]{t, (float)Math.abs(rSum), Math.abs(mn), stdv};
        }
        
        // TODO: might want to allow a tolerance to be passed in
        double limit = stdv;
        // if it's highly inclined line, tolerance has to allow for step width
        if ((Math.abs(AngleUtil.getAngleDifference(t, 22.5f)) < 12) || 
            (Math.abs(AngleUtil.getAngleDifference(t, 67.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 112.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 157.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 202.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 247.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 292.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 337.5f)) < 12)
            ) {
            //TODO:
            // this picks up the meeting of lines at a corner, so
            // may need revision
            limit = 2*stdv;
        }
        if (limit < 1) {
            limit = 1;
        }
        
        Set<Integer> remove = new HashSet<Integer>();
        for (int i = 0; i < d.length; ++i) {
            double diff = Math.abs(d[i] - mn);
            if (diff > limit) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        // --- recalc mean and stdev w/o outliers ----
  //TODO: consider recalculating the slope and intercept w/ thiel sen
        sum = 0;
        count = 0;
        for (int i = 0; i < d.length; ++i) {
            if (remove.contains(Integer.valueOf(i))) {
                linePoints.remove(points.get(i));
                continue;
            }
            sum += d[i];
            count++;
        }
        sum /= (double)count;
        mn = (float)sum;
        count = 0;
        sum = 0;
        for (int i = 0; i < d.length; ++i) {
            if (remove.contains(Integer.valueOf(i))) {
                continue;
            }
            double diff = Math.abs(d[i] - mn);
            sum += (diff * diff);
            count++;
        }
        stdv = (float)Math.sqrt(sum/((double)count - 1.));        
        return new float[]{t, (float)Math.abs(rSum), Math.abs(mn), stdv};
    }
    
    public class HoughTransformLines {
        
        private final Map<PairInt, PairInt> pixelToPolarCoordMap;
        
        private final List<Set<PairInt>> sortedLineGroups;
        
        public HoughTransformLines(Map<PairInt, PairInt> pixToTRMap,
            List<Set<PairInt>> sortedGroups) {
            this.pixelToPolarCoordMap = pixToTRMap;
            this.sortedLineGroups = sortedGroups;
        }
        
        public Map<PairInt, PairInt> getPixelToPolarCoordMap() {
            return pixelToPolarCoordMap;
        }
        
        public List<Set<PairInt>> getSortedLineGroups() {
            return sortedLineGroups;
        }
    }
    
    /**
     * runtime complexity is approx O(N_pix * lg_2(N_pix)).
     * @param sortedTRKeys
     * @param thetaRadiusPixMap
     * @param thetaTol
     * @param radiusTol
     * @return 
     */
    public HoughTransformLines createPixTRMapsFromSorted(List<PairInt> sortedTRKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, 
        int thetaTol, int radiusTol) {
        
        // using a DFS visitor pattern and a tolerance for contiguous neighbor
        // grouping to aggretate the (theta, radius) solutions of adjacent
        // points
        DFSSimilarThetaRadiusGroupsFinder groupFinder = new 
            DFSSimilarThetaRadiusGroupsFinder();
        
        boolean allowGaps = true;
        
        // note that thetaRadiusPixMap is altered by this method
        Map<PairInt, PairInt> pixToTRMap = groupFinder.findConnectedPointGroups(
            sortedTRKeys, thetaRadiusPixMap, thetaTol, radiusTol, allowGaps);
        
        List<Set<PairInt>> sortedGroups = groupFinder.getSortedGroupsOfPoints();
        
        HoughTransformLines htl = new HoughTransformLines(pixToTRMap, sortedGroups);
        
        return htl;
    }
    
}
