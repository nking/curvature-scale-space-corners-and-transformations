package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.LinearRegression;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class LinesFinder {
    
    /**
     * array of x coordinates of points found as lines.  note that some of the
     * short segments are parts of curves, indistinguishable unless larger 
     * minimum line lengths are used to filter them out.
     */
    private TIntList xs = new TIntArrayList();
    
    /**
     * array of y coordinates of points found as lines.  note that some of the
     * short segments are parts of curves, indistinguishable unless larger 
     * minimum line lengths are used to filter them out.
     */
    private TIntList ys = new TIntArrayList();
    
    /**
     * map with key = unique segment index number, value =
     * indexes of the xs and ys arrays which are a contiguous line segment.
     */
    private TIntObjectMap<TIntList> segmentIndexes = new TIntObjectHashMap<TIntList>();
    
    /**
     * map with key = unique index number, value = the polar theta in degrees
     * and distance from origin for the line segment.
     */
    private TIntObjectMap<PairInt> segmentTRMap = new TIntObjectHashMap<PairInt>();
   
    /**
     * map with key = the polar theta in degrees
     * and distance from origin for the line segment,
     * value = segment indexes of lines with this combination of
     * theta and radius.
     */
    private Map<PairInt, TIntList> trSegmentIndexesMap = new HashMap<PairInt, TIntList>();
    
    private List<PairInt> orderedTRList = null;
    private List<TIntList> orderedTRXYIndexes = null;
    
    private int lastSegIdx = -1;
    
    private int lastCol = -1;
    private int lastRow = -1;
    
    private int minLineLength = 20;
    
    private boolean debug = false;
    
    private float thresh = (float)(1.e-7);
    
    /**
     * if this is set, vertical lines found at polar radius 0 and width from
     * origin are removed and so are horizontal lines found at polar radius
     * and height from origin.
     */
    public void setToRemoveBorderLines(int lastColumn, int lastRow) {
        this.lastCol = lastColumn;
        this.lastRow = lastRow;
    }
    
    /**
     * override the default minimum line length of 30 to the given value
     * @param length 
     */
    public void overrideMinimumLineLength(int length) {
        this.minLineLength = length;
    }
    
    /**
     * override default minimum length of 20, to this value
     * @param length 
     */
    public void overrideMinimumLength(int length) {
        minLineLength = length;
    }
   
    /**
     * override the default threshold of 1.e7.  the absolute average differernce
     * of chords over a window size must be less than the threshold in order
     * for the segment to be considered a line.  Increasing this number
     * too much will possibly result in including curves.
     * @param threshold 
     */
    public void overrideThreshold(float threshold) {
        this.thresh = threshold;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void find(List<Set<PairInt>> listOfContigousLabels) {
                
        // -- extract the boundaries of the sets
        // -- find the lines around each boundary using shape fitting
        // -- store results as member variables
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<PairIntArray> listOfBounds = new ArrayList<PairIntArray>();
        
        // extract ordered bounds
        PerimeterFinder2 pFinder = new PerimeterFinder2();
        List<PairIntArray> extractedBounds = new ArrayList<PairIntArray>();
        for (int i = 0; i < listOfContigousLabels.size(); ++i) {
            Set<PairInt> set = listOfContigousLabels.get(i);
            if (set.size() < 3) {
                continue;
            }
         
            PairIntArray b = pFinder.extractOrderedBorder(
                new HashSet<PairInt>(set));
            if (b != null && b.getN() > 2) {
                listOfBounds.add(b);            
            }
        }
        
        assert(assertUniquePoints(listOfBounds));
        
Image dbg = new Image(256, 192);
ImageIOHelper.addAlternatingColorCurvesToImage(
listOfBounds, dbg, 0);
MiscDebug.writeImage(dbg, "_boundaries_");

        find1(listOfBounds);
    }
    
    public void find1(List<PairIntArray> listOfOrderedBounds) {
        
        // -- find the lines around each boundary using shape fitting
        // -- store results as lists of x, y 
        //    and map with key=segIdx, value=list of x,y indexes
        //    and map with key-segIdx, value = theta (float)
        //    and map with key=segIdx, value = polar radius (int)
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                
        // extract ordered bounds
        PerimeterFinder2 pFinder = new PerimeterFinder2();
        List<PairIntArray> extractedBounds = new ArrayList<PairIntArray>();
        for (int i = 0; i < listOfOrderedBounds.size(); ++i) {
            
            PairIntArray b = listOfOrderedBounds.get(i);
            if (b == null || b.getN() < 3) {
                continue;
            }
        
            LineFinder matcher = new LineFinder();
            if (debug) {
                matcher.setToDebug();
            }
            if (lastCol > -1) {
                matcher.setToRemoveBorderLines(lastCol, lastRow);
            }
            //matcher.overrideMinimumLineLength(minLineLength);
            //matcher._overrideToThreshhold(thresh);
            LineFinder.LineResult r = matcher.match(b);
            List<PairInt> lr = r.getLineIndexRanges();
            //System.out.println("reading nRanges=" + lr.size());
            
            for (int ii = 0; ii < lr.size(); ++ii) {
                int startIdx = lr.get(ii).getX(); 
                int stopIdx = lr.get(ii).getY(); 
                
                if (debug) {
                    System.out.println("indexes: " + startIdx + ":" + stopIdx 
                        + "   " + " segIdx=" + ii +
                        String.format(" (%d,%d) to (%d,%d) ",
                        b.getX(startIdx), b.getY(startIdx),
                        b.getX(stopIdx), b.getY(stopIdx))
                    );
                }
                
                int lineX0 = b.getX(startIdx);
                int lineY0 = b.getY(startIdx);
                int lineX1 = b.getX(stopIdx);
                int lineY1 = b.getY(stopIdx);

                if (debug) {
                    System.out.println("coords: (" + lineX0 + "," + lineY0 + ") "
                        + " (" + lineX1 + "," + lineY1 + ") ");
                }
                
                double polarR = curveHelper.distanceFromPointToALine(
                    lineX0, lineY0, lineX1, lineY1, 0, 0);
                int radius = (int)Math.round(polarR);
                
                // don't store lines on image boundaries if this is set
                if (lastCol > -1) {
                    if (radius < 5) {
                        continue;
                    }
                }
                
                // -180 to 180 are the ranges.  also need an offset by 90 perpendicular
                double theta = Math.atan2(lineY1 - lineY0, lineX1 - lineX0);
                int thetaDeg = (int)Math.round(theta * 180./Math.PI);
                // perpendicular to slipe:
                thetaDeg -= 90;
                // correction to place within range 0 to 180:
                if (thetaDeg < 0) {
                    while (thetaDeg < 0) {
                        // reverse the line direction
                        thetaDeg += 180;
                    }
                } else if (thetaDeg == 180) {
                    thetaDeg = 0;
                }
                
                if (debug) {
                    System.out.println("theta, radius: (" + thetaDeg + "," 
                        + radius + ")");
                }

                // don't store lines on image bundaries if this is set               
                if (lastCol > -1) {
                    if ((thetaDeg == 0 || thetaDeg == 180) && 
                        (radius > (lastCol - 5))) {
                        continue;
                    } else if ((thetaDeg == 90 || thetaDeg == 270) 
                        && (radius > (lastRow - 5))) {
                        continue;
                    }
                }
                
                lastSegIdx++;
                
                PairInt tr = new PairInt(thetaDeg, radius);
                segmentTRMap.put(lastSegIdx, tr);
                
                if (debug) {
                    System.out.println("*coords: (" + lineX0 + "," + lineY0 + ") "
                        + " (" + lineX1 + "," + lineY1 + ") ");
                }
                
                int xsIdx;
                TIntList idxs = new TIntArrayList();
                for (int j = startIdx; j <= stopIdx; ++j) {
                    int x = b.getX(j);
                    int y = b.getY(j);
                    xsIdx = xs.size();
                    xs.add(x);
                    ys.add(y);
                    idxs.add(xsIdx);
                }
                segmentIndexes.put(lastSegIdx, idxs);
                
                TIntList segIdxs = trSegmentIndexesMap.get(tr);
                if (segIdxs == null) {
                    segIdxs = new TIntArrayList();
                    trSegmentIndexesMap.put(tr, segIdxs);
                }
                segIdxs.add(lastSegIdx);
            }
        }
    }
    
    /**
     * for best results, consider the deltae2000 gradient
     * which is used in ImageSegmentation.objectSegmentation
     * @param gradient 
     */
    public void correctLinesWithGradient(GreyscaleImage gradient) {
         
        if (orderedTRList == null) {
            groupWithinTolerance();
        }
        
        MiscDebug.writeImage(gradient, "_gradient_");
    
        // for each line 
        //    extract a region + and minus 10 pixels
        //    to a side in the gradient image and fit a
        //    line to those pixels.
        //    can use thiel sen, but there are probably
        //    faster algorithms which also remove outliers.
        //    -- can use the intensities for weights, but the
        //       deltae gradient largely has same intensity
        //       in most pixels...
        
        boolean changed = false;
        
        for (int i = 0; i < orderedTRList.size(); ++i) {
            PairInt tr = orderedTRList.get(i);
            TIntList xyIdxs = orderedTRXYIndexes.get(i);
            //82, 48
            int[] xyMinMax = findXYBounds(xyIdxs);
            //System.out.println("syminmax=" + Arrays.toString(xyMinMax));    
            //System.out.println("line: " + tr);
            double thetaR = tr.getX() * Math.PI/180.;
            // extract bounds += 10 pix from gradient
            // perpendicular to endpoints
            int d = 8;
            int dX = (int)Math.round(d * Math.cos(thetaR));
            int dY = (int)Math.round(d * Math.sin(thetaR));
            int x0 = xyMinMax[0] - dX;
            if (x0 < 0) { x0 = 0;}
            int x1 = xyMinMax[1] + dX;
            if (x1 > (xyMinMax[1] - 1)) { x1 = xyMinMax[1] - 1;}
            int y0 = xyMinMax[2] - dY;
            if (y0 < 0) { y0 = 0;}
            int y1 = xyMinMax[3] + dY;
            if (y1 > (xyMinMax[3] - 1)) { y1 = xyMinMax[3] - 1;}

            TFloatList xG = new TFloatArrayList();
            TFloatList yG = new TFloatArrayList();

            for (int x = x0; x <= x1; ++x) {
                for (int y = y0; y <= y1; ++y) {
                    int v = gradient.getValue(x, y);
                    if (v > 0) {
                        yG.add(y);
                        xG.add(x);
                    }
                }
            }

            if (xG.size() < 4) {
                continue;
            }

            LinearRegression lr = new LinearRegression();
            float[] yinterceptSlope 
                = lr.calculateTheilSenEstimatorParams(
            //lr.plotTheLinearRegression(
                xG.toArray(new float[xG.size()]), 
                yG.toArray(new float[yG.size()]));

            System.out.println("yin and slope=" + Arrays.toString(yinterceptSlope));
            int t = Math.round(yinterceptSlope[1]) - 90;
            if (t < 0) {
                t += 180;
            }
            int r = Math.round(yinterceptSlope[0]);

            float dt = Math.abs(t - tr.getX());
            float dr = Math.abs(r - tr.getY());
            
            if ((dt > 25) || (dr > 25)) {
                continue;
            }

            if ((dt > 0.5) || (dr > 0.5)) {
                changed = true;
            }
            
            orderedTRList.set(i, new PairInt(t, r));
        }
        
        if (changed) {
            // assuming not many lines and using
            // an O(N^2) approach for the merging for now
            for (int i = (orderedTRList.size() - 1); i > -1; --i) {
                PairInt tr = orderedTRList.get(i);
                for (int j = 0; j < (i - 1); ++j) {
                    PairInt tr2 = orderedTRList.get(j);
                    if (tr2.equals(tr)) {
                        TIntList xyi = 
                            orderedTRXYIndexes.get(i);
                        orderedTRXYIndexes.get(j)
                            .addAll(xyi);
                        orderedTRList.remove(i);
                        orderedTRXYIndexes.remove(i);
                        break;
                    }
                }
            }
            // sort
            sortOrderedLists();
        }
    }
    
    private void sortOrderedLists() {
        
        if (orderedTRList == null) {
        
            groupWithinTolerance();
        
        } else {
            
            int[] indexes = new int[orderedTRList.size()];
            int[] np = new int[indexes.length];
            
            for (int i = 0; i < np.length; ++i) {
                indexes[i] = i;
                np[i] = orderedTRXYIndexes.get(i).size();
            }
            QuickSort.sortBy1stArg(np, indexes);
            
            List<PairInt> otl = new ArrayList<PairInt>();
            List<TIntList> otiL = new ArrayList<TIntList>();
            
            for (int i = (np.length - 1); i > -1; --i) {
                int idx = indexes[i];
                otl.add(orderedTRList.get(idx));
                otiL.add(orderedTRXYIndexes.get(idx));
            }
            orderedTRList = otl;
            orderedTRXYIndexes = otiL;
        }
    }
     
    /**
     * combine the entries for a theta and radius within theta tolerance
     * and radius tolerance into lists.
     */
    public void groupWithinTolerance() {
        
        orderedTRList = new ArrayList<PairInt>();
        orderedTRXYIndexes = new ArrayList<TIntList>();
    
        Map<PairInt, TIntList> trXYIndexesMap = new
            HashMap<PairInt, TIntList>();
        
        int n = trSegmentIndexesMap.size();
        PairInt[] trs = new PairInt[n];
        int[] nPoints = new int[n];
        int[] lIdxs = new int[n];
        
        int count = 0;
        for (Entry<PairInt, TIntList> entry : trSegmentIndexesMap.entrySet()) {
            PairInt tr = entry.getKey();            
            TIntList segIdxs = entry.getValue();
        
            int np = 0;
            for (int j = 0; j < segIdxs.size(); ++j) {
                int segIdx = segIdxs.get(j);
                TIntList idxs = segmentIndexes.get(segIdx);
                np += idxs.size();
                
                TIntList a = trXYIndexesMap.get(tr);
                if (a == null) {
                    a = new TIntArrayList();
                    trXYIndexesMap.put(tr, a);
                }
                a.addAll(idxs);
              
                /*
                for (int k = 0; k < idxs.size(); ++k) {
                    int idx = idxs.get(k);
                    System.out.println(
                        String.format("-- (%d,%d) ",
                        xs.get(idx), ys.get(idx))
                    );
                }
                */
            }
            trs[count] = tr;
            nPoints[count] = np;
            lIdxs[count] = count;
            count++;
        }
        QuickSort.sortBy1stArg(nPoints, lIdxs);
                
        Set<PairInt> skip = new HashSet<PairInt>();
        
        for (int i = (count - 1); i > -1; --i) {
            int lIdx = lIdxs[i];
            int np = nPoints[i];
            
            PairInt tr = trs[lIdx];
            
            if (skip.contains(tr)) {
                continue;
            }

            TIntList xyList = new TIntArrayList();

            int sumT = 0;
            int sumR = 0;
            int sumN = 0;
           
            for (int t0 = tr.getX() - 2; t0 <= tr.getX() + 2; ++t0) {
                for (int r0 = tr.getY() - 2; r0 <= tr.getY() + 2; ++r0) {
                    PairInt tr0 = new PairInt(t0, r0);
                    TIntList a = trXYIndexesMap.get(tr0);
                    if (a == null) {
                        continue;
                    }
                    
                    /*for (int k = 0; k < a.size(); ++k) {
                        int idx = a.get(k);
                        System.out.println(
                        String.format("-- (%d,%d) tr=%s",
                        xs.get(idx), ys.get(idx), tr0)
                        );
                    }*/
                    
                    xyList.addAll(a);
                    trXYIndexesMap.remove(tr0);
                    skip.add(tr0);
                    sumT += tr0.getX();
                    sumR += tr0.getY();
                    sumN++;
                }
            }
            sumR /= sumN;
            sumT /= sumN;
            PairInt tr0 = new PairInt(sumT, sumR);
            orderedTRList.add(tr0);
            orderedTRXYIndexes.add(xyList);
        }
    }
    
    public void debugDraw(algorithms.imageProcessing.Image img) {
        
        if (orderedTRList == null) {
            throw new IllegalStateException("groupWithinTolerance must be "
                + " invoked first");
        }
        
        sortOrderedLists();
        
        int w = img.getWidth();
        int h = img.getHeight();
          
        boolean drawLines = true;
             
        // TODO: consider revising this by nPoints
        int end = 10;
        if (end > (orderedTRList.size() - 1)) {
            end = orderedTRList.size();
        }
        for (int i = 0; i < end; ++i) {
            
            PairInt tr = orderedTRList.get(i);
            TIntList xyIdxs = orderedTRXYIndexes.get(i);
            int np = xyIdxs.size();
            
            int clr = ImageIOHelper.getNextColorRGB(i);
            
            System.out.println("np=" + np + " tr=" + tr);
            
            if (drawLines) {
                
                int[] eps = LinesAndAngles.calcPolarLineEndPoints(
                    tr.getX(), tr.getY(), img.getWidth(), img.getHeight());

                //System.out.println("tr=" + tr.toString() + " eps=" +
                //    Arrays.toString(eps) + " w=" + img.getWidth() + 
                //    " h=" + img.getHeight());
               
                ImageIOHelper.drawLineInImage(
                    eps[0], eps[1], eps[2], eps[3], img, 1, 
                    clr);
                
            } else {                
                for (int k = 0; k < xyIdxs.size(); ++k) {
                    int idx = xyIdxs.get(k);
                    int x = xs.get(idx);
                    int y = ys.get(idx);
                    ImageIOHelper.addPointToImage(x, y, img, 1, 
                        clr);
                    
                }
                //System.out.println("  tr=" + tr + " n=" 
                //    + xyIdxs.size() + " i=" + i + " clr=" + clr);
            }
        }
    }

    private int[] findXYBounds(TIntList xyIdxs) {

        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (int i = 0; i < xyIdxs.size(); ++i) {
            int idx = xyIdxs.get(i);
            int x = xs.get(idx);
            int y = ys.get(idx);
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        
        return new int[]{xMin, xMax, yMin, yMax};
    }

    private boolean assertUniquePoints(
        List<PairIntArray> listOfBounds) {

        Set<PairInt> exists = new HashSet<PairInt>();
        for (int i = 0; i < listOfBounds.size(); ++i) {
            PairIntArray a = listOfBounds.get(i);
            for (int j = 0; j < a.getN(); ++j) {
                PairInt p = new PairInt(a.getX(j), a.getY(j));
                if (exists.contains(p)) {
                    return false;
                }
                exists.add(p);
            }
        }
        
        return true;
    }
}
