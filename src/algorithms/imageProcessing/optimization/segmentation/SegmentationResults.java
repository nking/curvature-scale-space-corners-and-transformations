package algorithms.imageProcessing.optimization.segmentation;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.KNearestNeighbors;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class SegmentationResults {
    
    private final int[] xCentroids;
    private final int[] yCentroids;
    private final int[] nPoints;
    private List<Set<PairInt>> perimeters;
    
    private int maxY = Integer.MIN_VALUE;
    private int maxX = Integer.MIN_VALUE;
    
    public SegmentationResults(List<Set<PairInt>> segmentedSets) {
        
        int n = segmentedSets.size();
        
        int yMax = Integer.MIN_VALUE;
        int xMax = Integer.MIN_VALUE;
        
        {//debug
            for (int i = 0; i < n; ++i) {
                Set<PairInt> set = segmentedSets.get(i);
                int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(set);            
                if (xMinMaxYMinMax[1] > xMax) {
                    xMax = xMinMaxYMinMax[1];
                }
                if (xMinMaxYMinMax[3] > yMax) {
                    yMax = xMinMaxYMinMax[3];
                }
            }
            /*
            if (n > 0) {
                Image img = new Image(xMax + 1, yMax + 1);
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeAlternatingColor(
                    img, segmentedSets, "seg_" + ts);
            }*/
        }
                
        xCentroids = new int[n];
        yCentroids = new int[n];
        nPoints = new int[n];
        perimeters = new ArrayList<Set<PairInt>>();
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = segmentedSets.get(i);
            
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            xCentroids[i] = (int)Math.round(xyCen[0]);
            
            yCentroids[i] = (int)Math.round(xyCen[1]);
            
            nPoints[i] = set.size();
        
            Set<PairInt> border = finder.extractBorder(set);
            
            perimeters.add(border);
            
            //int[]{xMin, xMax, yMin, yMax}
            int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(border);
            
            if (xMinMaxYMinMax[1] > xMax) {
                xMax = xMinMaxYMinMax[1];
            }
            if (xMinMaxYMinMax[3] > yMax) {
                yMax = xMinMaxYMinMax[3];
            }
        }
    }
    
    public int getNumberOfItems() {
        return xCentroids.length;
    }
    
    public List<Set<PairInt>> getPointsList() {
        return this.perimeters;
    }
    
    /**
     * see notes in BenchmarkMeasurer.java
     * 
     * @param expected
     * @return 
     */
    public double evaluate(SegmentationResults expected,
        int dMax) {
            
        if (perimeters.size() == 0) {
            return 0;
        }
        
        BenchmarkMeasurer measurer = new BenchmarkMeasurer();
        
        float fMeasure = measurer.evaluate(this, expected, dMax);
        
        return fMeasure;
    }
    
    public int sumNPerimeters() {

        int n = 0;
        
        for (Set<PairInt> perimeter : perimeters) {
            n += perimeter.size();
        }
        
        return n;
    }    

    public Set<PairInt> getAllPoints() {

        Set<PairInt> output = new HashSet<PairInt>();
        
        for (Set<PairInt> set : perimeters) {
            output.addAll(set);
        }
        
        return output;
    }
    
    public List<Set<PairInt>> getPerimeters() {
        return perimeters;
    }
    
    public KNearestNeighbors createKNN() {
        
        int n = sumNPerimeters();
        
        if (n == 0) {
            throw new IllegalStateException("perimeters "
                + "cannot be empty");
        }
        
        float[] x = new float[n];
        float[] y = new float[n];
        int count = 0;
        for (Set<PairInt> perimeter : perimeters) {
            for (PairInt p : perimeter) {
                x[count] = p.getX();
                y[count] = p.getY();
                count++;
            }
        }
        KNearestNeighbors kNN = new KNearestNeighbors(x, y);
        return kNN;
    }

}
