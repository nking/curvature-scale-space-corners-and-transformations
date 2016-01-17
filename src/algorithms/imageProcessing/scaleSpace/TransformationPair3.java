package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.features.BlobPerimeterRegion;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class TransformationPair3 {
    
    private final int startIdx1;
    
    private final int startIdx2;
    
    private final int nMaxMatchable;
        
    private double cost = Double.MAX_VALUE;
    
    private List<BlobPerimeterRegion> matchedContourRegions1 = new ArrayList<BlobPerimeterRegion>();
    
    private List<BlobPerimeterRegion> matchedContourRegions2 = new ArrayList<BlobPerimeterRegion>();
    
    private List<FeatureComparisonStat> matchedCompStats = new ArrayList<FeatureComparisonStat>();
    
    private int contourIndex1 = -1;
    
    private int contourIndex2 = -1;
           
    public TransformationPair3(int startRegionIndex1, int startRegionIndex2,
        int size) {
        
        this.startIdx1 = startRegionIndex1;
        
        this.startIdx2 = startRegionIndex2;
        
        this.nMaxMatchable = size;
    }

    /**
     * @return the startRegionIndex1
     */
    public int getStartRegionIndex1() {
        return startIdx1;
    }

    /**
     * @return the startRegionIndex2
     */
    public int getStartRegionIndex2() {
        return startIdx2;
    }

    public void addMatched(BlobPerimeterRegion region1, BlobPerimeterRegion region2,
        FeatureComparisonStat stat) {
        
        matchedContourRegions1.add(region1);
        matchedContourRegions2.add(region2);
        matchedCompStats.add(stat);
    }

    /**
     * @return the cost
     */
    public double getCost() {
        return cost;
    }

    /**
     * @param cost the cost to set
     */
    public void setCost(double cost) {
        this.cost = cost;
    }

    /**
     * @return the matchedContourRegions1
     */
    public List<BlobPerimeterRegion> getMatchedContourRegions1() {
        return matchedContourRegions1;
    }

    /**
     * @return the matchedContourRegions2
     */
    public List<BlobPerimeterRegion> getMatchedContourRegions2() {
        return matchedContourRegions2;
    }

    /**
     * @return the matchedCompStats
     */
    public List<FeatureComparisonStat> getMatchedCompStats() {
        return matchedCompStats;
    }

    public void setContourIndex1(int idx) {
        this.contourIndex1 = idx;
    }

    public void setContourIndex2(int idx) {
        this.contourIndex2 = idx;
    }
    
}
