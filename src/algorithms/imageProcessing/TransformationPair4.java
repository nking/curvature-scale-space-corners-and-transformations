package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class TransformationPair4 implements Comparable<TransformationPair4> {
    
    private final int startIdx1;
    
    private final int startIdx2;
            
    private double cost = Double.MAX_VALUE;
    
    private List<CornerRegion> matchedCornerRegions1 = new ArrayList<CornerRegion>();
    
    private List<CornerRegion> matchedCornerRegions2 = new ArrayList<CornerRegion>();
    
    private List<FeatureComparisonStat> matchedCompStats = new ArrayList<FeatureComparisonStat>();
    
    private int cornerListIndex1 = -1;
    
    private int cornerListIndex2 = -1;
           
    public TransformationPair4(int startRegionIndex1, int startRegionIndex2,
        int size) {
        
        this.startIdx1 = startRegionIndex1;
        
        this.startIdx2 = startRegionIndex2;        
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

    public void addMatched(CornerRegion region1, CornerRegion region2,
        FeatureComparisonStat stat) {
        
        matchedCornerRegions1.add(region1);
        matchedCornerRegions2.add(region2);
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
    public List<CornerRegion> getMatchedCornerRegions1() {
        return matchedCornerRegions1;
    }

    /**
     * @return the matchedContourRegions2
     */
    public List<CornerRegion> getMatchedCornerRegions2() {
        return matchedCornerRegions2;
    }

    /**
     * @return the matchedCompStats
     */
    public List<FeatureComparisonStat> getMatchedCompStats() {
        return matchedCompStats;
    }

    public void setCornerListIndex1(int idx) {
        this.cornerListIndex1 = idx;
    }

    public void setCornerListIndex2(int idx) {
        this.cornerListIndex2 = idx;
    }

    @Override
    public int compareTo(TransformationPair4 other) {
        
        int tN = matchedCompStats.size();
        double tCost = cost;
        
        int oN = other.getMatchedCompStats().size();
        double oCost = other.getCost();
        
        if (((tN >= oN) && (tCost < oCost)) || ((oCost == tCost) && (tN > oN))) {
            return -1;
        } else if ((tN == oN) && (oCost == tCost)) {
            return 0;
        } else {
            return 1;
        }
    }

    /**
     * @return the cornerListIndex1
     */
    public int getCornerListIndex1() {
        return cornerListIndex1;
    }

    /**
     * @return the cornerListIndex2
     */
    public int getCornerListIndex2() {
        return cornerListIndex2;
    }
    
}
