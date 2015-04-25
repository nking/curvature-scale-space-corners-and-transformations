package algorithms.imageProcessing.optimization;

/**
 *
 * @author nichole
 */
public class SkylineFits implements Comparable<SkylineFits> {

    protected ANDedClauses[] clauses;
    protected SetComparisonResults results;

    public SetComparisonResults getResults() {
        return results;
    }
    
    public ANDedClauses[] getClauses() {
        return clauses;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof SkylineFits)) {
            return false;
        }
        return this.results.equals(((SkylineFits) obj).getResults());
    }

    @Override
    public int compareTo(SkylineFits other) {
        return this.results.compareTo(other.results);
    }
}
