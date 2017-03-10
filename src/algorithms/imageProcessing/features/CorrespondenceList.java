package algorithms.imageProcessing.features;

import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author nichole
 */
public class CorrespondenceList {
    
    private final float scale;
    private final int rotationInDegrees;
    private final int translationX;
    private final int translationY;
    private final int rangeRotation;
    private final int rangeTranslationX;
    private final int rangeTranslationY;
    
    private final List<PairInt> points1;
    private final List<PairInt> points2;
    
    private final TDoubleList cost12 = new TDoubleArrayList();
    
    //TODO: refactor so matches1 and 2 are constructed in class, not injected
    public CorrespondenceList(float scale, int rotationInDegrees, 
        int translationX, int translationY, int rangeRotation,
        int rangeTranslationX, int rangeTranslationY, 
        List<PairInt> matched1, List<PairInt> matched2) {
        
        this.scale = scale;
        this.rotationInDegrees = rotationInDegrees;
        this.translationX = translationX;
        this.translationY = translationY;
        
        this.rangeRotation = rangeRotation;
        this.rangeTranslationX = rangeTranslationX;
        this.rangeTranslationY = rangeTranslationY;
        
        this.points1 = matched1;
        this.points2 = matched2;
    }
    
    public CorrespondenceList(TransformationParameters params, 
        PairInt[] matched1, PairInt[] matched2) {
        
        this.scale = params.getScale();
        this.rotationInDegrees = Math.round(params.getRotationInDegrees());
        this.translationX = Math.round(params.getTranslationX());
        this.translationY = Math.round(params.getTranslationY());
        
        this.rangeRotation = Integer.MAX_VALUE;
        this.rangeTranslationX = Integer.MAX_VALUE;
        this.rangeTranslationY = Integer.MAX_VALUE;
        
        this.points1 = new ArrayList<PairInt>();
        this.points2 = new ArrayList<PairInt>();
        for (int i = 0; i < matched1.length; ++i) {
            points1.add(matched1[i]);
            points2.add(matched2[i]);
        }
    }
    
    public CorrespondenceList(List<QuadInt> matched12) {
        
        this.scale = Float.POSITIVE_INFINITY;
        this.rotationInDegrees = Integer.MIN_VALUE;
        this.translationX = Integer.MIN_VALUE;
        this.translationY = Integer.MIN_VALUE;
        this.rangeRotation = Integer.MAX_VALUE;
        this.rangeTranslationX = Integer.MAX_VALUE;
        this.rangeTranslationY = Integer.MAX_VALUE;
        
        this.points1 = new ArrayList<PairInt>();
        this.points2 = new ArrayList<PairInt>();
        for (int i = 0; i < matched12.size(); ++i) {
            QuadInt q = matched12.get(i);
            points1.add(new PairInt(q.getA(), q.getB()));
            points2.add(new PairInt(q.getC(), q.getD()));
        }
    }
    
    public void addMatch(PairInt p1, PairInt p2, double cost) {
        points1.add(p1);
        points2.add(p2);
        cost12.add(cost);
    }
    
    /*
     * sort point pairs by increasing cost.  NOTE that until the class
     * is refactored, it is up to the user to ensure that the numbers
     * of point and costs are the same.
    public void sortByAscendingCost() {
        
        if (points1.size() != points2.size() || points1.size() != cost12.size()) {
            throw new IllegalStateException("points lists and costs list must be same size");
        }
        int n = points1.size();
        int[] indexes = new int[n];
        for (int i = 0; i < n; ++i) {
            indexes[i] = i;
        }
        
        QuickSort.sortBy1stArg(cost12, indexes);
        
        List<PairInt> m1 = new ArrayList<PairInt>(n);
        List<PairInt> m2 = new ArrayList<PairInt>(n);
        
        for (int i = 0; i < indexes.length; ++i) {
            int idx = indexes[i];
            m1.add(points1.get(idx));
            m2.add(points2.get(idx));
        }
    }
    */
    
    /**
     * @return the scale
     */
    public float getScale() {
        return scale;
    }

    /**
     * @return the rotationInDegrees
     */
    public int getRotationInDegrees() {
        return rotationInDegrees;
    }

    /**
     * @return the translationX
     */
    public int getTranslationX() {
        return translationX;
    }

    /**
     * @return the translationY
     */
    public int getTranslationY() {
        return translationY;
    }

    /**
     * @return the rangeRotation
     */
    public int getRangeRotation() {
        return rangeRotation;
    }

    /**
     * @return the rangeTranslationX
     */
    public int getRangeTranslationX() {
        return rangeTranslationX;
    }

    /**
     * @return the rangeTranslationY
     */
    public int getRangeTranslationY() {
        return rangeTranslationY;
    }

    /**
     * @return the points1
     */
    public List<PairInt> getPoints1() {
        return points1;
    }

    /**
     * @return the points2
     */
    public List<PairInt> getPoints2() {
        return points2;
    }
   
    public TransformationParameters getParameters() {
        
        TransformationParameters params =
            new TransformationParameters();
        params.setRotationInDegrees(this.rotationInDegrees);
        params.setScale(this.scale);
        params.setTranslationX(this.translationX);
        params.setTranslationY(this.translationY);
        
        return params;
    }
}
