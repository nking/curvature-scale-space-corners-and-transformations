package algorithms.imageProcessing.features;

import algorithms.util.PairInt;
import java.util.List;

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
    
}
