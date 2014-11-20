package algorithms.imageProcessing;
/**
 *
 * @author nichole
 */
public class TransformationParameters {
        
    //NOTE:  these should be 0's as default
    
    private float rotationInRadians = 0;
    
    private float translationX = 0;
    
    private float translationY = 0;
    
    private float scale = 1;

    /**
     * @return the rotation in units of degrees
     */
    public float getRotationInDegrees() {
        return (float)(rotationInRadians * 180./Math.PI);
    }

    /**
     * @param rotationInDegrees the rotation to set in units of degrees
     */
    public void setRotationInDegrees(float rotationInDegrees) {
        this.rotationInRadians = (float) (rotationInDegrees * Math.PI/180.);
    }
    
    /**
     * @param theRotationInRadians the rotation to set in units of radians
     */
    public void setRotationInRadians(float theRotationInRadians) {
        this.rotationInRadians = theRotationInRadians;
    }

    /**
     * @return the rotation in units of radians
     */
    public float getRotationInRadians() {
        return rotationInRadians;
    }
    
    /**
     * @return the translationX
     */
    public float getTranslationX() {
        return translationX;
    }

    /**
     * @param translationX the translationX to set
     */
    public void setTranslationX(float translationX) {
        this.translationX = translationX;
    }

    /**
     * @return the translationY
     */
    public float getTranslationY() {
        return translationY;
    }

    /**
     * @param translationY the translationY to set
     */
    public void setTranslationY(float translationY) {
        this.translationY = translationY;
    }

    /**
     * @return the scale
     */
    public float getScale() {
        return scale;
    }

    /**
     * @param scale the scale to set
     */
    public void setScale(float scale) {
        this.scale = scale;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("rotationInRadians=").append(rotationInRadians)
            .append("\nrotationInDegrees=").append(rotationInRadians*180./Math.PI)
            .append("\nscale=").append(scale)
            .append("\ntranslationX=").append(translationX)
            .append("\ntranslationY=").append(translationY);
        
        return sb.toString();
    }
    
}
