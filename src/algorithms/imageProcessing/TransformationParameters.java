package algorithms.imageProcessing;
/**
 *
 * @author nichole
 */
public class TransformationParameters {
            
    private float rotationInRadians = 0;
    
    private float translationX = 0;
    
    private float translationY = 0;
    
    private float scale = 1;

    /* the projective parameters representation may change.
    for now, using the denominator homography*/
    private float p20 = 0;
    private float p21 = 0;
    private float p22 = 1;
    
    public void setP20(float h31) {
        p20 = h31;
    }
    public void setP21(float h32) {
        p21 = h32;
    }
    public void setP22(float h33) {
        p22 = h33;
    }
    public float getP20() {
        return p20;
    }
    public float getP21() {
        return p21;
    }
    public float getP22() {
        return p22;
    }
    
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
            .append(" rotationInDegrees=").append(rotationInRadians*180./Math.PI)
            .append(" scale=").append(scale)
            .append(" translationX=").append(translationX)
            .append(" translationY=").append(translationY);
        
        return sb.toString();
    }

}
