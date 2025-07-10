package algorithms.imageProcessing.transform;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class TransformationParameters {
            
    private float rotationInRadians = 0;
    
    private float translationX = 0;
    
    private float translationY = 0;
    
    private float scale = 1;
    
    private float originX = 0;
    
    private float originY = 0;
    
    /**
     * the number of points that went into the solution
     */
    private int n = -1;
    
    private float[] standardDeviationsScaleRotTransXY = null;

    /**
     * @return the rotation in units of degrees
     */
    public float getRotationInDegrees() {
        return (float)(rotationInRadians * 180./Math.PI);
    }
    
    public void setStandardDeviations(float[] stDevsScaleRotTransXY) {
        standardDeviationsScaleRotTransXY = new float[4];
        System.arraycopy(stDevsScaleRotTransXY, 0, 
            standardDeviationsScaleRotTransXY, 0, stDevsScaleRotTransXY.length);
    }
    
    public float[] getStandardDeviations() {
        return standardDeviationsScaleRotTransXY;
    }
    
    public void setNumberOfPointsUsed(int numberOfPoints) {
        n = numberOfPoints;
    }
    
    public int getNumberOfPointsUsed() {
        return n;
    }

    /**
     * @param rotationInDegrees the rotation to set in units of degrees
     */
    public void setRotationInDegrees(float rotationInDegrees) {
        if (rotationInDegrees < 0) {
            rotationInDegrees += 360;
        } else if (rotationInDegrees > 360) {
            int m = (int)(rotationInDegrees/360) * 360;
            rotationInDegrees -= m;
        }
        this.rotationInRadians = (float) (rotationInDegrees * Math.PI/180.);
    }
    
    /**
     * @param theRotationInRadians the rotation to set in units of radians
     */
    public void setRotationInRadians(float theRotationInRadians) {
        float twoPi = (float)(2. * Math.PI);
        if (theRotationInRadians < 0) {
            theRotationInRadians += twoPi;
        } else if (theRotationInRadians > twoPi) {
            float m = (int)(theRotationInRadians/twoPi) * twoPi;
            theRotationInRadians -= m;
        }
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
    
    public void setOriginX(float theXOrigin) {
        originX = theXOrigin;
    }
    
    public void setOriginY(float theYOrigin) {
        originY = theYOrigin;
    }
   
    public float getOriginX() {
        return originX;
    }
    
    public float getOriginY() {
        return originY;
    }
    
    public TransformationParameters copy() {
        
        TransformationParameters cp = new TransformationParameters();
        cp.setRotationInRadians(rotationInRadians);
        cp.setTranslationX(translationX);
        cp.setTranslationY(translationY);
        cp.setScale(scale);
        cp.setOriginX(originX);
        cp.setOriginY(originY);
        cp.setNumberOfPointsUsed(n);
        if (standardDeviationsScaleRotTransXY != null) {
            cp.setStandardDeviations(standardDeviationsScaleRotTransXY);
        }
        
        return cp;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof TransformationParameters)) {
            return false;
        }
        
        TransformationParameters other = (TransformationParameters)obj;
        
        boolean paramsEqual = ((other.rotationInRadians == rotationInRadians) &&
            (other.translationX == translationX) && 
            (other.translationY == translationY) &&
            (other.scale == scale) &&
            (other.originX == originX) && (other.originY == originY)
            );
        
        if (!paramsEqual) {
            return false;
        }
                
        if ((standardDeviationsScaleRotTransXY != null) && 
            (other.getStandardDeviations() == null)) {
            return false;
        } else if ((standardDeviationsScaleRotTransXY == null) && 
            (other.getStandardDeviations() != null)) {
            return false;
        } else if ((standardDeviationsScaleRotTransXY != null) && 
            (other.getStandardDeviations() != null)) {
            boolean sEquals = Arrays.equals(standardDeviationsScaleRotTransXY, 
                other.getStandardDeviations());
            if (!sEquals) {
                return false;
            }
        }

        //TODO: revist hashCode and equals w.r.t. standardDeviationsScaleRotTransXY and n
        
        return (n == other.getNumberOfPointsUsed());
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 31 * hash + Float.floatToIntBits(this.rotationInRadians);
        hash = 31 * hash + Float.floatToIntBits(this.translationX);
        hash = 31 * hash + Float.floatToIntBits(this.translationY);
        hash = 31 * hash + Float.floatToIntBits(this.scale);
        hash = 31 * hash + Float.floatToIntBits(this.originX);
        hash = 31 * hash + Float.floatToIntBits(this.originY);
        return hash;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("rotationInRadians=").append(rotationInRadians)
            .append(" rotationInDegrees=").append(rotationInRadians*180./Math.PI)
            .append(" scale=").append(scale)
            .append("\ntranslationX=").append(translationX)
            .append(" translationY=").append(translationY)
            .append(" originX=").append(originX)
            .append(" originY=").append(originY);
        
        if (standardDeviationsScaleRotTransXY != null) {
            sb.append("\nstDevScale=").append(standardDeviationsScaleRotTransXY[0]);
            sb.append(" stDevRot=").append(standardDeviationsScaleRotTransXY[1]);
            sb.append(" stDevTransX=").append(standardDeviationsScaleRotTransXY[2]);
            sb.append(" stDevTransY=").append(standardDeviationsScaleRotTransXY[3]);
        }
        
        if (n != -1) {
            sb.append("\nn=").append(n);
        }
        
        return sb.toString();
    }

}
