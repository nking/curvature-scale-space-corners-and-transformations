package algorithms.search;

public class KDTreeNodeFloat {

	KDTreeNodeFloat right = null;
	KDTreeNodeFloat left = null;
	KDTreeNodeFloat parent = null;
	
	float x = -1;
	float y = -1;
	float key = -1; // median value
	int nChildren = 0;
    int depth = 0;
    
    /**
     * @return the x
     */
    public float getX() {
        return x;
    }
    
    /**
     * @param theX
     */
    public void setX(float theX) {
        x = theX;
    }

    /**
     * @return the y
     */
    public float getY() {
        return y;
    }
    
    /**
     * @param theY
     */
    public void setY(float theY) {
        y = theY;
    }

    /**
     * @return the key
     */
    public float getKey() {
        return key;
    }
}
