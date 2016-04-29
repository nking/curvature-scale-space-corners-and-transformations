package algorithms.search;

public class KDTreeNode {

	KDTreeNode right = null;
	KDTreeNode left = null;
	KDTreeNode parent = null;
	
	int x = -1;
	int y = -1;
	int key = -1; // median value
	int nChildren = 0;

    /**
     * @return the x
     */
    public int getX() {
        return x;
    }
    
    /**
     * @param theX
     */
    public void setX(int theX) {
        x = theX;
    }

    /**
     * @return the y
     */
    public int getY() {
        return y;
    }
    
    /**
     * @param theY
     */
    public void setY(int theY) {
        y = theY;
    }

    /**
     * @return the key
     */
    public int getKey() {
        return key;
    }
}
