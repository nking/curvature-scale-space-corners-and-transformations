package algorithms.graphs;

import algorithms.QuickSort;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;

/**
 * customized DAG for the Watershed class to store the connections between
 * a pixel and it's lower intensity neighbors, ordered by steepness.
 * 
 * Unused functions have been removed, for example, the nodes do not store
 * incoming connections.
 * 
 * @author nichole
 */
public class CustomWatershedDAG {

    /*TODO: could consider using the image internal index here and converting
    the (x,y) coordinate into it to use an array if usage will always be
    on the entire image.  protected CustomWatershedNode[] vertices;
    
    Because Watershed is being designed for use on Sets of points rather than 
    the entire image, will use a Map instead.
    */
    protected final Map<PairInt, CustomWatershedNode> vertices;
    
    protected int nVertices = 0;

    public CustomWatershedDAG() {
        vertices = new HashMap<PairInt, CustomWatershedNode>();
    }
    
    public CustomWatershedDAG(int expectedNVertices) {
        vertices = new HashMap<PairInt, CustomWatershedNode>(expectedNVertices);
    }
    
    /**
     * insert into DAG for key, the list of nodes ordered such that the steepest
     * is at smallest indexes.  It's the invoker's responsibility to insure
     * that the order is correct.
     * @param key
     * @param orderedSLN 
     */
    public void insert(PairInt key, CustomWatershedNode orderedSLN) {
        
        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        if (orderedSLN == null) {
            throw new IllegalStateException("orderedSLN cannot be null");
        }
        
        vertices.put(key, orderedSLN);
    }
    
    /**
     * insert into DAG for key, the list of nodes after sorting them here such 
     * that the steepest is at smallest indexes.
     * @param key
     * @param intensityDifferences an array of lenUsable items holding the
     * difference in intensity of the pixel at key to the intensity of the
     * pixel located in points.  Note that the order of this array is 
     * altered here.
     * @param points an array of lenUsable items holding the neighbors of
     * key which had intensities lower than the key's intensity.  Note that the 
     * order of this array is altered here.
     * @param lenUsable the number of items in the arrays 
     * intensityDifferences and points to read and store.
     */
    public void orderAndInsert(PairInt key, int[] intensityDifferences, 
        PairInt[] points, int lenUsable) {
        
        if (key == null) {
            throw new IllegalArgumentException("key cannot be null");
        }
        if (intensityDifferences == null) {
            throw new IllegalArgumentException("intensityDifferences cannot be null");
        }
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (intensityDifferences.length < lenUsable) {
            throw new IllegalArgumentException(
            "intensityDifferences is smaller than lenUsable");
        }
        if (points.length < lenUsable) {
            throw new IllegalArgumentException("points is smaller than lenUsable");
        }
        
        CustomWatershedNode orderedSLN;
        
        if (lenUsable == 0) {
            orderedSLN = new CustomWatershedNode(key, 0);
        } else if (lenUsable == 1) {
            orderedSLN = new CustomWatershedNode(key, 1);
            orderedSLN.insertOutgoing(points[0]);
        } else {
            // ascending sort, so smallest are at index 0
            QuickSort.sortBy1stArg(intensityDifferences, points, 0, lenUsable - 1);
            orderedSLN = new CustomWatershedNode(key, lenUsable);
            for (int i = (lenUsable - 1); i > -1; --i) {
                PairInt p = points[i];
                orderedSLN.insertOutgoing(p);
            }
        }
                
        vertices.put(key, orderedSLN);
    }
    
    public boolean isEmpty() {
        return (nVertices == 0);
    }

    public void setToResolved(PairInt key, PairInt resolution) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        if (resolution == null) {
            throw new IllegalStateException("resolution cannot be null");
        }
        
        CustomWatershedNode node = vertices.get(key);

        if (node != null) {
            node.setToResolved(resolution);
        }
    }
    
    public PairInt getResolved(PairInt key) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        CustomWatershedNode node  = vertices.get(key);

        if (node != null) {
            return node.getResolved();
        }
        
        return null;
    }
    
    /**
     * return whether the node at key has been processed as resolved yet.
     * Note that an IllegalArgumentException is thrown if key is not in the dag.
     * @param key
     * @return 
     */
    public boolean isResolved(PairInt key) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        CustomWatershedNode node  = vertices.get(key);

        if (node == null) {
            throw new IllegalArgumentException("key was not found in dag");
        }
        
        return node.isResolved();
    }
    
    /**
     * returned the connected number for the node at key.  Note that an 
     * IllegalArgumentException is thrown if key is not in the dag.
     * @param key
     * @return 
     */
    public int getConnectedNumber(PairInt key) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        CustomWatershedNode node  = vertices.get(key);

        if (node == null) {
            throw new IllegalArgumentException("key was not found in dag");
        }
        
        return node.getConnectedNumber();
    }
    
    /**
     * returned the connected number for the node at key.  Note that an 
     * IllegalArgumentException is thrown if key is not in the dag.
     * @param key
     * @return 
     */
    public PairInt getConnectedNode(PairInt key, int nodeNumber) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        CustomWatershedNode node  = vertices.get(key);

        if (node == null) {
            throw new IllegalArgumentException("key was not found in dag");
        }
        
        return node.get(nodeNumber);
    }
    
    /**
     * returned the connected number for the node at key.  Note that an 
     * IllegalArgumentException is thrown if key is not in the dag.
     * @param key
     */
    public void resetConnectedNode(PairInt key, int nodeNumber, PairInt nodeValue) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        CustomWatershedNode node  = vertices.get(key);

        if (node == null) {
            throw new IllegalArgumentException("key was not found in dag");
        }
        
        node.reset(nodeNumber, nodeValue);
    }
    
    /**
     * returned the connected number for the node at key.  Note that an 
     * IllegalArgumentException is thrown if key is not in the dag.
     * @param key
     * @return 
     */
    public boolean contains(PairInt key) {

        if (key == null) {
            throw new IllegalStateException("key cannot be null");
        }
        
        return vertices.containsKey(key);
    }

}
