package algorithms.search;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 * k-dimension tree is a binary tree used to store coordinates for quick nearest
 * neighbor or range searches.
 *    -- values are always stored in leaves
 *    -- the meaning of an internal node depends upon the depth of the
 *       node within the tree.
 *       
 * Note that learning the true medium while constructing the tree takes more
 * time, but leads to a better balanced tree.
 * 
 * adapted from pseudocode from
 * http://ldots.org/kdtree which licenses the content as:
 * http://creativecommons.org/licenses/by-sa/2.0/
 * 
 * 
 *  Memory requirements:
 *    
 *        array arguments allocated elsewhere and modifed in this class (but not copied):
 *            x:  8 bytes for a reference to array location, then 4*N bytes for items
 *            y:  8 bytes for a reference to array location, then 4*N bytes for items
 *            arrays total = 2 * (8 + 4*N)
 *                         = 8N + 16 bytes   where N is the number of points (point is one x,y pair)
 *                         
 *        KDTreeNodes:  16 bytes for each object + 2 int = 8 bytes = 24 bytes.
 *                      * Each node has 3 references = 24 bytes, for total = 48 bytes
 *                      * A leaf will have one reference = 8 bytes, for total = 32 bytes.
 *                      * A root will have 2 references = 16 bytes, for total = 40 bytes.
 *                      (the total for a class has to round up to an 8 byte multiple too if needed.)
 *                      A tree has 1 root, and <= 2^(height-1) leaves.
 *                      where height = log_2(N+1).
 *                      
 *                      *
 *                  *       *
 *                *  *    *   *
 *               **  **  **   **    
 *                   height = 3, leafs = 4
 *                   height = 4, leafs = 8  <=== 2^(height-1)
 *                      
 *                   nleafs = 2^(height-1) = 2**( log_2(N+1) - 1) = N/2
 *                   
 *                   256 nodes, have height = 8 and nleafs = 128
 *                   
 *            tree memory total = 40 + nleafs*32  + nnodes*48 
 *                              = 40 + (N/2)*32 + (N-N/2-1)*32
 *                              = 40 + 16*N + (N/2 - 1)*32
 *                              = 40 + 32N - 32
 *                              = 32N + 8 bytes
 *                              
 *            
 *            Tree + arguments given to tree:
 *                  8N + 16 bytes + 32N + 8 = 40*N + 24 bytes
 *        
 * @author nichole
 */
public class KDTree {
	
	protected KDTreeNode root = null;
		
	public KDTree(int[] x, int[] y) {
		
		if (x == null) {
			throw new IllegalArgumentException("x cannot be null");
		}
		if (y == null) {
			throw new IllegalArgumentException("y cannot be null");
		}
		if (x.length < 2) {
			throw new IllegalArgumentException("x must be larger than 2");
		}
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y must have same number of points");
		}
	    
		int lastUsableIndex = reduceToUniqueWithMoveUp(x, y);

        this.root = buildTree(0, x, y, 0, lastUsableIndex);
	}

	public KDTreeNode getRoot() {
		return root;
	}
	
	/**
	 * remove unique values by moving up items underneath them.  returns
	 * the last index which should be used in the modified arrays given as arguments.
	 * @param x
	 * @param y
	 * @return
	 */
	static int reduceToUniqueWithMoveUp(int[] x, int[] y) {
		// reduce points to unique by moving them up if they already exist.
            
        Set<PairInt> added = new HashSet<PairInt>();
        
		int count = 0;
		// [0]  0,0      count=0  0,0
		// [1]  1,1            1  1,1
		// [2]  2,2            2  2,2
		// [3]  2,2            *
		// [4]  3,3            3  3,3 
		// [5]  4,4            4  4,4
		// [6]  4,4            *
		// [7]  5,5            5  5,5
		boolean moveUp = false;
		for (int i = 0; i < x.length; i++) {
			if (moveUp) {
				x[count] = x[i];
				y[count] = y[i];
			}
            PairInt p = new PairInt(x[i], y[i]);
			if (!added.contains(p)) {
				added.add(p);
				count++;
			} else {
				moveUp = true;
			}
		}
		return (count-1);
	}
	
	/**
	 * sort by x and return index of median x value for which all points 
	 * below are less than median x and >= for all points above median x index,
	 * in other words 
	 * for indexes : 0 < index have values x < x[index] where x[index] is median x.
	 * for indexes : index <= n have values x >= x[index] where x[index] is median x.
	 * @param x
	 * @param y
	 * @return index that divides the arrays as x<median value and x >= median value
	 */
	int partitionByX(int[] x, int[] y, int startSortRange, int stopSortRangeExclusive) {
				
        QuickSort.sortBy1stArg(x, y, startSortRange, stopSortRangeExclusive);
        
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1);
		int xMedian = x[index];
		while ((index+1) < stopSortRangeExclusive) {
			if (x[index + 1] == xMedian) {
				index++;
			} else {
				break;
			}
		}
		return index;
	}
	
	/**
	 * sort by y and return index of median y value for which all points 
	 * below are less than median y and >= for all points above median y index,
	 * in other words 
	 * for indexes : 0 < index have values y < y[index] where y[index] is median y.
	 * for indexes : index <= n have values y >= y[index] where y[index] is median y.
	 * @param x
	 * @param y
	 * @return index that divides the arrays as y < median value and y >= median value
	 */
	int partitionByY(int[] x, int[] y, int startSortRange, int stopSortRangeExclusive) {
				
        QuickSort.sortBy1stArg(y, x, startSortRange, stopSortRangeExclusive);
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1); // rounds towards zero
		int yMedian = y[index];
		while ((index+1) < stopSortRangeExclusive) {
			if (x[index + 1] == yMedian) {
				index++;
			} else {
				break;
			}
		}
		return index;
	}
        		
	/**
	 * build a tree at from given depth for x, y points using a depth-first algorithm.
	 * Exits from the recursion when npoints = 0, 1 or a leaf.
	 * 
	 * @param depth
	 * @param x
	 * @param y
	 * 
	 * @return
	 */
	protected KDTreeNode buildTree(int depth, int[] x, int[] y, int startSortRange, int stopSortRangeExclusive) {
				
		if (x == null || y == null || x.length == 0 || y.length == 0) {
			return null;
		}
				
		if (stopSortRangeExclusive == startSortRange) {
			// return a leaf of 2 points only
			KDTreeNode leaf = new KDTreeNode();
			leaf.x = x[startSortRange];
			leaf.y = y[startSortRange];			
			return leaf;
		}
		
		int medianIndex = -1;
		int median = 1;
		
	    // if depth of tree is even, partition the points by x, else y
		if (depth % 2 == 0) {
			medianIndex = partitionByX(x, y, startSortRange, stopSortRangeExclusive);
			median = x[medianIndex];
		} else {
			medianIndex = partitionByY(x, y, startSortRange, stopSortRangeExclusive);
			median = y[medianIndex];
		}
		
		depth++;
		
		// left points are  startSortRange through medianX
		KDTreeNode leftChildren = buildTree(depth, x, y, startSortRange, medianIndex);
		
		// right points are medianIndex    through stopSortRangeExclusive
		KDTreeNode rightChildren = buildTree(depth, x, y, medianIndex+1, stopSortRangeExclusive);
		
		KDTreeNode parent = new KDTreeNode();
		
		parent.key = median;
		parent.left = leftChildren;
		parent.right = rightChildren;
		
		leftChildren.parent = parent;
	    rightChildren.parent = parent;
	    
	    if (parent.left != null) {
	    	parent.nChildren += 1 + parent.left.nChildren;
	    }
	    if (parent.right != null) {
	    	parent.nChildren += 1 + parent.right.nChildren;
	    }
	    
	    return parent;
	}
	
	public KDTreeNode findNearestNeighbor(int x, int y) {
		return nearestNeighborSearch(root, x, y);
	}
	
	protected KDTreeNode nearestNeighborSearch(KDTreeNode tree, int leftValue, int rightValue) {
		if ( tree.nChildren == 0 ) {
			return tree;
		}
		
		int medianValue = tree.key;
		
		KDTreeNode subTree;
		
		if (leftValue < medianValue) {
			subTree = tree.left;
		} else {
			subTree = tree.right;
		}
		
		// swap left and right values
		KDTreeNode retVal = nearestNeighborSearch(subTree, rightValue, leftValue);
		
		if (retVal == null) {
			if (tree.left == null) {
				subTree = tree.right;
			} else if (tree.right == null) {
				subTree = tree.left;
			}
			retVal = nearestNeighborSearch(subTree, rightValue, leftValue);
		}
		
		return retVal;
	}
        		
	public void printTree() {
		printTree(root, " ");
	}
	private void printTree(KDTreeNode node, String preString) {
		if (node.left != null) {
		    printTree(node.left, preString + ":LEFT " + node.key);
		}
		if (node.right != null) {
		    printTree(node.right, preString + ":RIGHT " + node.key);
		}
		if (node.left == null && node.right == null) {
			System.out.println(preString + node.key + "(" + node.x + "," + node.y + ")");
			return;
		}
	}
}
