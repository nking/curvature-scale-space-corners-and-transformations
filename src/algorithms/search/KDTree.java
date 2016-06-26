package algorithms.search;

import algorithms.QuickSort;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.Arrays;

/**
 * k-dimension tree is a binary tree used to store coordinates for quick nearest
 * neighbor or range searches.
 *    -- values are always stored in leaves
 *    -- the meaning of an internal node depends upon the depth of the
 *       node within the tree.
 *       
 * Note that learning the true medium while 
 * constructing the tree takes more
 * time, but leads to a better balanced tree.
 * 
 * adapted from pseudocode from
 * http://ldots.org/kdtree which licenses the content as:
 * http://creativecommons.org/licenses/by-sa/2.0/
 * 
 * useful reading regarding best distances in nearest neighbor search:
 *    http://web.stanford.edu/class/cs106l/handouts/assignment-3-kdtree.pdf
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

        if (lastUsableIndex < (x.length - 1)) {
            x = Arrays.copyOf(x, lastUsableIndex + 1);
            y = Arrays.copyOf(y, lastUsableIndex + 1);
        }

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
			if (y[index + 1] == yMedian) {
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
    
    private KDTreeNode bestNode = null;
    private double bestDist = Double.MAX_VALUE;
        
	public KDTreeNode findNearestNeighbor(int x, int y) {
        bestNode = null;
        bestDist = Double.MAX_VALUE;
        
        KDTreeNode node = nearestNeighborSearch(root, x, y, 0);
        return node;
    }
	
	protected KDTreeNode nearestNeighborSearch(KDTreeNode tree, int leftValue, 
        int rightValue, int depth) {

        if (tree.nChildren == 0 ) {
			return tree;
		}
		
		int medianValue = tree.getKey();
		
		float diffMedValSq;
		
		KDTreeNode subTree1, subTree2;

        if ((depth & 1) == 0) {
            diffMedValSq = medianValue - leftValue;
            if (leftValue <= medianValue) {
                subTree1 = tree.left;
                subTree2 = tree.right;
            } else {
                subTree1 = tree.right;
                subTree2 = tree.left;
            }
        } else {
            diffMedValSq = medianValue - rightValue;
            if (rightValue <= medianValue) {
                subTree1 = tree.left;
                subTree2 = tree.right;
            } else {
                subTree1 = tree.right;
                subTree2 = tree.left;
            }
        }
        diffMedValSq *= diffMedValSq;
	
        KDTreeNode retVal1 = nearestNeighborSearch(
            subTree1, leftValue, rightValue, depth + 1);
		
        double dist1 = Double.MAX_VALUE;
        if (retVal1 != null) {
            dist1 = distanceSq(retVal1, leftValue, rightValue);
            // TODO: consider a tolerance
            if (dist1 == 0) {
                // this is the point
                bestDist = dist1;
                bestNode = retVal1;
                return retVal1;
            }
        }
        
        System.out.println("dist1=" + dist1 
            + "  med-key="+ diffMedValSq + 
            "  bestdist=" + bestDist);
        
        //TODO: this may need to be revised for a radius.
        //   basically, if (leftValue, rightValue) is closer to
        //      the median than it is to retVal1,
        //      search the other tree too.
        
		if ((2*diffMedValSq) < dist1) {
            
			KDTreeNode retVal2 = nearestNeighborSearch(
                subTree2, leftValue, rightValue, depth + 1);
            
            double dist2 = Double.MAX_VALUE;
            if (retVal2 != null) {
                dist2 = distanceSq(retVal2, leftValue, rightValue);
                // TODO: consider a tolerance
                if (dist2 == 0) {
                    // this is the point
                    bestDist = dist2;
                    bestNode = retVal2;
                    return bestNode;
                }
                if (dist2 < dist1) {
                    dist1 = dist2;
                    retVal1 = retVal2;
                }
            }
        }
        
        if (dist1 < bestDist && retVal1 != null) {
            bestDist = dist1;
            bestNode = retVal1;
        }
        
		return bestNode;
    }
    
    private double distanceSq(KDTreeNode tree, 
        float leftValue, float rightValue) {

        float diffX = tree.getX() - leftValue;
        float diffY = tree.getY() - rightValue;
        
        return (diffX * diffX) + (diffY * diffY);
    }
    
	public void printTree() {
		printTree(root, " ");
	}
	private void printTree(KDTreeNode node, String preString) {
		if (node.left != null) {
		    printTree(node.left, preString + ":LEFT " + node.getKey());
		}
		if (node.right != null) {
		    printTree(node.right, preString + ":RIGHT " + node.getKey());
		}
		if (node.left == null && node.right == null) {
			System.out.println(preString + node.getKey() + "(" + node.getX() + "," + node.getY() + ")");
			return;
		}
	}
}
