package algorithms.search;

import algorithms.MultiArrayMergeSort;
import algorithms.util.PairFloat;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashSet;
import java.util.Set;

/**
 * k-dimension tree is a binary tree used to store coordinates for quick nearest
 * neighbor or range searches.
 * This one is a 2D tree only.
 * 
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
public class KDTreeFloat {
	
	protected KDTreeNodeFloat root = null;
		
	public KDTreeFloat(float[] xPoints, float[] yPoints, boolean alreadySorted) {
		
		if (xPoints == null) {
			throw new IllegalArgumentException("xPoints cannot be null");
		}
		if (yPoints == null) {
			throw new IllegalArgumentException("yPoints cannot be null");
		}
		if (xPoints.length < 2) {
			throw new IllegalArgumentException("xPoints must be larger than 2");
		}
		if (xPoints.length != yPoints.length) {
			throw new IllegalArgumentException("xPoints and yPoints must have same number of points");
		}
        
        int lastUsableIndex = reduceToUniqueWithMoveUp(xPoints, yPoints);
	    
        if (lastUsableIndex < (xPoints.length - 1)) {
            xPoints = Arrays.copyOf(xPoints, lastUsableIndex + 1);
            yPoints = Arrays.copyOf(yPoints, lastUsableIndex + 1);
        }
        
        //TODO: to better handle space for large number of points,
        //  could change out the merge sorts for quick sorts
        
        if (!alreadySorted) {
            MultiArrayMergeSort.sortBy1stArgThen2nd(xPoints, yPoints);
        }
                
        this.root = buildTree(0, xPoints, yPoints, 0, lastUsableIndex);
	}
    
	/**
	 * remove unique values by moving up items underneath them.  returns
	 * the last index which should be used in the modified arrays given as arguments.
	 * @param x
	 * @param y
	 * @return
	 */
	static int reduceToUniqueWithMoveUp(float[] x, float[] y) {
		// reduce points to unique by moving them up if they already exist.
            
        Set<PairFloat> added = new HashSet<PairFloat>();
        
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
            PairFloat p = new PairFloat(x[i], y[i]);
			if (!added.contains(p)) {
				added.add(p);
				count++;
			} else {
				moveUp = true;
			}
		}
        
		return added.size() - 1;
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
	protected KDTreeNodeFloat buildTree(int depth, float[] x, float[] y, 
        int startSortRange, int stopSortRangeExclusive) {

		if (x == null || y == null || x.length == 0 || y.length == 0) {
			return null;
		}
				
		if (stopSortRangeExclusive == startSortRange) {
			KDTreeNodeFloat leaf = new KDTreeNodeFloat();
            leaf.depth = depth;
            leaf.key = Integer.MIN_VALUE;
			leaf.x = x[startSortRange];
			leaf.y = y[startSortRange];			
			return leaf;
		}
		
		int medianIndex = -1;
		float median = 1;
		
	    // if depth of tree is even, partition the points by x, else y
		if ((depth & 1) == 0) {
			medianIndex = partitionByX(
                x, y, startSortRange, stopSortRangeExclusive);
			median = x[medianIndex];
		} else {
			medianIndex = partitionByY(
                x, y, startSortRange, stopSortRangeExclusive);
			median = y[medianIndex];
		}
                
		depth++;
       
		// left points are  startSortRange through medianX
		KDTreeNodeFloat leftChildren = buildTree(depth, x, y, 
            startSortRange, medianIndex);
		
		// right points are medianIndex    through stopSortRangeExclusive
		KDTreeNodeFloat rightChildren = buildTree(depth, x, y, 
            medianIndex+1, stopSortRangeExclusive);
		
		KDTreeNodeFloat parent = new KDTreeNodeFloat();
		
        parent.depth = depth - 1;
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
	int partitionByX(float[] x, float[] y, int startSortRange, int stopSortRangeExclusive) {
				
        MultiArrayMergeSort.sortBy1stArg(
            x, y, startSortRange, stopSortRangeExclusive);
        
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1);
		float xMedian = x[index];
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
	int partitionByY(float[] x, float[] y, int startSortRange, int stopSortRangeExclusive) {
				
        MultiArrayMergeSort.sortBy1stArg(y, x, startSortRange, stopSortRangeExclusive);
		int n = (stopSortRangeExclusive - startSortRange);
		int index = startSortRange + (n >> 1); // rounds towards zero
		float yMedian = y[index];
		while ((index+1) < stopSortRangeExclusive) {
			if (y[index + 1] == yMedian) {
				index++;
			} else {
				break;
			}
		}
		return index;
	}	
    
    private KDTreeNodeFloat bestNode = null;
    private double bestDist = Double.MAX_VALUE;
        
	public KDTreeNodeFloat findNearestNeighbor(float x, float y) {
        bestNode = null;
        bestDist = Double.MAX_VALUE;
        return nearestNeighborSearch(root, x, y, 0);
	}
	
	protected KDTreeNodeFloat nearestNeighborSearch(
        KDTreeNodeFloat tree, float leftValue, float rightValue,
        int depth) {
        
		if (tree.nChildren == 0 ) {
			return tree;
		}
        
		float medianValue = tree.getKey();
        
        float diffMedValSq;
		
		KDTreeNodeFloat subTree1, subTree2;

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
	 
		KDTreeNodeFloat retVal1 = nearestNeighborSearch(
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
        
        //System.out.println("dist1=" + dist1 
        //    + "  med-key="+ diffMedValSq + 
        //    "  bestdist=" + bestDist);
        
        //TODO: this may need to be revised for a radius.
        //   basically, if (leftValue, rightValue) is closer to
        //      the median than it is to retVal1,
        //      search subtree2 too.
        
		if ((2*diffMedValSq) < dist1) {
			KDTreeNodeFloat retVal2 = nearestNeighborSearch(
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
        
        if (dist1 < bestDist) {
            bestDist = dist1;
            bestNode = retVal1;
        }
        
		return bestNode;
	}
        		
	public void printTree() {
		printTree(root, " ");
	}
	private void printTree(KDTreeNodeFloat node, String preString) {
		if (node == null) {
            return;
        }
        
        Deque<KDTreeNodeFloat> q0 = new ArrayDeque<KDTreeNodeFloat>();
        Deque<KDTreeNodeFloat> q1 = new ArrayDeque<KDTreeNodeFloat>();
        q0.offer(node);

        int count = 0;
        boolean skip = true;
        while(!q0.isEmpty()) {
            while(!q0.isEmpty()) {
                node = q0.poll();
                System.out.println("level=" + node.depth + " (med=" + node.getKey() +
                    " x=" + node.getX() + " y=" + node.getY() + ")");
                if (node.left != null) {
                    q1.offer(node.left);
                }
                if (node.right != null) {
                    q1.offer(node.right);
                }
            }
            if (!skip) {
                count++;
            } else {
                skip = false;
            }
            q0.addAll(q1);
            q1.clear();
        }
	}
    
	public KDTreeNodeFloat getRoot() {
		return root;
	}

    private double distanceSq(KDTreeNodeFloat tree, 
        float leftValue, float rightValue) {

        float diffX = tree.getX() - leftValue;
        float diffY = tree.getY() - rightValue;
        
        return (diffX * diffX) + (diffY * diffY);
    }
	
}
