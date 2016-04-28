package algorithms.search;

public class KDTreeNode {

	KDTreeNode right = null;
	KDTreeNode left = null;
	KDTreeNode parent = null;
	
	int x = -1;
	int y = -1;
	int key = -1; // median value
	int nChildren = 0;
}
