package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 <pre>
  A holder for two-point identities, where the identities are the indexes
  of the indexer internal arrays.  N is the size of the dataset, that is
  the indexer.nXY.
  
  Runtime complexity:
     inserts are     O(lg₂(N)) at best and O(N) at worse.
     comparisons are O(lg₂(N))
 
  Space complexity:
     O(N)
 
  Note, could implement a balanced tree to make inserts at worse O(lg₂(N)).
 </pre>
 
 * @author nichole
 */
class TwoPointBinarySearchTree implements ITwoPointIdentity {

    protected int n = 0;

    protected Node root = null;

    /**
     * constructor
     */
    TwoPointBinarySearchTree() {
    }

    @Override
    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int arrayRefBits = 32;

        int overheadBytes = 16;

        // each Node:  overhead + int + int + 3 references
        int oneNodeInBits = 2*nbits * 3*arrayRefBits;
        int oneNodeInBytes = (oneNodeInBits/8) + overheadBytes;
        oneNodeInBytes += (oneNodeInBytes % 8);

        int nNodesInBytes = n * oneNodeInBytes;

        //                              root ref        nodes              1 ints
        long sumBytes = overheadBytes + (nbits/8)     + nNodesInBytes  + (nbits/8);

        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    /**
     * check if combination is already stored, if not add it and return true, else
     * return false
     *
     * @param index0
     * @param index1
     * @return
     */
    @Override
    public boolean storeIfDoesNotContain(int index0, int index1) {

        // order the indexes to avoid double counting.
        int i0, i1;
        if (index0 < index1) {
            i0 = index0;
            i1 = index1;
        } else {
            i0 = index1;
            i1 = index0;
        }

        Node node = search(i0, i1);

        if (node != null) {
            return false;
        }

        insert(i0, i1);

        return true;
    }

    public Node search(int i0, int i1) {
        Node x = root;
        while ((x != null) && (x.compare(i0, i1) != 0)) {
            if (x.compare(i0, i1) > 0) {
                x = x.getLeft();
            } else {
                x = x.getRight();
            }
        }
        return x;
    }

    protected void insert(int i0, int i1) {
        Node y = null;
        Node x = root;
        while (x != null) {
            y = x;
            if (x.compare(i0, i1) > 0) {
                x = x.getLeft();
            } else {
                x = x.getRight();
            }
        }

        Node z = new Node(i0, i1);
        z.parent = y;
        if (y == null) {
            root = z;
        } else if (y.compare(z) > 0) {
            y.setLeft(z);
        } else {
            y.setRight(z);
        }
        n++;
    }

    protected class Node {
        protected final int a0;
        protected final int a1;
        protected Node left;
        protected Node right;
        protected Node parent;
        public Node(int i0, int i1) {
            this.a0 = i0;
            this.a1 = i1;
        }
        public int compare(Node other) {
            return compare(other.a0, other.a1);
        }
        public int compare(int other0, int other1) {
            // for now, using n0.a0 < n1.a0 then n0.a1 < n1.a1
            if (a0 < other0) {
                return -1;
            } else if (a0 > other0) {
                return 1;
            }
            if (a1 < other1) {
                return -1;
            } else if (a1 > other1) {
                return 1;
            }
            return 0;
        }
        public Node getLeft() {
            return left;
        }
        public Node getRight() {
            return right;
        }
        public Node getParent() {
            return parent;
        }
        public void setLeft(Node lft) {
            left = lft;
        }
        public void setRight(Node rght) {
            right = rght;
        }
        public void setParent(Node prnt) {
            parent = prnt;
        }
    }

}
