package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 <pre>
  A holder for two-point identities, where the identities are the indexes
  of the indexer internal arrays.  N is the size of the dataset, that is
  the indexer.nXY.
  
  Runtime complexity:
     inserts is  O(1)
     search  is O(1)
 
  Space complexity:
     O(N)
 
 </pre>
 
 * @author nichole
 */
class TwoPointHash implements ITwoPointIdentity {

    protected int n = 0;
    
    protected int capacity = 1000;
    
    protected final Node[] nodes;
    
    public TwoPointHash(int indexerNXY) {
        
        //TODO: refine a capacity based upon indexerNXY
        
        if (indexerNXY > 100000) {
            capacity = 10000;
        } else if (indexerNXY > 1000) {
            capacity = indexerNXY;
        }
        
        nodes = new Node[capacity];
    }    

    @Override
    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int arrayRefBits = 32;

        int overheadBytes = 16;

        // each Node:  overhead + int + int + 3 references
        int oneNodeInBits = 2*nbits * arrayRefBits;
        int oneNodeInBytes = (oneNodeInBits/8) + overheadBytes;
        oneNodeInBytes += (oneNodeInBytes % 8);

        int nNodesInBytes = n * oneNodeInBytes;

        //                              nodes ref        nodes            2 ints
        long sumBytes = overheadBytes + arrayRefBits + nNodesInBytes  + (2*nbits/8);

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

    /**
     * search is O(1). (unless there are many collisions)
     * @param i0
     * @param i1
     * @return 
     */
    public Node search(int index0, int index1) {
        
        int i0, i1;
        if (index0 < index1) {
            i0 = index0;
            i1 = index1;
        } else {
            i0 = index1;
            i1 = index0;
        }
        
        int h = hash(i0, i1);
        
        Node p = nodes[h];
        
        while ((p != null) && (p.compare(i0, i1) != 0)) {
            p = p.getNext();
        }
        
        return p;
    }

    /**
     * insert is O(1). (unless there are many collisions)
     * @param i0
     * @param i1 
     */
    protected void insert(int i0, int i1) {
        
        int h = hash(i0, i1);
        
        Node p = nodes[h];
        
        n++;
        
        if (p == null) {
            nodes[h] = new Node(i0, i1);
            return;
        }
        
        while ((p != null) && (p.getNext() != null)) {
            p = p.getNext();
        }
        
        p.next = new Node(i0, i1);        
    }

    protected class Node {
        protected final int a0;
        protected final int a1;
        protected Node next;
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
        public Node getNext() {
            return next;
        }
        public void setNext(Node node) {
            next = node;
        }
    }
    
    protected int hash(int i0, int i1) {
        
        int fnv = fnvHashCode(i0, i1);
        
        int h = fnv % (capacity - 1);
        if (fnv < 0) {
            h += (capacity - 1);
        }
        
        return h;
    }

    protected static int fnv321aInit = 0x811c9dc5;
    protected static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(int i0, int i1) {

        /*
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
         */

        int hash = 0;

        int sum = fnv321aInit;

        // xor the bottom with the current octet.
        sum ^= i0;

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        sum ^= i1;
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

}
