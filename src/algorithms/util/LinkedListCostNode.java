package algorithms.util;

/**
 * extends SimpleLinkedListNode to hold an integer for the cost (which is
 * by default Integer.MAX_VALUE for use with min heaps and priority queues).
 * (NOTE, could edit the code to make defaultCost modifiable.)
 * 
 * @author nichole
 */
public class LinkedListCostNode extends SimpleLinkedListNode {
    
    public final static int DEFAULT_COST = Integer.MAX_VALUE;
    
    protected int cost = DEFAULT_COST;
    
    public LinkedListCostNode() {
        super();
    }
    
    public LinkedListCostNode(int insertKey, int cost) {
        super(insertKey);
        this.cost = cost;
    }
    
    public LinkedListCostNode(int insertKey) {
        super(insertKey);
    }
    
    public int getCost() {
        return cost;
    }
    
    @Override
    public SimpleLinkedListNode insert(int insertKey) {
        return insert(insertKey, DEFAULT_COST);
    }
    
    public LinkedListCostNode insert(int insertKey, int insertCost) {
         if (insertKey == -1) {
            throw new IllegalArgumentException(
            "insertKey must be larger than -1");
        }
        n++;
        if (this.key == -1) {
            key = insertKey;
            cost = insertCost;
            return this;
        }
        
        LinkedListCostNode node = new LinkedListCostNode(key, cost);
        
        key = insertKey;
        cost = insertCost;

        if (next == null) {
            next = node;
            return this;
        }
        
        node.next = next;
        
        next = node;
        
        return node;
    }
    
    @Override
    public SimpleLinkedListNode insertIfDoesNotAlreadyExist(int insertKey) {
        return insertIfDoesNotAlreadyExist(insertKey, DEFAULT_COST);
    }
    
    public LinkedListCostNode insertIfDoesNotAlreadyExist(int insertKey, int insertCost) {
        
        if (insertKey == -1) {
            throw new IllegalArgumentException(
            "insertKey must be larger than -1");
        }
        if (insertKey == this.key) {
            return null;
        }
        if (this.key == -1) {
            key = insertKey;
            cost = insertCost;
            return this;
        }
        
        SimpleLinkedListNode node = search(insertKey);
        
        if (node != null) {
            return null;
        }
        
        return insert(insertKey, insertCost);
    }
    
    public void delete(LinkedListCostNode node) {

        if (key == -1) {
            return;
        }
        
        // if its the first node, we have to transfer data
        if (this.equals(node)) {
            if (this.next == null) {
                this.key = -1;
                this.cost = DEFAULT_COST;
            } else {
                this.key = next.key;
                this.cost = ((LinkedListCostNode)next).cost;
                this.next = next.next;
            }
            n--;
            return;
        }

        // start w/ 2nd node because we've already searched the first
        LinkedListCostNode last = this;
        
        LinkedListCostNode current = last;

        while (current.next != null) {
            current = (LinkedListCostNode)current.next;
            if (current.equals(node)) {
                last.next = current.next; 
                n--;
                break;
            }
            last = current;            
        }
    }
    
    /**
     * delete the first node found with key == deleteKey.
     * 
     * @param deleteKey 
     */
    @Override
    public void delete(int deleteKey) {

        if (deleteKey == -1) {
            return;
        }
        
        // if its the first node, we have to transfer data
        if (this.key == deleteKey) {
            if (this.next == null) {
                this.key = -1;
                this.cost = DEFAULT_COST;
            } else {
                this.key = next.key;
                this.cost = ((LinkedListCostNode)next).getCost();
                this.next = next.next;
            }
            n--;
            return;
        }
        
        // start w/ 2nd node because we've already searched the first
        LinkedListCostNode last = this;
        
        LinkedListCostNode current = last;

        while (current.next != null) {
            current = (LinkedListCostNode)next;
            if (current.key == deleteKey) {
                last.next = current.next;
                n--;
                break;
            }
            last = current;
        }
    }
    
    public static long approximateMemoryUsed() {
        
        String arch = System.getProperty("sun.arch.data.model");
    
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
    
        int nbits = (is32Bit) ? 32 : 64;
    
        int overheadBytes = 16;
    
        int intBytes = (is32Bit) ? 4 : 8;
        // 4 ints:
        intBytes *= 4;
        
        int refBytes = nbits/8;

        long sumBytes = intBytes + refBytes;
       
        sumBytes += overheadBytes;
        
        long padding = (sumBytes % 8);
        
        sumBytes += padding;
        
        return sumBytes;
    }
    
    /**
     * only the key is used for this equals identity
     * @param arg0
     * @return 
     */
    @Override
    public boolean equals(Object arg0) {
        if (!(arg0 instanceof LinkedListCostNode)) {
            return false;
        }
        LinkedListCostNode other = (LinkedListCostNode)arg0;
        
        return (other.key == this.key);
    }
    
    @Override
    public int hashCode() {
        // even if same keys, want different hashcodes
        return super.hashCode(); 
    }
}
