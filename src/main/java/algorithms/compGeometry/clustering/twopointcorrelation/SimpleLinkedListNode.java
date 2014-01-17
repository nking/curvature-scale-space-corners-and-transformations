package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;

/**
 * a node holding only a integer key and the next reference.  the key must be larger than -1.
 *
 * @author nichole
 */
public class SimpleLinkedListNode {

    protected int key = -1;

    protected SimpleLinkedListNode next = null;

    public SimpleLinkedListNode() {}
    
    public SimpleLinkedListNode(int insertKey) {
        this.key = insertKey;
    }
    
    public int getKey() {
        return key;
    }
    
    public SimpleLinkedListNode getNext() {
        return next;
    }

    /**
     * set next to nextNode.  note that if this.next is not null, it is lost.
     * @param nextNode
     */
    public void setNext(SimpleLinkedListNode nextNode) {
        this.next = nextNode;
    }
    
    public SimpleLinkedListNode insert(int insertKey) {
        if (insertKey == -1) {
            throw new IllegalArgumentException("insertKey must be larger than -1");
        }
        if (this.key == -1) {
            key = insertKey;
            return this;
        }
        SimpleLinkedListNode node = new SimpleLinkedListNode(insertKey);

        SimpleLinkedListNode last = this;
        while (last.next != null) {
            last = last.next;
        }
        last.next = node;

        return node;
    }

    public SimpleLinkedListNode insertIfDoesNotAlreadyExist(int insertKey) {
        if (insertKey == -1) {
            throw new IllegalArgumentException("insertKey must be larger than -1");
        }
        if (insertKey == this.key) {
            return null;
        }
        if (this.key == -1) {
            key = insertKey;
            return this;
        }
        
        SimpleLinkedListNode last = this;
        SimpleLinkedListNode current = next;
        while (current != null) {
            if (current.key == insertKey) {
                return null;
            }
            last = current;
            current = current.next;
        }
        
        last.next = new SimpleLinkedListNode(insertKey);
        return last.next;
        
    }

    public void delete(SimpleLinkedListNode node) {

        if (key == -1) {
            return;
        }
        
        // if its the first node, we have to transfer data
        if (this.equals(node)) {
            if (this.next == null) {
                this.key = -1;
            } else {
                this.key = next.key;
                this.next = next.next;
            }
            return;
        }

        // start w/ 2nd node because we've already searched the first
        SimpleLinkedListNode last = this;
        SimpleLinkedListNode current = this.next;

        while (current != null) {
            if (current.equals(node)) {
                last.next = current.next;               
                break;
            }
            last = current;
            current = current.next;
        }
    }

    public void delete(int deleteKey) {

        if (deleteKey == -1) {
            return;
        }
        
        // if its the first node, we have to transfer data
        if (this.key == deleteKey) {
            if (this.next == null) {
                this.key = -1;
            } else {
                this.key = next.key;
                this.next = next.next;
            }
            return;
        }
        
        // start w/ 2nd node because we've already searched the first
        SimpleLinkedListNode last = this;
        SimpleLinkedListNode current = this.next;

        while (current != null) {
            if (current.key == deleteKey) {
                last.next = current.next;               
                break;
            }
            last = current;
            current = current.next;
        }
    }

    public SimpleLinkedListNode search(int searchKey) {

        SimpleLinkedListNode latest = this;

        while (latest != null) {
            if (latest.key == searchKey) {
                return latest;
            }
            latest = latest.next;
        }
        return null;
    }

    public boolean contains(int searchKey) {
        SimpleLinkedListNode node = search(searchKey);
        return (node != null);
    }

    public int[] getKeys() {
        if (key == -1) {
            return new int[0];
        }
        int n = 10;
        int[] nodeKeys = new int[n];
        int count = 0;

        SimpleLinkedListNode latest = this;
        while (latest != null && latest.key != -1) {
            if ((count + 1) > n) {
                n = 2*n;
                nodeKeys = Arrays.copyOf(nodeKeys, n);
            }
            nodeKeys[count] = latest.key;
            count++;
            latest = latest.next;
        }
        return Arrays.copyOf(nodeKeys, count);
    }

    public int getNumberOfKeys() {
        if (key == -1) {
            return 0;
        }
        int count = 0;

        SimpleLinkedListNode latest = this;
        while (latest != null && latest.key != -1) {
            count++;
            latest = latest.next;
        }
        return count;
    }
    
    public static long approximateMemoryUsed() {
            
        String arch = System.getProperty("sun.arch.data.model");
    
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
    
        int nbits = (is32Bit) ? 32 : 64;
    
        int overheadBytes = 16;
    
        int intBytes = (is32Bit) ? 4 : 8;
        
        int refBytes = nbits/8;

        long sumBytes = intBytes + refBytes;
       
        sumBytes += overheadBytes;
        
        long padding = (sumBytes % 8);
        
        sumBytes += padding;
        
        return sumBytes;
    }
    
    @Override
    public boolean equals(Object arg0) {
        if (!(arg0 instanceof SimpleLinkedListNode)) {
            return false;
        }
        SimpleLinkedListNode other = (SimpleLinkedListNode)arg0;
        return (other.key == this.key);
    }
}
