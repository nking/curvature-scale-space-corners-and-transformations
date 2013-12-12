package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;

/**
 * a node holding only a integer key and the next reference
 *
 * @author nichole
 */
public class SimpleLinkedListNode {

    public int key = -1;

    public SimpleLinkedListNode next = null;

    public SimpleLinkedListNode() {}
    public SimpleLinkedListNode(int insertKey) {
        this.key = insertKey;
    }
    
    public int getKey() {
        return key;
    }

    public SimpleLinkedListNode insert(int insertKey) {
        if (insertKey == -1) {
            throw new IllegalArgumentException("key must be larger than -1");
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

    /*
    SimpleLinkedListNode latest = this;
    while (latest != null) {
        if (latest.key == key) {
            return latest;
        }
        latest = latest.next;
    }
    */

    public SimpleLinkedListNode insertIfDoesNotAlreadyExist(int insertKey) {
        if (insertKey == -1) {
            throw new IllegalArgumentException("key must be larger than -1");
        }
        if (insertKey == this.key) {
            return null;
        }
        if (this.key == -1) {
            key = insertKey;
            return this;
        }

        if (next == null) {
            next = new SimpleLinkedListNode(insertKey);
            return next;
        } else {
        	SimpleLinkedListNode last = next;
            if (last.key == insertKey) {
                return null;
            }
            while (last.next != null) {
                if (last.key == insertKey) {
                    return null;
                }
                last = last.next;
            }
            last.next = new SimpleLinkedListNode(insertKey);
            return last.next;
        }
    }

    public void delete(SimpleLinkedListNode node) {

        if (key == -1) {
            return;
        }

        SimpleLinkedListNode last = null;
        SimpleLinkedListNode current = this;

        while (current != null) {
            if (current.equals(node)) {
                if (last == null) {
                    if (next != null) {
                        // this is the first node in the list.  reassign field values to next
                        this.key = next.key;
                        this.next = next.next;
                    } else {
                        // this is the first and only node
                        this.key = -1;
                    }
                } else {
                    last.next = current.next;
                }
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

        SimpleLinkedListNode last = null;
        SimpleLinkedListNode current = this;

        while (current != null && current.key != -1) {
            if (current.key == deleteKey) {
                if (last == null) {
                    if (next != null) {
                        // this is the first node in the list.  reassign field values to next
                        this.key = next.key;
                        this.next = next.next;
                    } else {
                        // this is the first and only node
                        this.key = -1;
                    }
                } else {
                    // else, this node is not the first in the list, we can skip over it to delete it
                    last.next = current.next;
                }
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
    
}
