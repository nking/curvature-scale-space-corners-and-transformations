package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.SimpleLinkedListNode;

/**
 * a Last In First Out structure with
 * 
 *   insert: O(1)
 *   
 *   pop: O(1)
 *   
 *   isEmpty: O(1)
 *   
 * @author nichole
 *
 */
public class Stack extends SimpleLinkedListNode {

    /**
     * insert key at top of list
     */
    public SimpleLinkedListNode insert(int insertKey) {
        
        if (insertKey == -1) {
            throw new IllegalArgumentException("key must be larger than -1");
        }
        if (this.key == -1) {
            key = insertKey;
            return this;
        }
        
        // we have to move this.key into a new node which becomes this.next
        SimpleLinkedListNode node = new SimpleLinkedListNode(key);
        
        node.setNext(next);
        
        this.setNext(node);
        
        this.key = insertKey;
      
        return this;
    }
    
    public SimpleLinkedListNode pop() {
        
        if (this.key == -1) {
            
            return null;
        }
        
        SimpleLinkedListNode node = new SimpleLinkedListNode(key);
        
        if (next == null) {
            
            this.key = -1;
            
        } else {
            
            this.key = next.getKey();
            
            this.next = next.getNext();
        }
        
        return node;        
    }
    
    public int peek() {
        return this.key;
    }
    
    public boolean isEmpty() {
        
        return (this.key == -1);
    }
    
}
