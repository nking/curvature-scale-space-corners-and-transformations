package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.SimpleLinkedListNode;

/**
 * <pre>
 * a Last In First Out structure with
 * 
 *   insert: O(1)
 *   
 *   pop: O(1)
 *   
 *   isEmpty: O(1)
 * </pre>  
 * @author nichole
 *
 */
public class Stack extends SimpleLinkedListNode {

    /**
     *
     * @return
     */
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
    
    /**
     *
     * @return
     */
    public int peek() {
        return this.key;
    }
    
    /**
     *
     * @return
     */
    public boolean isEmpty() {
        
        return (this.key == -1);
    }
    
}
