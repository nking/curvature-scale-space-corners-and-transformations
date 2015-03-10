package algorithms.util;

/**
 * <pre>
 * a Last In First Out structure with
 * 
 *   insert: O(1)
 *   
 *   pop: O(1)
 *   
 *   isEmpty: O(1)
 * 
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/util/Stack.java
 * under MIT License (MIT), Nichole King 2013
 * 
 * </pre>  
 * @author nichole
 *
 */
public class StackInt extends SimpleLinkedListNode {

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
