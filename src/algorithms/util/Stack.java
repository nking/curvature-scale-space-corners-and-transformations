package algorithms.util;

/**
 * data structure in which last inserted is first out (LIFO).
 * it's primary methods are push and pop, both of which have
 * runtime complexity of O(1).
 * 
 * The structure has a space complexity of O(N).
 * 
 * @author nichole
 *
 */
public class Stack<T> {

    private Node<T> list = null;
    
    private int n = 0;
    
    public void push(T obj) {
        if (obj == null) {
            return;
        }
        
        Node<T> nd = new Node<T>();
        nd.data = obj;
        
        if (list == null) {
            list = nd;
        } else {
            nd.next = list;
            this.list = nd;
        }
        n++;
    }
    
    public T pop() {
        if (list == null) {
            return null;
        }
        Node<T> top = list;
        list = top.next;
        n--;
        return top.data;
    }
    
    public T peek() {
        if (list == null) {
            return null;
        }
        return list.data;
    }
    
    public T peekPopNext() {
        if (list == null || list.next == null) {
            return null;
        }
        return list.next.data;
    }
    
    public boolean isEmpty() {
        return (list == null);
    }
    
    public int size() {
        return n;
    }
    
    protected static class Node<T> {
        Node<T> next = null;
        T data = null;
    }
}
