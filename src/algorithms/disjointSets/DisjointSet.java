package algorithms.disjointSets;

/**
 * based upon pseudocode from "Introduction to Algorithms" by Cormen et al.
 * 
 * @author nichole
 * @param <T> 
 */
public class DisjointSet<T> {

    //1st obj in each linked list is the set representative
    protected DisjointSetNode<T> head;
    protected DisjointSetNode<T> tail;
    protected int numberOfNodes;
    
    public DisjointSetNode<T> getHead() {
        return head;
    }
    public void setHead(DisjointSetNode<T> head) {
        this.head = head;
    }
    public DisjointSetNode<T> getTail() {
        return tail;
    }
    public void setTail(DisjointSetNode<T> tail) {
        this.tail = tail;
    }
    public int getNumberOfNodes() {
        return numberOfNodes;
    }
    public void setNumberOfNodes(int numberOfNodes) {
        this.numberOfNodes = numberOfNodes;
    }
    
}
