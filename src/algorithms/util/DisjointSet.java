package algorithms.util;

public class DisjointSet {

  //1st obj in each linked list is the set representative
    DisjointSetNode head;
    DisjointSetNode tail;
    int numberOfNodes;
    
    public DisjointSetNode getHead() {
        return head;
    }
    public void setHead(DisjointSetNode head) {
        this.head = head;
    }
    public DisjointSetNode getTail() {
        return tail;
    }
    public void setTail(DisjointSetNode tail) {
        this.tail = tail;
    }
    public int getNumberOfNodes() {
        return numberOfNodes;
    }
    public void setNumberOfNodes(int numberOfNodes) {
        this.numberOfNodes = numberOfNodes;
    }
    
}
