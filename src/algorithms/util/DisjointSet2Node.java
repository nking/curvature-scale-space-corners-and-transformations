package algorithms.util;

/**
 * a node for the forest implementation of the disjoint set.
 * 
 * @author nichole
 * @param <T> 
 */
public class DisjointSet2Node<T> {

    // the set member data to be held in a disjoint set
    protected T member;
    
    // pointer to the set representative
    protected DisjointSet2Node<T> parent;

    // upper bound of edges between this node and the longest path to it's descendants
    protected int rank = 0;

    public T getMember() {
        return member;
    }

    public void setMember(T member) {
        this.member = member;
    }

    public DisjointSet2Node<T> getParent() {
        return parent;
    }

    public void setParent(DisjointSet2Node<T> theParent) {
        this.parent = theParent;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }
}
