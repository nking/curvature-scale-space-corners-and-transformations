package algorithms.disjointSets;

/**
 * a node for the forest implementation of the disjoint set.
 * 
 * @author nichole
 * @param <T> 
 */
public class DisjointSet2Node<T> {

    /**
     * the set member data to be held in a disjoint set
     */
    protected T member = null;
    
    /**
     * pointer to the set representative
     */
    protected DisjointSet2Node<T> parent = null;

    /**
    upper bound of edges between this node and the longest path to it's descendants
    */
    protected int rank = 0;
    
    protected Object data = null;

    public DisjointSet2Node(T member) {
        this.member = member;
    }
    
    public DisjointSet2Node() {
    }
   
    /**
     * get the member data
     * @return 
     */
    public T getMember() {
        return member;
    }

    public void setMember(T member) {
        this.member = member;
    }
    
    public void setDeta(Object data) {
        this.data = data;
    }
    
    public Object getObject() {
        return data;
    }

    /**
     * get the set representative
     * @return 
     */
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

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("[member=");
        if (member != null) {
            sb.append(member.toString());
        }
        sb.append("; ");
        sb.append("rank=").append(Integer.toString(rank)).append("; ");
        
        sb.append("data=");
        if (data != null) {
            sb.append(data.toString());
        }
        sb.append("; ");
        
        sb.append("parent=");
        if (parent != null) {
            if (parent.equals(this)) {
                sb.append("self");
            } else {
                sb.append(parent.toString());
            }
        }
        sb.append(";] ");
        
        return sb.toString();
    }
  
}