package algorithms.disjointSets;

/**
 * a disjoint set implemented with doubly linked lists based upon pseudocode
 * from "Introduction to Algorithms" by Cormen et al.
 *
 * Each set has a single representative which all members point to.   This
 * makes comparing whether two objects are in the same set fast by checking
 * that their representatives are the same.
 *
 * @author nichole
 */
public class DisjointSetHelper {

    /**
     * make a set out of the given node.
     * runtime complexity is O(1).
     *
     * @param x
     * @return
     */
    public <T> DisjointSet<T> makeSet(DisjointSetNode<T> x) {
        x.setRepresentative(x);
        DisjointSet<T> list = new DisjointSet<T>();
        list.setHead(x);
        list.setTail(x);
        list.setNumberOfNodes(1);
        return list;
    }

    /**
     * find the set representative for the given node.
     * runtime complexity is O(1).
     * @param x
     * @return
     */
    public <T> DisjointSetNode<T> findSet(DisjointSetNode<T> x) {
        return x.getRepresentative();
    }

    /**
     * append the shorter list onto the end of the longer's list ("link by size").
     * runtime complexity is  O(N_shorter).
     * @param x
     * @param y
     * @return
     */
    public <T> DisjointSet<T> union(DisjointSet<T> x, DisjointSet<T> y) {

        if (x.equals(y)) {
            return x;
        }
        if (x.getHead().getRepresentative().equals(y.getHead().getRepresentative())) {
            return x;
        }

        DisjointSet<T> longer, shorter;

        if (x.getNumberOfNodes() >= y.getNumberOfNodes()) {
            longer = x;
            shorter = y;
        } else {
            longer = y;
            shorter = x;
        }

        //note that the doubly linked list node characteristic of "previous"
        //is served by the "representative"

        // add next references to longer

        // longer.tail.next might not be pointing to last of next, so walk to end
        if (longer.getTail().getNext() != null) {
            DisjointSetNode<T> tmp = longer.getTail().getNext();
            while (tmp.getNext() != null) {
               tmp = tmp.getNext();
            }
            longer.setTail(tmp);
        }

        longer.getTail().setNext(shorter.getHead());

        DisjointSetNode<T> latest = shorter.getHead();
        while (latest != null) {
            latest.setRepresentative(longer.getHead());
            latest = latest.getNext();
        }
        longer.setTail(shorter.getTail());

        longer.setNumberOfNodes(longer.getNumberOfNodes() + shorter.getNumberOfNodes());

        return longer;
    }

    public static <T> String print(DisjointSet<T> x) {

        DisjointSetNode<T> current = x.getHead();

        int nIter = 0;

        StringBuilder sb = new StringBuilder();
        while (current != null && (nIter < x.getNumberOfNodes())) {
            if (sb.length() > 0) {
                sb.append(", ");
            }
            sb.append("[").append(current.getMember().toString()).append("] ")
            .append("repr=").append(current.getRepresentative().getMember().toString())
            ;
            current = current.getNext();
            nIter++;
        }
        sb.append(", tail=").append(x.getTail().getMember().toString());

        return sb.toString();
    }

}
