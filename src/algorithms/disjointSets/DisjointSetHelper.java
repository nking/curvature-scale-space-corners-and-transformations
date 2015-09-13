package algorithms.disjointSets;

/**
 * a disjoint set implemented with linked lists.
 * each set is a linked list.
 *
 * based upon pseudocode from "Introduction to Algorithms" by Cormen et al.
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
     * append the shorter list onto the end of the longer's list.
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

        StringBuilder sb = new StringBuilder();
        while (current != null) {
            if (sb.length() > 0) {
                sb.append(", ");
            }
            sb.append(current.getMember().toString());
            current = current.getNext();
        }

        return sb.toString();
    }

}
