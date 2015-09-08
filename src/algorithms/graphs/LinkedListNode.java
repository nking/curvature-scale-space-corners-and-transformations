package algorithms.graphs;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class LinkedListNode {

    protected int key = -1;

    protected int[] incoming = null;
    protected int nIncoming = 0;
    protected int[] outgoing = null;
    protected int nOutgoing = 0;

    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeNumber
     */
    public LinkedListNode(int nodeNumber) {
        this.key = nodeNumber;
        incoming = new int[1];
        outgoing = new int[1];
    }
    /**
     * construct a node with a start capacity for edges of edgeCapacity.
     * @param nodeNumber
     * @param edgeCapacity
     */
    public LinkedListNode(int nodeNumber, int edgeCapacity) {
        this.key = nodeNumber;
        incoming = new int[edgeCapacity];
        outgoing = new int[edgeCapacity];
    }

    public void insertIncoming(int nodeNumber) {

        if ((nIncoming + 1) > incoming.length) {
            int expand = 2*nIncoming;
            incoming = Arrays.copyOf(incoming, expand);
        }
        incoming[nIncoming] = nodeNumber;
        nIncoming++;
    }

    public void insertOutgoing(int nodeNumber) {

        if ((nOutgoing + 1) > outgoing.length) {
            int expand = 2*nOutgoing;
            outgoing = Arrays.copyOf(outgoing, expand);
        }
        outgoing[nOutgoing] = nodeNumber;
        nOutgoing++;
    }

    public void removeIncoming(int nodeNumber) {
        int indexNumber = -1;
        for (int i = 0; i < incoming.length; i++) {
            if (incoming[i] == nodeNumber) {
                indexNumber = i;
                break;
            }
        }
        if (indexNumber > -1) {
            if (indexNumber < (nIncoming - 1)) {
                int nLen = nIncoming - indexNumber - 1;
                System.arraycopy(incoming, (indexNumber + 1), incoming, indexNumber, nLen);
            }
            nIncoming--;
        }
    }

    public void removeOutgoing(int nodeNumber) {
        int indexNumber = -1;
        for (int i = 0; i < outgoing.length; i++) {
            if (outgoing[i] == nodeNumber) {
                indexNumber = i;
                break;
            }
        }
        if (indexNumber > -1) {
            if (indexNumber < (nOutgoing - 1)) {
                int nLen = nOutgoing - indexNumber - 1;
                System.arraycopy(outgoing, (indexNumber + 1), outgoing, indexNumber, nLen);
            }
            nOutgoing--;
        }
    }
}
