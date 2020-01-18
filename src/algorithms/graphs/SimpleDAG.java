package algorithms.graphs;

/**
 *
 * @author nichole
 */
public class SimpleDAG {

    protected LinkedListNode[] vertices;
    protected int nVertices = 0;

    public SimpleDAG(boolean[][] connected) {
        createDAG(connected);
    }
    
    private void createDAG(boolean[][] connected) {
    
        vertices = new LinkedListNode[connected.length];

        for (int i = 0; i < connected.length; i++) {

            if (vertices[i] == null) {
                vertices[i] = new LinkedListNode(i);
            }

            for (int j = 0; j < connected[i].length; j++) {
                boolean con = connected[i][j];
                if (con) {

                    if (vertices[j] == null) {
                        vertices[j] = new LinkedListNode(j);
                    }

                    vertices[i].insertOutgoing(j);
                    vertices[j].insertIncoming(i);
                }
            }
        }

        nVertices = connected.length;
    }
    
    public LinkedListNode[] getNodes() {
        return vertices;
    }

    public boolean isEmpty() {
        return (nVertices == 0);
    }

    public void removeNode(int nodeNumber) {

        LinkedListNode node = vertices[nodeNumber];

        if (node != null) {
            // get incomingEdges list and remove this nodeNumber from each of their outgoing lists
            // get outgoingEdges list and remove this nodeNumber from each of their incoming lists

            for (int i = 0; i < node.incoming.length; i++) {
                int num = node.incoming[i];
                if (vertices[num] != null) {
                    vertices[num].removeOutgoing(nodeNumber);
                }
            }

            for (int i = 0; i < node.outgoing.length; i++) {
                int num = node.outgoing[i];
                if (vertices[num] != null) {
                    vertices[num].removeIncoming(nodeNumber);
                }
            }

            vertices[nodeNumber] = null;

            nVertices--;
        }
    }

}
