package thirdparty.edu.princeton.cs.algs4;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

/**
 * adapted from QuadTree implementation from Sedgewick and Wayne
 * from 
 * http://algs4.cs.princeton.edu/92search/QuadTree.java.html
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
 */
public class QuadTreeInterval2D<T extends Comparable<T>, Value>  {
    
    private Node<T> root;

    // helper node data type
    private class Node<T extends Comparable<T>> {
        Interval2D<T> xy;
        Node<T> NW, NE, SE, SW;   // four subtrees
        Value value;           // associated data

        Node(Interval2D<T> box, Value value) {
            this.xy = box;
            this.value = value;
        }
    }


  /***********************************************************************
    *  Insert (x, y) into appropriate quadrant
    ***************************************************************************/
    public void insert(Interval2D<T> box, Value value) {
        root = insert(root, box, value);
    }

    private Node<T> insert(Node<T> h, Interval2D<T> box, Value value) {
        
        if (h == null) {
            return new Node<T>(box, value);
        }
        
        int cX = h.xy.intervalX.compareTo(box.intervalX);
        int cY = h.xy.intervalY.compareTo(box.intervalY);
        
        if ((cX < 0) && (cY < 0)) { 
            h.SW = insert(h.SW, box, value);
        } else if ((cX < 0) && (cY >= 0)) {
            h.NW = insert(h.NW, box, value);
        } else if ((cX >= 0) && (cY < 0)) {
            h.SE = insert(h.SE, box, value);
        } else if ((cX >= 0) && (cY >= 0)) { 
            h.NE = insert(h.NE, box, value);
        }
        
        return h;
    }

    public void remove(Interval2D<T> box) {
        List<Node<T>> parents = new ArrayList<Node<T>>();
        remove(root, box, parents);
    }
    
    private void remove(Node<T> h, Interval2D<T> box, 
        List<Node<T>> parents) {
        
        //TODO: need to simplify this pattern
        
        boolean isH = (h != null) && h.xy.equals(box);
        boolean isRoot = isH && h.equals(root);
        
        if ((h == null) || isRoot || isH) {
            
            if (parents.isEmpty() && !isRoot) {
                return;
            } else if (parents.size() > 1) {
                Node<T> parentParent = parents.get(
                    parents.size() - 2);
                Node<T> parent = parents.get(parents.size() - 1);
                if (parent.xy.equals(box)) {
                    removeNodeReattachChildren(parentParent, box);
                    return;
                } else if (isH) {
                    removeNodeReattachChildren(parent, box);
                    return;
                }
            }
            // for this case, parent is root node,
            // - if the box is the xy,
            //   then have to remove all children and
            //   reinsert them
            //   
            Node<T> parent = isRoot ? root : 
                parents.get(parents.size() - 1);
            if (parent.xy.equals(box)) {
                assert(parent.equals(root));
                List<Interval2D<T>> boxes = new
                    ArrayList<Interval2D<T>>();
                List<Value> values = new ArrayList<Value>();
                
                extractAllChildren(parent, boxes, values);
                root = null;
                for (int i = 0; i < boxes.size(); ++i) {
                    if (!boxes.get(i).equals(box)) {
                        insert(boxes.get(i), values.get(i));
                    }
                }
                return;
            }
        }
        
        parents.add(h);
        
        int cX = -1;
        int cY = -1;
        try {
            cX = h.xy.intervalX.compareTo(box.intervalX);
            cY = h.xy.intervalY.compareTo(box.intervalY);
        } catch (Throwable t) {
            int z = 1;
        }
        if ((cX < 0) && (cY < 0)) { 
            remove(h.SW, box, parents);
        } else if ((cX < 0) && (cY >= 0)) {
            remove(h.NW, box, parents);
        } else if ((cX >= 0) && (cY < 0)) {
            remove(h.SE, box, parents);
        } else if ((cX >= 0) && (cY >= 0)) { 
            remove(h.NE, box, parents);
        }
    }
    
    private void removeNodeReattachChildren(Node<T> parent, 
        Interval2D<T> rmBox) {

        if (parent.NW != null && parent.NW.xy.equals(rmBox)) {
            Node<T> node = parent.NW;
            parent.NW = getSingleNonNullChild(node);
            if (parent.NW != null) {
                return;
            }
            /*
               NW       NE
               nwe
               swe
            
               SW       SE
            */
            if (node.SE != null) {
                parent.NW = node.SE;
            }
            if (node.NW != null) {
                insert(parent, node.NW);
            }
            if (node.SW != null) {
                insert(parent, node.SW);
            }
            if (node.NE != null) {
                insert(parent, node.NE);
            }
        } else if (parent.NE != null && parent.NE.xy.equals(rmBox)) {
            Node<T> node = parent.NE;
            parent.NE = getSingleNonNullChild(node);
            if (parent.NE != null) {
                return;
            }
            /*
               NW       NE
                         nwe
                         swe
            
               SW       SE
            */
            if (node.SW != null) {
                parent.NE = node.SW;
            }
            if (node.NW != null) {
                insert(parent, node.NW);
            }
            if (node.NE != null) {
                insert(parent, node.NE);
            }
            if (node.SE != null) {
                insert(parent, node.SE);
            }
        } else if (parent.SW != null && parent.SW.xy.equals(rmBox)) {
            Node<T> node = parent.SW;
            parent.SW = getSingleNonNullChild(node);
            if (parent.SW != null) {
                return;
            }
            /*
               NW       NE
                         
               SW       SE
            nwe
            swe
            */
            if (node.NE != null) {
                parent.SW = node.NE;
            }
            if (node.NW != null) {
                insert(parent, node.NW);
            }
            if (node.SE != null) {
                insert(parent, node.SE);
            }
            if (node.SW != null) {
                insert(parent, node.SW);
            }
        } else if (parent.SE != null && parent.SE.xy.equals(rmBox)) {
            Node<T> node = parent.SE;
            parent.SE = getSingleNonNullChild(node);
            if (parent.SE != null) {
                return;
            }
            /*
               NW       NE
                         
               SW       SE
                         nwe
                         swe
            */
            if (node.NW != null) {
                parent.SE = node.NW;
            }
            if (node.NE != null) {
                insert(parent, node.NE);
            }
            if (node.SE != null) {
                insert(parent, node.SE);
            }
            if (node.SW != null) {
                insert(parent, node.SW);
            }
        } else if (parent.xy != null && parent.xy.equals(rmBox)) {
            int z = 1;
        }
    }
    
    private Node<T> insert(Node<T> h, Node<T> insNode) {
        
        if (h == null) {
            return insNode;
        }
        
        int cX = h.xy.intervalX.compareTo(insNode.xy.intervalX);
        int cY = h.xy.intervalY.compareTo(insNode.xy.intervalY);
        
        if ((cX < 0) && (cY < 0)) { 
            h.SW = insert(h.SW, insNode);
        } else if ((cX < 0) && (cY >= 0)) {
            h.NW = insert(h.NW, insNode);
        } else if ((cX >= 0) && (cY < 0)) {
            h.SE = insert(h.SE, insNode);
        } else if ((cX >= 0) && (cY >= 0)) { 
            h.NE = insert(h.NE, insNode);
        }
        
        return h;
    }
    
    /**
     * if there is only one non-null child, return that,
     * else return null.
     * @param node
     * @return 
     */
    private Node<T> getSingleNonNullChild(Node<T> node) {
        Node<T> nodeC = null;
        if (node.NW != null) {
            nodeC = node.NW;
        }
        if (node.NE != null) {
            if (nodeC != null) {
                return null;
            }
            nodeC = node.NE;
        }
        if (node.SW != null) {
            if (nodeC != null) {
                return null;
            }
            nodeC = node.SW;
        }
        if (node.SE != null) {
            if (nodeC != null) {
                return null;
            }
            nodeC = node.SE;
        }
        return nodeC;
    }
    
    /**
     * using pre-order traversal, find and extract all children 
     * from parent node.
     * @param node parent node
     * @param boxes
     * @param values 
     */
    private void extractAllChildren(Node<T> node, 
        List<Interval2D<T>> boxes, List<Value> values) {
        
        Set<Interval2D<T>> added = new HashSet<Interval2D<T>>();
        
        Stack<Node<T>> stack = new Stack<Node<T>>();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                if (!added.contains(node.xy)) {
                    boxes.add(node.xy);
                    values.add(node.value);
                }
                stack.push(node);
                if (node.NW != null) {
                    Node<T> node2 = node.NW;
                    node.NW = null;
                    node = node2;
                } else if (node.NE != null) {
                    Node<T> node2 = node.NE;
                    node.NE = null;
                    node = node2;
                } else if (node.SW != null) {
                    Node<T> node2 = node.SW;
                    node.SW = null;
                    node = node2;
                } else if (node.SE != null) {
                    Node<T> node2 = node.SE;
                    node.SE = null;
                    node = node2;
                } else {
                    node = null;
                }
            } else {
                node = stack.pop();
                if (node.NW != null) {
                    Node<T> node2 = node.NW;
                    node.NW = null;
                    node = node2;
                } else if (node.NE != null) {
                    Node<T> node2 = node.NE;
                    node.NE = null;
                    node = node2;
                } else if (node.SW != null) {
                    Node<T> node2 = node.SW;
                    node.SW = null;
                    node = node2;
                } else if (node.SE != null) {
                    Node<T> node2 = node.SE;
                    node.SE = null;
                    node = node2;
                } else {
                    node = null;
                }
            }
        }
    }

  /***********************************************************************
    *  Range search.
    ***************************************************************************/

    public List<Interval2D<T>> query2D(Interval2D<T> rect) {
        
        List<Interval2D<T>> output = new ArrayList<Interval2D<T>>();
        
        query2D(root, rect, output);
        
        return output;
    }

    private void query2D(Node<T> h, Interval2D<T> srch,
        List<Interval2D<T>> output) {
      
        /*
        TODO:
        consider improvements that lead to a balanced
        tree and hence faster queries.
        
        a search returns this which I havent read:
        "Improving the Performance of Region
        Quadtrees" by Wolfensberger
        http://www.ifi.uzh.ch/dam/jcr:ffffffff-96c1-007c-ffff-fffff2d50548/ReportWolfensbergerFA.pdf
        
        TODO:
        consider adding other methods:
        http://www.cs.cmu.edu/~rcm/papers/thesis/ch4.pdf
        
        */
        
        if (h == null) return;
        
        int cX = h.xy.intervalX.compareTo(srch.intervalX);
        int cY = h.xy.intervalY.compareTo(srch.intervalY);
        
        if ((cX == 0) && (cY == 0)) {
            output.add(h.xy);
        }
        
        if (h.SW != null && (cX <= 0) && (cY <= 0)) 
            query2D(h.SW, srch, output);
        if (h.NW != null && (cX <= 0) && (cY >= 0)) 
            query2D(h.NW, srch, output);
        if (h.SE != null && (cX >= 0) && (cY <= 0)) 
            query2D(h.SE, srch, output);
        if (h.NE != null && (cX >= 0) && (cY >= 0)) 
            query2D(h.NE, srch, output);    
    }

}
