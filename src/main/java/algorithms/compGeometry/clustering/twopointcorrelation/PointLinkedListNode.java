package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * a linked list node tailored for 2-pt correlation class: it holds 3 integer
 *    fields and a reference
 *
 * @author nichole
 */
public class PointLinkedListNode {

    protected int secondPointIndex = -1;

    protected int groupIndex = -1;

    protected int distanceIndex = -1;

    protected PointLinkedListNode next = null;

    public PointLinkedListNode() {}
    public PointLinkedListNode(int indexOfSecondPoint, int indexWithinGroupArray, int indexWithinDistanceArray) {
        this.secondPointIndex = indexOfSecondPoint;
        this.groupIndex = indexWithinGroupArray;
        this.distanceIndex = indexWithinDistanceArray;
    }

    public PointLinkedListNode insert(int indexOfSecondPoint, int indexWithinGroupArray, int indexWithinDistanceArray) {
        if (indexOfSecondPoint == -1) {
            throw new IllegalArgumentException("indexOfSecondPoint must be larger than -1");
        }
        if (this.secondPointIndex == -1) {
            this.secondPointIndex = indexOfSecondPoint;
            this.groupIndex = indexWithinGroupArray;
            this.distanceIndex = indexWithinDistanceArray;
            return this;
        }
        PointLinkedListNode node = new PointLinkedListNode(indexOfSecondPoint, indexWithinGroupArray, indexWithinDistanceArray);
        if (next == null) {
            next = node;
        } else {
        	PointLinkedListNode last = next;
            while (last.next != null) {
                last = last.next;
            }
            last.next = node;
        }
        return node;
    }

    public boolean isEmpty() {
        return (secondPointIndex == -1);
    }

    public PointLinkedListNode insertIfDoesNotAlreadyExist(int indexOfSecondPoint, int indexWithinGroupArray, int indexWithinDistanceArray) {
        if (indexOfSecondPoint == -1) {
            throw new IllegalArgumentException("indexOfSecondPoint must be larger than -1");
        }
        if (indexOfSecondPoint == this.secondPointIndex) {
            return null;
        }
        if (this.secondPointIndex == -1) {
            this.secondPointIndex = indexOfSecondPoint;
            this.groupIndex = indexWithinGroupArray;
            this.distanceIndex = indexWithinDistanceArray;
            return this;
        }

        if (next == null) {
            next = new PointLinkedListNode(indexOfSecondPoint, indexWithinGroupArray, indexWithinDistanceArray);
            return next;
        } else {
        	PointLinkedListNode last = next;
            if (last.secondPointIndex == indexOfSecondPoint) {
                return null;
            }
            while (last.next != null) {
                if (last.secondPointIndex == indexOfSecondPoint) {
                    return null;
                }
                last = last.next;
            }
            last.next = new PointLinkedListNode(indexOfSecondPoint, indexWithinGroupArray, indexWithinDistanceArray);
            return last.next;
        }
    }

    public void delete(PointLinkedListNode node) {

        if (secondPointIndex == -1) {
            return;
        }

        PointLinkedListNode last = null;
        PointLinkedListNode current = this;

        while (current != null) {
            if (current.equals(node)) {
                if (last == null) {
                    if (next != null) {
                        // this is the first node in the list.  reassign field values to next
                        this.secondPointIndex = next.secondPointIndex;
                        this.next = next.next;
                    } else {
                        // this is the first and only node
                        this.secondPointIndex = -1;
                    }
                } else {
                    last.next = current.next;
                }
                break;
            }
            last = current;
            current = current.next;
        }
    }

    public void delete(int deleteKey) {

        if (deleteKey == -1) {
            return;
        }

        PointLinkedListNode last = null;
        PointLinkedListNode current = this;

        while (current != null && current.secondPointIndex != -1) {
            if (current.secondPointIndex == deleteKey) {
                if (last == null) {
                    if (next != null) {
                        // this is the first node in the list.  reassign field values to next
                        this.secondPointIndex = next.secondPointIndex;
                        this.next = next.next;
                    } else {
                        // this is the first and only node
                        this.secondPointIndex = -1;
                    }
                } else {
                    // else, this node is not the first in the list, we can skip over it to delete it
                    last.next = current.next;
                }
                break;
            }
            last = current;
            current = current.next;
        }
    }

    public PointLinkedListNode search(int searchKey) {

        PointLinkedListNode latest = this;

        while (latest != null) {
            if (latest.secondPointIndex == searchKey) {
                return latest;
            }
            latest = latest.next;
        }
        return null;
    }
}
