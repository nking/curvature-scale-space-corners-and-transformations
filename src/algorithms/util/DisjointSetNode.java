package algorithms.util;

public class DisjointSetNode {

 // the set member data to be held in a disjoint set
    Object member;
    DisjointSetNode next;

    // pointer to the set representative
    DisjointSetNode representative;

    public Object getMember() {
        return member;
    }

    public void setMember(Object member) {
        this.member = member;
    }

    public DisjointSetNode getNext() {
        return next;
    }

    public void setNext(DisjointSetNode next) {
        this.next = next;
    }

    public DisjointSetNode getRepresentative() {
        return representative;
    }

    public void setRepresentative(DisjointSetNode representative) {
        this.representative = representative;
    }
}
