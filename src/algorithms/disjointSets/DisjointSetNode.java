package algorithms.util;

public class DisjointSetNode<T> {

 // the set member data to be held in a disjoint set
    protected T member;
    protected DisjointSetNode<T> next;

    // pointer to the set representative
    protected DisjointSetNode<T> representative;

    public T getMember() {
        return member;
    }

    public void setMember(T member) {
        this.member = member;
    }

    public DisjointSetNode<T> getNext() {
        return next;
    }

    public void setNext(DisjointSetNode<T> next) {
        this.next = next;
    }

    public DisjointSetNode<T> getRepresentative() {
        return representative;
    }

    public void setRepresentative(DisjointSetNode<T> representative) {
        this.representative = representative;
    }
}
