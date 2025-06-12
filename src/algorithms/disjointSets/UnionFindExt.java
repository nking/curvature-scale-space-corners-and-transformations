package algorithms.disjointSets;

public class UnionFindExt extends UnionFind {
    public UnionFindExt(UnionFind uf) {
        super(uf.getParent().length);
        System.arraycopy(uf.getParent(), 0, this.getParent(), 0, uf.getParent().length);
        System.arraycopy(uf.rank, 0, this.rank, 0, uf.rank.length);
    }
}
