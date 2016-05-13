package algorithms.graphs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 *
 * @author nichole
 */
public class RAGCSubGraph {
    
    protected final List<NormalizedCutsNode> nodes;
    
    protected final FlexCompRowMatrix diffOrSum;
    
    protected final Map<Integer, Set<Integer>> adjacentIndexes;
    
    public RAGCSubGraph(List<NormalizedCutsNode> theNodes, FlexCompRowMatrix
        matrixOfDiffOrSim, Map<Integer, Set<Integer>> adjacencyMap) {
        
        this.nodes = theNodes;
        
        this.diffOrSum = matrixOfDiffOrSim;
        
        this.adjacentIndexes = adjacencyMap;
    }
    
    public int getNumberOfNodes() {
        return nodes.size();
    }
    
    public RAGCSubGraph[] partition(boolean[] partitionCut) {
        
        if (partitionCut.length != nodes.size()) {
            throw new IllegalArgumentException(
                "partitionCute must be same length as nodes.size");
        }
        
        List<NormalizedCutsNode> nodes1 = new ArrayList<NormalizedCutsNode>();
        List<NormalizedCutsNode> nodes2 = new ArrayList<NormalizedCutsNode>();
        
        Map<Integer, Integer> map1 = new HashMap<Integer, Integer>();
        Map<Integer, Integer> map2 = new HashMap<Integer, Integer>();
        
        for (int i = 0; i < partitionCut.length; ++i) {
            
            Integer index1 = Integer.valueOf(i);
            
            if (partitionCut[i]) {
                nodes1.add(nodes.get(i));
                map1.put(index1, Integer.valueOf(map1.size()));
            } else {
                nodes2.add(nodes.get(i));
                map2.put(index1,Integer.valueOf(map2.size()));
            }
        }
        
        int n1 = map1.size();
        int n2 = map2.size();
        
        FlexCompRowMatrix w1 = new FlexCompRowMatrix(n1, n1);
        FlexCompRowMatrix w2 = new FlexCompRowMatrix(n2, n2);
        
        for (MatrixEntry entry : diffOrSum) {
            int row = entry.row();
            int col = entry.column();
            if (row > col) {
                continue;
            }
            double v = entry.get();
            Integer rKey1 = Integer.valueOf(row);
            Integer cKey1 = Integer.valueOf(col);
            
            if (map1.containsKey(rKey1)) {
                if (!map1.containsKey(cKey1)) {
                    continue;
                }
                Integer rKey2 = map1.get(rKey1);
                Integer cKey2 = map1.get(cKey1);
                w1.set(rKey2.intValue(), cKey2.intValue(), v);
                w1.set(cKey2.intValue(), rKey2.intValue(), v);
            } else {
                if (!map2.containsKey(cKey1)) {
                    continue;
                }
                Integer rKey2 = map2.get(rKey1);
                Integer cKey2 = map2.get(cKey1);
                w2.set(rKey2.intValue(), cKey2.intValue(), v);
                w2.set(cKey2.intValue(), rKey2.intValue(), v);
            }
        }
                
        Map<Integer, Set<Integer>> adjMap1 = new HashMap<Integer, Set<Integer>>();
        Map<Integer, Set<Integer>> adjMap2 = new HashMap<Integer, Set<Integer>>();
        for (int i = 0; i < 2; ++i) {
            Map<Integer, Integer> map;
             Map<Integer, Set<Integer>> adjMap;
            if (i == 0) {
                map = map1;
                adjMap = adjMap1;
            } else {
                map = map2;
                adjMap = adjMap2;
            }
            for (Entry<Integer, Integer> entry : map.entrySet()) {
                Integer indexOld = entry.getKey();
                Integer indexNew = entry.getValue();
                Set<Integer> indexes2Old = adjacentIndexes.get(indexOld);
                Set<Integer> indexes2New = new HashSet<Integer>();
                for (Integer index2Old : indexes2Old) {
                    Integer index2New = map.get(index2Old);
                    if (index2New != null) {
                        indexes2New.add(index2New);
                    }
                }
                adjMap.put(indexNew, indexes2New);
            }
        }        
        
        RAGCSubGraph g1 = new RAGCSubGraph(nodes1, w1, adjMap1);
        RAGCSubGraph g2 = new RAGCSubGraph(nodes2, w2, adjMap2);
        
        return new RAGCSubGraph[]{g1, g2};
    }

    public FlexCompRowMatrix getEdgeMatrix() {
        return diffOrSum;
    }
    
    public Set<Integer> getAdjacentIndexes(Integer index) {
        return adjacentIndexes.get(index);
    }
    
    public List<NormalizedCutsNode> getNodes() {
        return nodes;
    }
}
