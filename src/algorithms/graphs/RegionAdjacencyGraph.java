package algorithms.graphs;

import algorithms.imageProcessing.ImageExt;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * class to hold a list of region nodes and calculate an adjacency map.
 * The classes that extend it hold the edge values.
 * 
 * @author nichole
 */
public class RegionAdjacencyGraph {
    
    //NOTE: may change o use more compact structures in future
    protected final List<Region> regions;
    
    protected final Map<Integer, Set<Integer>> adjacencyMap;
        
    protected final int imageWidth;
    protected final int imageHeight;
    
    //NOTE: these are transcribed to format [row][col]
    protected final int[][] labels;
        
    /**
     * constructor
     * @param img 
     * @param labels double array of labels for each pixel using the convention
     * labels[pixelIndex].  Note that the largest label must be less than 
     * the number of pixels in the image.
     * Also note that labels isn't copied and will be modified as the graph changes.
     */
    public RegionAdjacencyGraph(ImageExt img, int[] labels1D) {
        
        imageWidth = img.getWidth();
        imageHeight = img.getHeight();
        
        this.labels = new int[imageHeight][];
        for (int i = 0; i < imageHeight; ++i) {
            labels[i] = new int[imageWidth];
            for (int j = 0; j < imageWidth; ++j) {
                int pixIdx = img.getInternalIndex(j, i);
                labels[i][j] = labels1D[pixIdx];
            } 
        } 
        
        this.regions = createRegionsList(img, labels);
        
        this.adjacencyMap = createAdjacencyMap(this.regions);
    }
    
    /*
    public void mergeRegions(int regionIndex1, int regionIndex2) {

        Region region1 = regions.get(regionIndex1);
        
        Region region2 = regions.get(regionIndex2);
        
        // update region2 pixel labels to regionIndex1
        for (PairInt p : region2.getPoints()) {
            labels[p.getX()][p.getY()] = regionIndex1;
        }
        
        // update the regions
        region1.mergeIntoThis(region2);
        
        Integer index1 = Integer.valueOf(regionIndex1);
        Integer index2 = Integer.valueOf(regionIndex2);
        
        // update the adjacency map
        Set<Integer> indexes1 = adjacencyMap.get(index1);
        Set<Integer> indexes2 = adjacencyMap.get(index2);
        indexes1.addAll(indexes2);
        indexes1.remove(index2);
        for (Integer index3 : indexes2) {
            int idx3 = index3.intValue();
            if (regionIndex1 == idx3 || regionIndex2 == idx3) {
                continue;
            }
            Set<Integer> indexes4 = adjacencyMap.get(index3);
            if (indexes4 != null) {
                indexes4.remove(index2);
                indexes4.add(index1);
            }
        }
        adjacencyMap.remove(index2);
    }
    */
 
    public Map<Integer, Set<Integer>> createAdjacencyMap(List<Region> aRegion) {
        
        Map<Integer, Set<Integer>> map = new HashMap<Integer, Set<Integer>>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Map<PairInt, Integer> pToIMap = createPerimetersPointToIndexMap(aRegion);
        
        for (Map.Entry<PairInt, Integer> entry : pToIMap.entrySet()) {
            
            PairInt p = entry.getKey();
            
            Integer index = entry.getValue();
            
            Set<Integer> indexes = map.get(index);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                map.put(index, indexes);
            }
            
            int x = p.getX();
            int y = p.getY();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = pToIMap.get(p2);
                if ((index2 != null) && !index.equals(index2)) {
                    indexes.add(index2);
                }
            }
        }
        
        return map;
    }
    
    private Map<PairInt, Integer> createPerimetersPointToIndexMap(List<Region> regionsList) {
        
        Map<PairInt, Integer> map = new HashMap<PairInt, Integer>();
        
        for (int i = 0; i < regionsList.size(); ++i) {
            
            Integer index = Integer.valueOf(i);
            
            Region region = regionsList.get(i);
            
            for (PairInt p : region.getPerimeter()) {
                map.put(p, index);
            }
        }
        
        return map;
    }
    
    /**
     * @param img
     * @param labels array of format [xcoord][ycoord] = label where label 
     * is less than the number of pixels in an image.
     * @return 
     */
    private List<Region> createRegionsList(ImageExt img, int[][] labels) {
        
        int nPix = img.getNPixels();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int maxLabel = Integer.MIN_VALUE;
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = labels[j][i];
                if (label > maxLabel) {
                    maxLabel = label;
                }
            }
        }
        
        Map<Integer, Set<PairInt>> map = createRegionsMap(img, labels);
        
        List<Region> regionList = new ArrayList<Region>();
        for (int i = 0; i <= maxLabel; ++i) {
            Integer label = Integer.valueOf(i);
            Set<PairInt> set = map.get(label);
            if (set == null) {
                regionList.add(new Region(new HashSet<PairInt>()));
            } else {
                regionList.add(new Region(set));
            }
        }
        
        return regionList;
    }
    
    private Map<Integer, Set<PairInt>> createRegionsMap(ImageExt img, int[][] labels) {
        
        Map<Integer, Set<PairInt>> map = new HashMap<Integer, Set<PairInt>>();
        
        int w = img.getWidth();
        int h = img.getHeight();
                
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                int label = labels[j][i];
                Integer index = Integer.valueOf(label);
                
                Set<PairInt> set = map.get(index);
                
                if (set == null) {
                    set = new HashSet<PairInt>();
                    map.put(index, set);
                }
                set.add(new PairInt(i, j));
            }
        }
        
        return map;
    }

    public Set<Integer> getAdjacentIndexes(Integer index) {
        return adjacencyMap.get(index);
    }
    
    public int getNumberOfRegions() {
        return regions.size();
    }
}
