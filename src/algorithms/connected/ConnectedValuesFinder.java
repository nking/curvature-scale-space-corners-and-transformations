package algorithms.connected;

import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * given a set of points, finds the connected among them and
 * places them into groups.
 * note that connected here means adjacent to one another and
 * adjacent is defined by the default "4 neighbor" offsets,
 * but can be overridden to use all 8 neighbors.
 * 
 * The runtime complexity is essentially O(N_points).
 * 
 * @author nichole
 */
public class ConnectedValuesFinder {
    
    private final DisjointSet2Helper disjointSetHelper;
    
    // key = pixIdx, value = disjoint set node with key pixIdx
    private final TIntObjectMap<DisjointSet2Node<Integer>> pixNodes;

    protected boolean notValue = false;
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    /**
     * a list to hold each group as an item.
     */
    protected List<TIntSet> groupList = null;
    
    protected int minimumNumberInCluster = 3;
    
    protected final GreyscaleImage img;
    
    protected final TIntSet exclude;
    
    private final int imgWidth;
    private final int imgHeight;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;

    public ConnectedValuesFinder(final GreyscaleImage input) {
        
        this.img = input;
        
        this.log = Logger.getLogger(this.getClass().getName());
        
        this.exclude = new TIntHashSet();
        
        imgWidth = input.getWidth();
        
        imgHeight = input.getHeight();
        
        pixNodes = new TIntObjectHashMap<DisjointSet2Node<Integer>>();
    
        disjointSetHelper = new DisjointSet2Helper();
    }
    
    public ConnectedValuesFinder(final GreyscaleImage input, TIntSet mask) {
        
        this.img = input;
        
        this.log = Logger.getLogger(this.getClass().getName());
        
        this.exclude = new TIntHashSet(mask);
        
        imgWidth = input.getWidth();
        
        imgHeight = input.getHeight();
        
        pixNodes = new TIntObjectHashMap<DisjointSet2Node<Integer>>();
    
        disjointSetHelper = new DisjointSet2Helper();
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    /**
     * find the groups of connected points in pixIdxs where connected
     * means is adjacent to another point in the group and having this 
     * pixelValue.  The adjacency by default is using the 4 neighbor
     * pattern search unless the user has set that to 8 neighbors.
     * The runtime complexity is essentially O(pixIdxs.size()).
     * 
     * @param pixelValue 
     */
    public void findGroups(int pixelValue) {
    
        notValue = false;
        
        initMap();
        
        findClustersIterative(pixelValue);
        
        prune();        
    }
    
    /**
     * find the groups of connected points in pixIdxs where connected
     * means is adjacent to another point in the group and re any points that do
     * not have this pixelValue.  The adjacency by default is using the 4 neighbor
     * pattern search unless the user has set that to 8 neighbors.
     * The runtime complexity is essentially O(pixIdxs.size()).
     * 
     * @param pixelValue 
     */
    public void findGroupsNotThisValue(int pixelValue) {
    
        notValue = true;
        
        initMap();
        
        findClustersIterative(pixelValue);
        
        prune();        
    }

    protected void findClustersIterative(int pixelValue) {
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            /*
            for 8 neighbor, can use 4 offsets instead of 8 if visiting all pix
            
             2  *  *  *       (2,0) 1:2,1:1,2:1
             1  *  *  +       (1,1) 0:1,0:0,1:0,2:0,2:1,2:2,1:2,0:2
             0  *  *  *       (2,1) 1:1,1:0,2:0,3:0,3:1,3:2,2:2,1:2
                0  1  2             X: 1:1,1:0,2:0, 1:2
                                    use: +1,-1  +1,0  +1,+1  0:1
            */
            dxs = new int[]{1,  1, 1, 0};
            dys = new int[]{-1, 0, 1, 1};
        }        
                
        for (int uPoint = 0; uPoint < img.getNPixels(); ++uPoint) {
            
            if (exclude.contains(uPoint)) {
                continue;
            }
            
            int uPixValue = img.getValue(uPoint);
            
            if ((notValue && (uPixValue == pixelValue)) ||
                (!notValue && (uPixValue != pixelValue))) {
                                
                continue;
            }
            
            int uY = uPoint/imgWidth;
            int uX = uPoint - (uY * imgWidth);
                                    
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
                
                if (vX < 0 || vY < 0 || (vX >= (imgWidth - 1))
                    || (vY >= (imgHeight - 1))) {
                    continue;
                }
            
                int vPoint = (vY * imgWidth) + vX;

                if (exclude.contains(vPoint)) {
                    continue;
                }
                
                int vPixValue = img.getValue(vPoint);

                if ((notValue && (vPixValue == pixelValue)) ||
                    (!notValue && (vPixValue != pixelValue))) {

                    continue;
                }

                processPair(uPoint, vPoint);
            }            
        }
    }
    
    /*
    private int[] debugCoords(int pixIdx) {
        int y = pixIdx/this.imgWidth;
        int x = pixIdx - (y * imgWidth);
        return new int[]{x, y};
    }*/
  
    protected void processPair(int uPoint, int vPoint) {
        
        DisjointSet2Node<Integer> uNode = pixNodes.get(uPoint); 
        DisjointSet2Node<Integer> uParentNode 
            = disjointSetHelper.findSet(uNode);
        assert(uParentNode != null);
        
        int uGroupId = uParentNode.getMember().intValue();

        DisjointSet2Node<Integer> vNode = pixNodes.get(vPoint); 
        DisjointSet2Node<Integer> vParentNode 
            = disjointSetHelper.findSet(vNode);
        assert(vParentNode != null);
        
        int vGroupId = vParentNode.getMember().intValue();
        
        DisjointSet2Node<Integer> merged = 
            disjointSetHelper.union(uParentNode, vParentNode);
        
        pixNodes.put(uGroupId, merged);
        
        pixNodes.put(vGroupId, merged);
    }

    public List<TIntSet> getGroupMembershipList() {
        return groupList;
    }

    public int getNumberOfGroups() {
        return groupList.size();
    }

    /**
     * gather groups and remove those smaller than minimumNumberInCluster
     */
    protected void prune() {
        
        // key = repr node index, value = set of pixels w/ repr
        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<DisjointSet2Node<Integer>> iter =
            pixNodes.iterator();
        for (int i = 0; i < pixNodes.size(); ++i) {
            
            iter.advance();
            
            int pixIdx = iter.key();
            DisjointSet2Node<Integer> node = iter.value();
            
            DisjointSet2Node<Integer> repr = disjointSetHelper.findSet(node);
            
            int reprIdx = repr.getMember().intValue();
            
            TIntSet set = map.get(reprIdx);
            if (set == null) {
                set = new TIntHashSet();
                map.put(reprIdx, set);
            }
            set.add(pixIdx);
        }
        
        log.finest("number of groups before prune=" + map.size());
        
        // rewrite the above into a list
        List<TIntSet> groups = new ArrayList<TIntSet>();
        
        TIntObjectIterator<TIntSet> iter2 = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            iter2.advance();
                        
            TIntSet idxs = iter2.value();
            
            if (idxs.size() >= minimumNumberInCluster) {
                groups.add(idxs);
            }
        }
        
        this.groupList = groups;
        
        log.finest("number of groups after prune=" + groups.size());
    }
    
    public TIntIntMap createPointIndexMap() {
        
        TIntIntMap ptIdxMap = new TIntIntHashMap();
        
        int n = getNumberOfGroups();
        for (int i = 0; i < n; ++i) {
            
            TIntSet set = getXY(i);
            
            TIntIterator iter = set.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ptIdxMap.put(pixIdx, i);
            }
        }
        
        return ptIdxMap;
    }


    public int getNumberofGroupMembers(int groupId) {
        if (groupList.isEmpty()) {
            return 0;
        }
        if (groupId > (groupList.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
                + " is outside of range of nGroups=" + groupList.size());
        }
        return groupList.get(groupId).size();
    }

    public TIntSet getXY(int groupId) {
        if (groupList.isEmpty()) {
            return new TIntHashSet();
        }
        if (groupId > (groupList.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
                + " is outside of range of nGroups=" + groupList.size());
        }
        TIntSet set = groupList.get(groupId);
        return set;
    }

    private void initMap() {
        
        pixNodes.clear();
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); ++pixIdx) {
                        
            if (exclude.contains(pixIdx)) {
                continue;
            }
            
            Integer index = Integer.valueOf(pixIdx);
            
            DisjointSet2Node<Integer> pNode =
                disjointSetHelper.makeSet(
                    new DisjointSet2Node<Integer>(index));
            
            pixNodes.put(pixIdx, pNode);
        }
    }
    
}
