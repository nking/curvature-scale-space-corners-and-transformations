package algorithms.connected;

import algorithms.misc.Misc;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
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
public class ConnectedPointsFinder {
    
    // key = groupIdx, value = pixels
    private TIntObjectMap<TIntSet> groupPixIdxsMap = new TIntObjectHashMap<TIntSet>();
    
    // key - pixIdx, value = groupIdx
    private TIntIntMap pixGroupMap = new TIntIntHashMap();
    
    protected boolean notValue = false;
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    /**
     * a list to hold each group as an item.
     */
    protected List<TIntSet> groupList = new ArrayList<TIntSet>();
    
    protected int minimumNumberInCluster = 3;
    
    private final int imgWidth;
    private final int imgHeight;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;

    public ConnectedPointsFinder(int imageWidth, int imageHeight) {
        
        imgWidth = imageWidth;
        
        imgHeight = imageHeight;
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
     * means is adjacent to another point in the group, making the group
     * contiguous.  The adjacency by default is using the 4 neighbor
     * pattern search unless the user has set that to 8 neighbors.
     * The runtime complexity is essentially O(pixIdxs.size()).
     * 
     * @param pixIdxs 
     */
    public void findConnectedPointGroups(TIntSet pixIdxs) {
            
        findClustersIterative(pixIdxs);
        
        prune();      
    }

    protected void findClustersIterative(TIntSet pixIdxs) {
        
        if (pixIdxs.isEmpty()) {
            return;
        }
                
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
    
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];
                
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {

            int uPoint = iter.next();
            
            ph.toPixelCoords(uPoint, imgWidth, xyout);
            
            int uY = xyout[1];
            int uX = xyout[0];
                        
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            boolean found = false;
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
                
                if (vX < 0 || vY < 0 || (vX >= imgWidth) || (vY >= imgHeight)) {
                    continue;
                }
            
                int vPoint = (int)ph.toPixelIndex(vX, vY, imgWidth);

                if (!pixIdxs.contains(vPoint)) {
                    continue;
                }
                
                //System.out.format("(%d,%d) (%d,%d)\n", uX, uY, vX, vY);

                processPair(uPoint, vPoint);
                
                found = true;
            }
            if (!found && minimumNumberInCluster == 1) {
                process(uPoint);
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
        
        int uGroupIdx = pixGroupMap.containsKey(uPoint) ?
            pixGroupMap.get(uPoint) : -1;
        int vGroupIdx = pixGroupMap.containsKey(vPoint) ?
            pixGroupMap.get(vPoint) : -1;
        int groupIdx;
        
        //System.out.format("u=%d v=%d sz=%d\n", uGroupIdx, vGroupIdx, 
        //    groupPixIdxsMap.size());
        
        if (uGroupIdx == -1 && vGroupIdx == -1) {
            TIntSet set = new TIntHashSet();
            set.add(uPoint);
            set.add(vPoint);
            groupIdx = pixGroupMap.size();
            pixGroupMap.put(uPoint, groupIdx);
            pixGroupMap.put(vPoint, groupIdx);
            groupPixIdxsMap.put(groupIdx, set);
        } else if (uGroupIdx == -1) {
            // add u to the vGroup
            groupIdx = vGroupIdx;
            TIntSet set = groupPixIdxsMap.get(groupIdx);
            set.add(uPoint);
            pixGroupMap.put(uPoint, groupIdx);
            groupPixIdxsMap.put(groupIdx, set);
        } else if (vGroupIdx == -1) {
            // add v to the uGroup
            groupIdx = uGroupIdx;
            TIntSet set = groupPixIdxsMap.get(groupIdx);
            set.add(vPoint);
            pixGroupMap.put(vPoint, groupIdx);
            groupPixIdxsMap.put(groupIdx, set);
        } else if (uGroupIdx != -1 && vGroupIdx != -1) {
            if (uGroupIdx == vGroupIdx) {
                return;
            }
            // merge the two
            TIntSet uSet = groupPixIdxsMap.get(uGroupIdx);
            TIntSet vSet = groupPixIdxsMap.get(vGroupIdx);
            if (vSet.size() < uSet.size()) {
                // add U to V
                groupPixIdxsMap.remove(uGroupIdx);
                TIntIterator iter = uSet.iterator();
                while (iter.hasNext()) {
                    int uIdx = iter.next();
                    pixGroupMap.put(uIdx, vGroupIdx);
                    vSet.add(uIdx);
                }
            } else {
                // add V to U
                groupPixIdxsMap.remove(vGroupIdx);
                TIntIterator iter = vSet.iterator();
                while (iter.hasNext()) {
                    int vIdx = iter.next();
                    pixGroupMap.put(vIdx, uGroupIdx);
                    uSet.add(vIdx);
                }
            }
        }
    }
    
    protected void process(int uPoint) {
        
        int uGroupIdx = pixGroupMap.containsKey(uPoint) ?
            pixGroupMap.get(uPoint) : -1;
        
        //System.out.format("u=%d sz=%d\n", uGroupIdx, groupPixIdxsMap.size());
        
        if (uGroupIdx == -1) {
            TIntSet set = new TIntHashSet();
            set.add(uPoint);
            int groupIdx = pixGroupMap.size();
            pixGroupMap.put(uPoint, groupIdx);
            groupPixIdxsMap.put(groupIdx, set);
        } else {
            // merge the two
            TIntSet uSet = groupPixIdxsMap.get(uGroupIdx);
            uSet.add(uPoint);
            pixGroupMap.put(uPoint, uGroupIdx);
        }
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
        
        groupList.clear();
                
        // key = groupIdx, value = set of pixels w/ group
        TIntObjectIterator<TIntSet> iter2 = groupPixIdxsMap.iterator();
       
        for (int i = 0; i < groupPixIdxsMap.size(); ++i) {
            
            iter2.advance();
            
            int gIdx = iter2.key();
            TIntSet pixIdxs = iter2.value();
            
            assert(pixIdxs != null);
            
            if (pixIdxs.size() >= minimumNumberInCluster) {
                groupList.add(pixIdxs);            
            }
        }
        
        groupPixIdxsMap.clear();
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
}
