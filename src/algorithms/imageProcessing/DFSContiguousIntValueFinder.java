package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class DFSContiguousIntValueFinder extends DFSContiguousValueFinder {
    
    protected final int[] imgValues;
    protected final int imgWidth;
    protected final int imgHeight;
    
    public DFSContiguousIntValueFinder(final GreyscaleImage input) {
        super(input);
        this.imgValues = null;
        imgWidth = input.getWidth();
        imgHeight = input.getHeight();
    }
    
    public DFSContiguousIntValueFinder(final GreyscaleImage input, 
        Set<PairInt> mask) {
        super(input, mask);
        this.imgValues = null;
        imgWidth = input.getWidth();
        imgHeight = input.getHeight();
    }
    
    public DFSContiguousIntValueFinder(int[] imgValues, 
        int imgWidth, int imgHeight) {
        super(null, null);
        
        this.imgValues = imgValues;
        this.imgWidth = imgWidth;
        this.imgHeight = imgHeight;
    }

    /**
     * NOTE: to keep the performance reasonable, the use of jvm stack has to
     * be reduced.  For default jvm settings, the size of the local array variable
     * within this method frame should be kept to having fewer that 128k items
     * in it, therefore the GreyscaleImage should be binned before use to keep
     * the number of pixels below 128k/1.5 (roughly keep below 87000 pixels).
     * If all pixels are connected, that limit has to be lowered.
     * 
     * @param pixelValue 
     */
    @Override
    protected void findClustersIterative(final int pixelValue) {
       
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
        
        java.util.Stack<Integer> stack = new java.util.Stack<Integer>();
        
        //O(N)
        for (int uIndex = (imgValues.length - 1); uIndex > -1; uIndex--) {
            int row = uIndex/imgWidth;
            int col = uIndex - (row * imgWidth);
                    
            PairInt p = new PairInt(col, row);
            if (!exclude.contains(p)) {
                
                int uPixValue = imgValues[uIndex];
            
                if ((notValue && (uPixValue != pixelValue)) ||
                    (!notValue && (uPixValue == pixelValue))) {
                    
                    stack.add(Integer.valueOf(uIndex));
                }
            }
        }
                
        while (!stack.isEmpty()) {

            int uIndex = stack.pop().intValue();
            
            Integer uKey = Integer.valueOf(uIndex);
            
            if (visited.contains(uKey)) {
                continue;
            }
            
            int uPixValue = imgValues[uIndex];
            
            if ((notValue && (uPixValue == pixelValue)) ||
                (!notValue && (uPixValue != pixelValue))) {
                
                visited.add(uKey);
                
                continue;
            }
            
            int uY = uIndex/imgWidth;
            int uX = uIndex - (uY * imgWidth);
            
            boolean foundANeighbor = false;

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
            
                if ((vX < 0) || (vX > (imgWidth - 1)) 
                    || (vY < 0) || (vY > (imgHeight - 1))) {
                    continue;
                }
                
                if (exclude.contains(new PairInt(vX, vY))) {
                    continue;
                }

                int vIndex = (vY * imgWidth) + vX;

                Integer vKey = Integer.valueOf(vIndex);

                int vPixValue = imgValues[vIndex];

                if ((notValue && (vPixValue == pixelValue)) ||
                    (!notValue && (vPixValue != pixelValue))) {

                    continue;
                }

                processPair(uKey, vKey, false);
                
                // inserting back at the top of the stack assures that the 
                // search continues next from an associated point
                stack.add(vKey);
                
                foundANeighbor = true;
            }
            
            if (!foundANeighbor && (minimumNumberInCluster == 1)) {                
                process(uKey, false);
            }
            
            visited.add(uKey);
        }
    }
    
    public PairIntArray getXY(int groupId) {
        
        if (groupMembership.size() == 0) {
            return new PairIntArray();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<Integer> indexes = getIndexes(groupId);
        
        PairIntArray xy = new PairIntArray(indexes.size());
                
        for (Integer index : indexes) {
                      
            int row = index.intValue()/imgWidth;
        
            int col = index.intValue() - (row * imgWidth);
        
            xy.add(col, row);
        }
        
        return xy;
    }
    
    public void getXY(final int groupId, final Set<PairInt> output) {
        
        if (groupMembership.isEmpty()) {
            return;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<Integer> indexes = groupMembership.get(groupId);

        for (Integer index : indexes) {
            
            int idx = index.intValue();
            
            int row = idx/imgWidth;
        
            int col = idx - (row * imgWidth);
            
            PairInt p = new PairInt(col, row);
                       
            output.add(p);
        }        
    }
    
    public Set<Integer> getIndexes(int groupId) {
        
        if (groupMembership.isEmpty()) {
            return new HashSet<Integer>();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
               
        return groupMembership.get(groupId);
    }
 
}
