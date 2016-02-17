package algorithms.imageProcessing.features;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * class to encapsulate methods to make a rough skeletonization of a set of
 * points (using line thinning) and methods to find the closest skeleton
 * point to a query point.  The class currently expects to handle a list of
 * such points to make an internal resulting list of skeletons 
 * (a.k.a. medial axes for each point set in a list).
 * 
 * (for example, this instance is constructed from points which are a larger
 * set of blobs in a hierarchy of segmentation.  this instance and the bounds
 * are used to merge the groups of points from a finer segmentation into 
 * the groups defined by larger bounds.  the resulting merged list is then
 * ordered to have same order as the internal data here and then both are
 * filtered to remove empty sets.  thereafter, the medial axes can be used
 * complementarily with the grouped point lists for functions that need to know
 * the direction of "inward" for the group of points.)
 * 
 * @author nichole
 */
public class BlobMedialAxes implements Serializable {
    
    private static final long serialVersionUID = 987123L;
    
    protected List<Map<Integer, List<Integer>>> skeletonXMapList = null;
    protected List<Map<Integer, List<Integer>>> skeletonYMapList = null;
    protected final List<PairInt> xyCentroids;
    protected final List<Double> labLColorList;
    protected final List<Double> labAColorList;
    protected final List<Double> labBColorList;
    protected final List<Double> rColorList;
    protected final List<Double> gColorList;
    protected final List<Double> bColorList;
    
    public BlobMedialAxes(final List<Set<PairInt>> blobs, 
        final List<Double> labLClrList, final List<Double> labAClrList,
        final List<Double> labBClrList, final List<Double> rClrList,
        final List<Double> gClrList, final List<Double> bClrList) {
        
        int n = blobs.size(); 
        xyCentroids = new ArrayList<PairInt>(n);
        skeletonXMapList = new ArrayList<Map<Integer, List<Integer>>>(n);
        skeletonYMapList = new ArrayList<Map<Integer, List<Integer>>>(n);
        
        this.labLColorList = new ArrayList<Double>(labLClrList);
        this.labAColorList = new ArrayList<Double>(labAClrList);
        this.labBColorList = new ArrayList<Double>(labBClrList);
        
        this.rColorList = new ArrayList<Double>(rClrList);
        this.gColorList = new ArrayList<Double>(gClrList);
        this.bColorList = new ArrayList<Double>(bClrList);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        
        // make skeleton for each blob for detailed perimeter "inward" directions
        for (int i = 0; i < blobs.size(); ++i) {
            
            Set<PairInt> blob = blobs.get(i);
            Set<PairInt> skeleton = new HashSet<PairInt>(blob);
            
            int[] minMaxXY = MiscMath.findMinMaxXY(skeleton);
            lt.applyLineThinner(skeleton, minMaxXY[0], minMaxXY[1], minMaxXY[2],
                minMaxXY[3]);
            int[] xPoints = new int[skeleton.size()];
            int[] yPoints = new int[skeleton.size()];

            int count = 0;
            for (PairInt p : skeleton) {
                xPoints[count] = p.getX();
                yPoints[count] = p.getY();
                count++;
            }
            
            double[] xyCen = curveHelper.calculateXYCentroids(blob);
            
            // order skeleton by x and then y
            Map<Integer, List<Integer>> xSkeletonMap = makeXMap(xPoints, yPoints);
            Map<Integer, List<Integer>> ySkeletonMap = makeYMap(xPoints, yPoints);
            
            xyCentroids.add(new PairInt((int)Math.round(xyCen[0]),
                (int)Math.round(xyCen[1])));
            skeletonXMapList.add(xSkeletonMap);
            skeletonYMapList.add(ySkeletonMap);
        }
    }
    
    public int getNumberOfItems() {
        return xyCentroids.size();
    }
    
    /**
     * find the closest skeleton point (that is, the medial axis point) to
     * the coordinates (x, y) in the blob at index.
     * @param index
     * @param x
     * @param y
     * @return 
     */
    public PairInt findClosestPoint(int index, int x, int y) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return findClosestPoint(x, y, skeletonXMapList.get(index), 
            skeletonYMapList.get(index));
    }
    
    protected PairInt findClosestPoint(int x, int y, 
        Map<Integer, List<Integer>> xSkeletonMap,
        Map<Integer, List<Integer>> ySkeletonMap) {
        
        List<Integer> ys = xSkeletonMap.get(Integer.valueOf(x));
        List<Integer> xs = ySkeletonMap.get(Integer.valueOf(y));
        
        if (xs == null && ys == null) {
            // TODO: need to replace w/ improved data structure.
            // brute force search for closest.
            int minDist = Integer.MAX_VALUE;
            int xSkel3 = -1;
            int ySkel3 = -1;
            for (Entry<Integer, List<Integer>> entry : xSkeletonMap.entrySet()) {
                Integer x0 = entry.getKey();
                for (Integer y0 : entry.getValue()) {
                    int diffX = x0.intValue() - x;
                    int diffY = y0.intValue() - y;
                    int distSq = (diffX * diffX + diffY * diffY);
                    if (distSq < minDist) {
                        minDist = distSq;
                        xSkel3 = x0.intValue();
                        ySkel3 = y0.intValue();
                    }
                }
            }
            return new PairInt(xSkel3, ySkel3);
        }
        
        final Integer xSkel1 = Integer.valueOf(x);
        Integer ySkel1 = null;
        if (ys != null) {
            int minDistY2 = Integer.MAX_VALUE;
            for (Integer y0 : ys) {
                int distYSq = y0.intValue() - y;
                distYSq *= distYSq;
                if (distYSq < minDistY2) {
                    minDistY2 = distYSq;
                    ySkel1 = y0;
                }
            }
        }
        
        Integer xSkel2 = null;
        final Integer ySkel2 = Integer.valueOf(y);
        if (xs != null) {
            int minDistX2 = Integer.MAX_VALUE;
            for (Integer x0 : xs) {
                int distXSq = x0.intValue() - x;
                distXSq *= distXSq;
                if (distXSq < minDistX2) {
                    minDistX2 = distXSq;
                    xSkel2 = x0;
                }
            }
        }
        
        if (ySkel1 != null && xSkel2 != null) {
            int diffX1 = xSkel1.intValue() - x;
            int diffY1 = ySkel1.intValue() - y;
            int diffX2 = xSkel2.intValue() - x;
            int diffY2 = ySkel2.intValue() - y;
            if (((diffX1 * diffX1) + (diffY1 * diffY1)) <
                ((diffX2 * diffX2) + (diffY2 * diffY2))) {
                return new PairInt(xSkel1, ySkel1);
            } else {
                return new PairInt(xSkel2, ySkel2);
            }
        } else if (ySkel1 != null) {
            return new PairInt(xSkel1, ySkel1);
        }
        
        return new PairInt(xSkel2, ySkel2);
    }
    
    public PairInt getOriginalBlobXYCentroid(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return xyCentroids.get(index).copy();
    }
    
    /**
     * get the LAB color space averaged colors for the points within the 
     * bounds at index in lists.
     * @param index
     * @return 
     */
    public float[] getLABColors(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        float[] c = new float[3];
        c[0] = labLColorList.get(index).floatValue();
        c[1] = labAColorList.get(index).floatValue();
        c[2] = labBColorList.get(index).floatValue();
        
        return c;
    }
    
    /**
     * get red color
     * @param index
     * @return 
     */
    public double getR(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return rColorList.get(index);
    }
    
    /**
     * get green color
     * @param index
     * @return 
     */
    public double getG(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return gColorList.get(index);
    }
    
    /**
     * get blue color
     * @param index
     * @return 
     */
    public double getB(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return bColorList.get(index);
    }
    
    /**
     * map w/ key = x and y = ascending ordered values for that x
     * @param x
     * @param y
     * @return 
     */
    protected Map<Integer, List<Integer>> makeXMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < x.length; ++i) {
            Integer key = Integer.valueOf(x[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(y[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> ys = entry.getValue();
            if (ys.size() > 1) {
                Collections.sort(ys);
            }
        }
        
        return map;
    }
    
    /**
     * map w/ key = y and x = ascending ordered values for that y
     * @param x
     * @param y
     * @return 
     */
    protected Map<Integer, List<Integer>> makeYMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < y.length; ++i) {
            Integer key = Integer.valueOf(y[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(x[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> xs = entry.getValue();
            if (xs.size() > 1) {
                Collections.sort(xs);
            }
        }
        
        return map;
    }

    /**
     * rewrite the internal data structures to remove the subset of indexes
     * given.  Note that the internal skeleton
     * maps are for now, not trimmed to the smaller subset of points as the
     * larger set - that is not harmful to the internal functions.
     * @param removeIndexes an ascending list of unique indexes to remove
     */
    public void removeIndexes(List<Integer> removeIndexes) {
        
        // NOTE: should consider retaining blobs in constructor so can
        // rebuild the skeleton maps here from the reduced set of points.
        // for now, deciding that the extra information in the skeleton maps 
        // doesn't harm the use of the structures.
        
        Set<Integer> exclude = new HashSet<Integer>(removeIndexes);
        
        int n2 = xyCentroids.size() - removeIndexes.size();
        
        List<Double> lLab = new ArrayList<Double>();
        List<Double> aLab = new ArrayList<Double>();
        List<Double> bLab = new ArrayList<Double>();
        List<Double> r = new ArrayList<Double>();
        List<Double> g = new ArrayList<Double>();
        List<Double> b = new ArrayList<Double>();
        List<PairInt> xyCentroid2 = new ArrayList<PairInt>();
        
        for (int i = 0; i < xyCentroids.size(); ++i) {
            
            Integer index = Integer.valueOf(i);
            if (!exclude.contains(index)) {
                xyCentroid2.add(xyCentroids.get(i));
                lLab.add(labLColorList.get(i));
                aLab.add(labAColorList.get(i));
                bLab.add(labBColorList.get(i));
                r.add(rColorList.get(i));
                g.add(gColorList.get(i));
                b.add(bColorList.get(i));
            }
        }   
        this.xyCentroids.clear();
        this.xyCentroids.addAll(xyCentroid2);
                
        this.labLColorList.clear();
        this.labLColorList.addAll(lLab);
        
        this.labAColorList.clear();
        this.labAColorList.addAll(aLab);
        
        this.labBColorList.clear();
        this.labBColorList.addAll(bLab);
        
        this.rColorList.clear();
        this.rColorList.addAll(r);
        
        this.gColorList.clear();
        this.gColorList.addAll(g);
        
        this.bColorList.clear();
        this.bColorList.addAll(b);
    }
    
    public void peristToFile(String filePath) throws IOException {
        FileOutputStream fos = null;
        try {
            fos = new FileOutputStream(filePath);
            peristToFile(fos);
        } catch (IOException e) {
            throw new IOException(e);
        } finally {
            if (fos != null) {
                fos.close();
            }
        }
    }
    
    public void peristToFile(FileOutputStream fos)throws IOException {
        
        ObjectOutputStream oos = null;
        try {
            oos = new ObjectOutputStream(fos);
            peristToFile(oos);
        } catch (IOException e) {
            throw new IOException(e);
        } finally {
            if (oos != null) {
                oos.close();
            }
        }
    }
    
    public void peristToFile(ObjectOutputStream oos)throws IOException {
        writeObject(oos);
        oos.flush();
    }
    
    public BlobMedialAxes(ObjectInputStream ois) throws IOException, ClassNotFoundException {
        
        this.labLColorList = new ArrayList<Double>();
        this.labAColorList = new ArrayList<Double>();
        this.labBColorList = new ArrayList<Double>();
        this.rColorList = new ArrayList<Double>();
        this.gColorList = new ArrayList<Double>();
        this.bColorList = new ArrayList<Double>();
        this.xyCentroids = new ArrayList<PairInt>();
        
        readObject(ois);
    }
    
    public BlobMedialAxes(String filePath) throws IOException {
        
        this.labLColorList = new ArrayList<Double>();
        this.labAColorList = new ArrayList<Double>();
        this.labBColorList = new ArrayList<Double>();
        this.rColorList = new ArrayList<Double>();
        this.gColorList = new ArrayList<Double>();
        this.bColorList = new ArrayList<Double>();
        this.xyCentroids = new ArrayList<PairInt>();
        
        FileInputStream fis = null;
        ObjectInputStream ois = null;
        try {
            fis = new FileInputStream(filePath);
            ois = new ObjectInputStream(fis);
            readObject(ois);
        } catch (IOException e) {
            throw new IOException(e);
        } catch (ClassNotFoundException ex) {
            throw new IOException(ex);
        } finally {
            if (fis != null) {
                fis.close();
            }
            if (ois != null) {
                ois.close();
            }
        }
    }
    
    private void writeObject(java.io.ObjectOutputStream out) throws IOException {
        
        int nMapX = (skeletonXMapList == null) ? 0 : skeletonXMapList.size();
        int nMapY = (skeletonYMapList == null) ? 0 : skeletonYMapList.size();
        int nXYC = (xyCentroids == null) ? 0 : xyCentroids.size();
        int nColorList = labLColorList.size();
        
        out.writeInt(nMapX);
        out.writeInt(nMapY);
        out.writeInt(nXYC);
        out.writeInt(nColorList);
        
        writeSkeletonMap(out, skeletonXMapList);
        writeSkeletonMap(out, skeletonYMapList);
        writeCentroids(out, xyCentroids);
        
        writeDoubleList(out, labLColorList);
        writeDoubleList(out, labAColorList);
        writeDoubleList(out, labBColorList);
        
        out.writeInt(rColorList.size());
        
        writeDoubleList(out, rColorList);
        writeDoubleList(out, gColorList);
        writeDoubleList(out, bColorList);
    }
    
    protected void readObject(java.io.ObjectInputStream in) throws IOException, 
        ClassNotFoundException {
                
        int nMapX = in.readInt();
        int nMapY = in.readInt();
        int nXYC = in.readInt();
        int nColorList = in.readInt();
        
        this.skeletonXMapList = readSkeletonMap(in, nMapX);
        this.skeletonYMapList = readSkeletonMap(in, nMapY);
        
        xyCentroids.clear();
        xyCentroids.addAll(readCentroids(in, nXYC));
        
        labLColorList.clear();
        labAColorList.clear();
        labBColorList.clear();
        
        List<Double> a = readDoubleList(in, nColorList);
        labLColorList.addAll(a);
        a = readDoubleList(in, nColorList);
        labAColorList.addAll(a);
        a = readDoubleList(in, nColorList);
        labBColorList.addAll(a);
        
        try {
            
            rColorList.clear();
            gColorList.clear();
            bColorList.clear();
            
            int nOs = in.readInt();

            a = readDoubleList(in, nOs);
            rColorList.addAll(a);
            a = readDoubleList(in, nOs);
            gColorList.addAll(a);
            a = readDoubleList(in, nOs);
            bColorList.addAll(a);
        } catch (EOFException e) {
            // some existing persisted files do not have the o1,o2,o3 colors
        }
    }
  
    private void writeSkeletonMap(java.io.ObjectOutputStream out,
        List<Map<Integer, List<Integer>>> mapList) throws IOException {
        
        if (mapList != null) {
            for (int i = 0; i < mapList.size(); ++i) {
                
                // write index i and map size
                Map<Integer, List<Integer>> map = mapList.get(i);
                out.writeInt(i);
                out.writeInt(map.size());
                
                for (Entry<Integer, List<Integer>> entry : map.entrySet()) {
                    
                    // write key and list size then list items
                    Integer key = entry.getKey();
                    List<Integer> value = entry.getValue();
                    
                    out.writeInt(key);
                    out.writeInt(value.size());
                    
                    for (Integer index : value) {
                        out.writeInt(index.intValue());
                    }
                }
            }
        }
    }
    
    private List<Map<Integer, List<Integer>>> readSkeletonMap(
        java.io.ObjectInputStream in, int listSize) throws IOException {
        
        List<Map<Integer, List<Integer>>> mapsList = 
            new ArrayList<Map<Integer, List<Integer>>>(listSize);
        
        if (listSize == 0) {
            return mapsList;
        }
        
        for (int i = 0; i < listSize; ++i) {
            
            // read index i and map size
            int idx = in.readInt();
            int mapSize = in.readInt(); 
            
            assert(idx == i);
                        
            Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>(mapSize);

            for (int j = 0; j < mapSize; ++j) {
                
                // read each key, each value list size, then each value list item
                int key = in.readInt();
                int valueSize = in.readInt();
                
                List<Integer> values = new ArrayList<Integer>(valueSize);
                map.put(Integer.valueOf(key), values);
                
                for (int k = 0; k < valueSize; ++k) {
                    int value = in.readInt();
                    values.add(Integer.valueOf(value));
                }
            }
            
            mapsList.add(map);
        }
        
        assert(mapsList.size() == listSize);
        
        return mapsList;
    }
   
    protected void writeCentroids(ObjectOutputStream out, List<PairInt> xyCen) 
        throws IOException {
        
        if (xyCen == null) {
            return;
        }
        
        int n = xyCen.size();
        
        for (int i = 0; i < n; ++i) {
            PairInt p = xyCen.get(i);
            out.writeInt(p.getX());
            out.writeInt(p.getY());
        }
        
    }
    
    private List<PairInt> readCentroids(java.io.ObjectInputStream in, int length) throws IOException {
        
        List<PairInt> xyCen = new ArrayList<PairInt>();
        
        // write the length of the array
        if (length == 0) {
            return xyCen;
        }
                
        for (int i = 0; i < length; ++i) {
            int x = in.readInt();
            int y = in.readInt();
            xyCen.add(new PairInt(x, y));
        }
        
        return xyCen;
    }

    private void writeDoubleList(ObjectOutputStream out, List<Double> list) throws IOException {
                
        for (Double d : list) {
            out.writeDouble(d.doubleValue());
        }
    }
    
    private List<Double> readDoubleList(java.io.ObjectInputStream out, int listSize) 
        throws IOException {

        List<Double> list = new ArrayList<Double>(listSize);
        
        for (int i = 0; i < listSize; ++i) {
            double d = out.readDouble();
            list.add(Double.valueOf(d));
        }
        
        return list;
    }
    
}
