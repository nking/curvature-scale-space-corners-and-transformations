package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.imageProcessing.Sky.SkyObject;
import algorithms.imageProcessing.features.mser.EllipseHelper;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * NOT READY FOR USE.
 * undergoing changes
 * 
 * @author nichole
 */
public class SunFinder {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public SunFinder() {
    }
    
    public SkyObject findSun(MSEREdges mserEdges) {
        
        ImageExt img = mserEdges.getClrImg();
        
        List<Region> greyscaleNegative = mserEdges.getOrigGsPtRegions().get(1);
        
        SunColors sunColors = new SunColors();
        
        //n=3, area=657, dens=0.0014 is not correct
        List<Set<PairInt>> listOfSets = new ArrayList<Set<PairInt>>();
        
        for (int rIdx = 0; rIdx < greyscaleNegative.size(); ++rIdx) {
            
            Region r = greyscaleNegative.get(rIdx);
           
            Set<PairInt> set1 = null;
                
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                PairInt p = new PairInt(x, y);
                int pixIdx = img.getInternalIndex(p);
               
                if (sunColors.isSunCenterColor(img, pixIdx)) {
                    if (set1 == null) {
                        set1 = new HashSet<PairInt>();
                    }
                    set1.add(p);
                }
            }
            if (set1 != null) {
                listOfSets.add(set1);
            }
        }
        
        int nList = listOfSets.size();
        List<EllipseHelper> ehs = new ArrayList<EllipseHelper>(nList);
        List<Set<PairInt>> listOfSets2 = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < listOfSets.size(); ++i) {
            Set<PairInt> points = listOfSets.get(i);
            EllipseHelper eh = createRegion(points);
            int n = points.size();
        
            double area = 2. * Math.PI * eh.getMajorTimesMinor();
            double dens = (double)n/area;
            if (dens > 0.1) {
                ehs.add(eh);
                listOfSets2.add(points);
            }
        }
        
        nList = listOfSets2.size();
        float[] densities = new float[nList];
        float[] ellipticities = new float[nList];
        int[] ns = new int[nList];
        int[] indexes = new int[nList];
        for (int i = 0; i < nList; ++i) {
            Set<PairInt> points = listOfSets.get(i);
            EllipseHelper eh = createRegion(points);
            int n = points.size();
        
            double area = 2. * Math.PI * eh.getMajorTimesMinor();
            double dens = (double)n/area;
            densities[i] = (float)dens;
            ellipticities[i] = (float)eh.getEccentricity();
            ns[i] = n;
            indexes[i] = i;
        }
        
        QuickSort.sortBy1stArg(densities, indexes);
        
        for (int i = (ns.length - 1); i > -1; --i) {
            
            Set<PairInt> points = listOfSets2.get(i);
            EllipseHelper eh = ehs.get(i);
            double area = 2. * Math.PI * eh.getMajorTimesMinor();
            double dens = (double)points.size()/area;
        
            //System.out.println(Arrays.toString(eh.getXYCenter()) 
            //    + " n=" + points.size()
            //    + " area=" + area + " dens=" + dens);
        
            SkyObject obj = new SkyObject();
            obj.points = points;
            obj.xyCenter = eh.getXYCenter();
            return obj;
        }
        
        //SkyObject obj = new SkyObject();
        //obj.points = listOfSets2.get(0);
        //obj.xyCenter = ehs.get(0).getXYCenter();
        //return obj;
        
        return null;        
    }

    private EllipseHelper createRegion(Set<PairInt> points) {

        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        int[] xyCen = ch.calculateRoundedXYCentroids(points);
    
        PairIntArray xy = Misc.convertWithoutOrder(points);
        
        EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], xy);
        
        return eh;
    }
  
}
