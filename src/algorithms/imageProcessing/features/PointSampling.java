package algorithms.imageProcessing.features;

import algorithms.compGeometry.PointPartitioner;
import algorithms.compGeometry.PointPartitioner.Bounds;
import algorithms.imageProcessing.DistanceTransform;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class PointSampling {
    
    /**
     * NOT READY FOR USE YET.  The results are not uniform yet - there are
     * irregularities for drawing the smallest indexes (less than 2*k) and
     * for the largest indexes (larger than nPoints - 2*k).
     * 
     * @param pv
     * @param sr
     * @param alreadySelected
     * @param outputPoints 
     */
    public void choose7RandomPoints(PointValueDistr pv, Random sr, 
        Set<BigInteger> alreadySelected, List<PairInt> outputPoints) {
        
        chooseKRandomPoints(pv, sr, alreadySelected, outputPoints, 7);
    }
    
    // keeping this private until have tested range of values for k.
    // needs to be correct for k=7
    private void chooseKRandomPoints(PointValueDistr pv, Random sr, 
        Set<BigInteger> alreadySelected, List<PairInt> outputPoints, int k) {
        
        // NOTE: ideally, want a way to forward index to gosper's hack
        // without calculating every one of the results,
        // and then choose randomly between 0 and the number of possible
        // combinations, then use that randomly chosen number to pick the 
        // combination (that is the bit vector) from the forward indexed 
        // gosper's hack
        //
        // Not having that yet, have implemented instead a possibly biased
        // (haven't tested yet) k-bit random number selector.
        
        /*
        choosing 7 indexes from an array with a single random number would
        be as follows:
        
        example:
            n=10 numbers in pv, k=7
        
        min bit vector = "0111 111"     = 127
        max bit vector = "11 1111 1000" = (127 << 3) = 1016
        
        randomly selecting between numbers 127 and 1016.
        if the number of set bits isn't k, find nearest next lowest value w/
            the required number of set bits, etc.
        
        ==> That would be if wanting to treat each point as equally valued. <==
        
        This method instead wants to represent each value that many number of
        times, so the PointValueDistr was created.
        In PointValueDistr is a field called maxValue which holds the cumulative 
        number of values without storing each one in an array.  
        
        For example, if the original values are [2, 3, 4], the probability 
        distribution based upon value would be [2, 2, 3, 3, 3, 4, 4, 4, 4].
        The cumulative indexes for [2, 3, 4] would then be [2, 5, 9].
        Choosing from a range between 0 and maxValue of 9 then selects in a 
        uniform manner from the value based probability distribution.
        
        ==> The next adjustment to that selection model is during the random
        selection of the number.  Have to transform from its cumulative bit 
        vector to the pv.points equivalent bit vector in order to make sure
        the result is uniquely chosen pv.points indexes.
        */
        
        // a bitvector in the reference frame of pv.points
        BigInteger randomlyChosen = randomlyChooseKBitNumber(pv, sr, k);
        
        while (alreadySelected.contains(randomlyChosen)) {
            randomlyChosen = randomlyChooseKBitNumber(pv, sr, k);
        }
        
        alreadySelected.add(randomlyChosen);
        
        // read off the set bits in randomlyChosen to get the pairints from pv
        
        outputPoints.clear();
        
        int bitLength = randomlyChosen.bitLength();
        for (int i = (bitLength - 1); i > -1; --i) {
            
            if (randomlyChosen.testBit(i)) {
                               
                PairInt p = pv.getPoints()[i];
                
                outputPoints.add(p);
            }
            if (outputPoints.size() == k) {
                break;
            }
        }
    }
    
    /**
     * randomly choose a k bit number that is a bit vector that has set bits
     * representing the pv.points indexes selected.
     * NOT READY FOR USE YET.  The results are not uniform yet - there are
     * irregularities for drawing the smallest indexes (less than 2*k) and
     * for the largest indexes (larger than nPoints - 2*k).
     * 
     * NOTE: keeping this nearly private until have tested range of values for k.
       needs to be correct for k=7.
     
     * @param maxValue
     * @param sr
     * @param k
     * @return 
     */
    BigInteger randomlyChooseKBitNumber(PointValueDistr pv, Random sr, int k) {
        
        if (k != 7) {
            throw new IllegalArgumentException(
            "algorithm is currently only tested for k=7");
        }
        
        int nPoints = pv.getPoints().length;
        
        int minBitsValue = (1 << k) - 1;
        
        /*
        wanting to choose 7 values from range 0 to nCumulativeValues using one 
        random number where nCumulativeValues is pv.maxValue.
        
        bit vector with each bit holding a meaning of chosen ("1") or not ("0").
        
        the minimum usable bit vector is 7 bits set, = (1 << 7) - 1 = 127
        the maximum usable bit vector is the top 7 bits set of
             nCumulativeValues number of bits.
             = 127 << (nCumulativeValues - 7)
        
        select random numbers between minimum and maximum usable bit vector.
        
        then convert that to the bit vector representing maxOriginalIndex
        where maxOriginalIndex is pv.points.size().
        
        and then if 7 bits exactly are not set, flip bits to the closest 
        lower number to arrive at 7 bits set, etc.
        */
                
        BigInteger maxCumulativeIndex = pv.getMaxValue();
                
        if (maxCumulativeIndex.bitLength() > 31) {
            throw new ArithmeticException("result does not fit within an int");
        }
        int vMax = maxCumulativeIndex.intValue();
        
        // --- randomly choose a number between 0 and 2^nBits. ---
        
        // the random algorithms and BigInteger algorithm appear to be biased
        // towards very large numbers when nBits is high, so
        // working around that by randomly selecting nBits, then using
        // the random from BigInteger or Random.
        
        int nBits = sr.nextInt(vMax - k) + k;
        BigInteger randomlyChosen = new BigInteger(nBits, sr);
        
        int bitLength = randomlyChosen.bitLength();
        
        //System.out.println("before: " + randomlyChosen.toString(2));
        
        // --- create a bitstring in the reference frame of pv.points ---
        BigInteger rc0 = BigInteger.ZERO;
        for (int i = 0; i < bitLength; ++i) {
            
            if (randomlyChosen.testBit(i)) {
                
                int oIdx = pv.getPointsIndexForCumulativeValue(
                    new BigInteger(MiscMath.writeToBigEndianBytes(i)));
          
                rc0 = rc0.setBit(oIdx);
            }
        }
        
        randomlyChosen = rc0;
                
        bitLength = randomlyChosen.bitLength();
        
        int bitCount = randomlyChosen.bitCount();
        
        if (bitCount > k) {
            
            //keep the top k bits
                        
            BigInteger sum = BigInteger.ZERO;
            int c = 0;
            for (int j = (bitLength - 1); j > -1; --j) {
                if (randomlyChosen.testBit(j)) {
                    BigInteger c2 = BigInteger.ONE;
                    c2 = c2.shiftLeft(j);
                    sum = sum.add(c2);
                    c++;
                    if (c == k) {
                        break;
                    }
                }
            }
            
            randomlyChosen = sum;
                        
        } else if (bitCount < k) {

            if (bitLength <= k) {
                
                randomlyChosen = new BigInteger(MiscMath.writeToBigEndianBytes(minBitsValue));
                
            } else {
                
                // plus one because unsetting a '1' and need to replace it, then 
                // setting '0's below it
                int nBitsToSet = k + 1 - bitCount;

                // flip the low 0's to 1's and flip the next left 1 to 0
                
                bitLength = randomlyChosen.bitLength();
                int c = 0;
                int last = 0;
                for (int j = 0; j < bitLength; ++j) {
                    if (!randomlyChosen.testBit(j)) {
                        randomlyChosen = randomlyChosen.setBit(j);
                        c++;
                        last = j;
                    }
                    if (c == nBitsToSet) {
                        break;
                    }
                }
                
                for (int j = (last + 1); j < bitLength; ++j) {
                    if (randomlyChosen.testBit(j)) {
                        if (j != (bitLength - 1)) {
                            randomlyChosen = randomlyChosen.clearBit(j);
                        } else {
                            randomlyChosen = randomlyChosen.clearBit(last);
                        }
                        break;
                    }
                }                
            }
        }

        //System.out.println("after:  " + randomlyChosen.toString(2) + "\n");
        
        return randomlyChosen;
    }
    
    /**
     * randomly choose a k bit number that is a bit vector that has set bits
     * representing the pv.points indexes selected.
     * NOT READY FOR USE YET.  The results are not uniform yet - there are
     * irregularities for drawing the smallest indexes (less than 2*k) and
     * for the largest indexes (larger than nPoints - 2*k).
     * 
     * NOTE: keeping this nearly private until have tested range of values for k.
       needs to be correct for k=7.
     
     * @param nIndexes
     * @param sr
     * @param k
     * @return 
     */
    BigInteger randomlyChooseKBitNumber(int nIndexes, Random sr, int k) {
        
        if (k != 7) {
            throw new IllegalArgumentException(
            "algorithm is currently only tested for k=7");
        }
        
        int minBitsValue = (1 << k) - 1;
        
        /*
        wanting to choose 7 values from range 0 to nIndexes using one 
        random number.
        
        bit vector with each bit holding a meaning of chosen ("1") or not ("0").
        
        the minimum usable bit vector is 7 bits set, = (1 << 7) - 1 = 127
        the maximum usable bit vector is 1 << nIndexes
        
        select random numbers between minimum and maximum usable bit vector.
        
        and then if 7 bits exactly are not set, flip bits to the closest 
        lower number to arrive at 7 bits set, etc.
        */
        
        // --- randomly choose a number between 0 and 2^nBits. ---
        
        // the random algorithms and BigInteger algorithm appear to be biased
        // towards very large numbers when nBits is high, so
        // working around that by randomly selecting nBits, then using
        // the random from BigInteger or Random.
        
        int nBits = sr.nextInt(nIndexes - k) + k;
        BigInteger randomlyChosen = new BigInteger(nBits, sr);
       
        //System.out.println("before: " + randomlyChosen.toString(2));
          
        int bitLength = randomlyChosen.bitLength();
        
        int bitCount = randomlyChosen.bitCount();
        
        if (bitCount > k) {
            
            //keep the top k bits
            
            BigInteger sum = BigInteger.ZERO;
            int c = 0;
            for (int j = (bitLength - 1); j > -1; --j) {
                if (randomlyChosen.testBit(j)) {
                    BigInteger c2 = BigInteger.ONE;
                    c2 = c2.shiftLeft(j);
                    sum = sum.add(c2);
                    c++;
                    if (c == k) {
                        break;
                    }
                }
            }

            randomlyChosen = sum;
                        
        } else if (bitCount < k) {

            if (bitLength <= k) {
                
                randomlyChosen = new BigInteger(MiscMath.writeToBigEndianBytes(minBitsValue));
                
            } else {
                
                // plus one because unsetting a '1' and need to replace it, then 
                // setting '0's below it
                int nBitsToSet = k + 1 - bitCount;

                // flip the low 0's to 1's and flip the next left 1 to 0
                
                bitLength = randomlyChosen.bitLength();
                int c = 0;
                int last = 0;
                for (int j = 0; j < bitLength; ++j) {
                    if (!randomlyChosen.testBit(j)) {
                        randomlyChosen = randomlyChosen.setBit(j);
                        c++;
                        last = j;
                    }
                    if (c == nBitsToSet) {
                        break;
                    }
                }
                
                for (int j = (last + 1); j < bitLength; ++j) {
                    if (randomlyChosen.testBit(j)) {
                        if (j != (bitLength - 1)) {
                            randomlyChosen = randomlyChosen.clearBit(j);
                        } else {
                            randomlyChosen = randomlyChosen.clearBit(last);
                        }
                        break;
                    }
                }                
            }
        }

        //System.out.println("after:  " + randomlyChosen.toString(2) + "\n");
        
        return randomlyChosen;
    }
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension of roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner, biased towards points closer to the 
     * seed centers.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * @return  
     */
    public PointValueDistr createSpatialInvDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions) {
        
        return createSpatialDistBasedValues(points, numCellsPerDimensions, true);
    }
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension of roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner, biased towards points further from 
     * the seed centers.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * @return  
     */
    public PointValueDistr createSpatialDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions) {
        
        return createSpatialDistBasedValues(points, numCellsPerDimensions, false);
    }
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension or roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * #param inverse if true, uses maxDist - distance for distances to weight
     * for points closer to the seeds.
     * @return  
     */
    protected PointValueDistr createSpatialDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions, boolean inverse) {
        
        if (numCellsPerDimensions < 1) {
            throw new IllegalArgumentException(
                "numCellsPerDimensions must be > 0");
        }
        
        PointPartitioner pp = new PointPartitioner();
        
        List<Bounds> bounds = pp.findCells(numCellsPerDimensions, points);
        
        int[] maxXY = findMaxXAndY(bounds);
        
        List<PairInt> seeds = centers(bounds);
        
        int w = maxXY[0] + 1;
        int h = maxXY[1] + 1;
        
        // when n_points * log_2(n_points) is less than n_pixels, prefer
        // vornoi instead of distance transform to find cells then calculate 
        // distance of points to seeds.  TODO: impl vornoi fortunes.
        
        Map<PairInt, Integer> pointDistMap = useDistanceTransform(points, w, h, 
            seeds);
        
        int maxDist = Integer.MIN_VALUE;
        if (inverse) {
            for (Entry<PairInt, Integer> entry : pointDistMap.entrySet()) {
                int dist = entry.getValue().intValue();
                if (dist > maxDist) {
                    maxDist = dist;
                }
            }
        }
        
        int n = pointDistMap.size();
        
        BigInteger[] cumulativeCount = new BigInteger[pointDistMap.size()];
        PairInt[] pointsArray = new PairInt[pointDistMap.size()];
        
        int count = 0;
        BigInteger prev = null;
        for (Entry<PairInt, Integer> entry : pointDistMap.entrySet()) {
            
            int value = entry.getValue().intValue();
            if (inverse) {
                value = maxDist + 1 - value;
            }
            
            pointsArray[count] = entry.getKey();
            if (prev == null) {
                prev = BigInteger.ZERO;
            }
            
            BigInteger current = new BigInteger(MiscMath.writeToBigEndianBytes(value));
            
            cumulativeCount[count] = prev.add(current);
            
            prev = cumulativeCount[count];
            
            count++;
        }
        
        PointValueDistr pv = new PointValueDistr(prev, pointsArray, cumulativeCount);
        
        return pv;
    }
    
    protected List<PairInt> centers(List<Bounds> bounds) {
        
        List<PairInt> centers = new ArrayList<PairInt>();
        
        for (Bounds b : bounds) {
            
            float avgX = (b.lowerLeft.getX() + b.lowerRight.getX() 
                + b.upperLeft.getX() + b.upperRight.getX())/4.f;
            
            float avgY = (b.lowerLeft.getY() + b.lowerRight.getY() 
                + b.upperLeft.getY() + b.upperRight.getY())/4.f;
            
            PairInt p = new PairInt(Math.round(avgX), Math.round(avgY));
            
            centers.add(p);
        }
        
        return centers;
    }

    private int[] findMaxXAndY(List<Bounds> bounds) {
        
        int maxX = Integer.MIN_VALUE;
        
        int maxY = Integer.MIN_VALUE;
        
        for (Bounds b : bounds) {
            if (b.lowerRight.getX() > maxX) {
                maxX = b.lowerRight.getX();
            }
            if (b.upperRight.getX() > maxX) {
                maxX = b.upperRight.getX();
            }
            if (b.upperLeft.getY() > maxY) {
                maxY = b.upperLeft.getY();
            }
            if (b.upperRight.getY() > maxY) {
                maxY = b.upperRight.getY();
            }
        }
        
        return new int[]{maxX, maxY};
    }

    /**
     * use a  
       density based clustering algorithm
       http://nking.github.io/two-point-correlation/
       which has an MIT license
      ---- begin nking copyright ----
      The MIT License (MIT)
      Copyright (c) 2013-* Nichole King
      http://nking.github.io/two-point-correlation/

        Permission is hereby granted, free of charge, to any person obtaining 
        a copy of this software and associated documentation files 
        (the "Software"), to deal in the Software without restriction, 
        including without limitation the rights to use, copy, modify, merge, 
        publish, distribute, sublicense, and/or sell copies of the Software, 
        and to permit persons to whom the Software is furnished to do so, 
        subject to the following conditions:

        The above copyright notice and this permission notice shall be included 
        in all copies or substantial portions of the Software.
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
        OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
        IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
        CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
        SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
     ---- end nking copyright ---- 
     * @param points
     * @param w
     * @param h
     * @param seeds
     * @return 
     */
    private Map<PairInt, Integer> useDistanceTransform(Set<PairInt> points, 
        int w, int h, List<PairInt> seeds) {
             
        DistanceTransform dt = new DistanceTransform();
        
        PixelHelper ph = new PixelHelper();
        
        TLongSet seedPixs = new TLongHashSet();
        
        for (PairInt p : seeds) {
            long pixIdx = ph.toPixelIndex(p, w);
            seedPixs.add(pixIdx);
        }
        
        int[][] distances = dt.applyMeijsterEtAl(seedPixs, w, h);
        
        Map<PairInt, Integer> distMap = new HashMap<PairInt, Integer>();
        for (PairInt p : points) {
            int distSq = distances[p.getX()][p.getY()];
            Integer dist = Integer.valueOf((int)Math.round(Math.sqrt(distSq)));
            distMap.put(p, dist);
        }
        
        return distMap;
    }    
    
}
