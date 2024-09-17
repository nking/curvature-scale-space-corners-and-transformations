package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.matrix.MatrixUtil;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.IntStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;

/**
 * class to read in SBA datasets

    caveat is that my BA code current needs a densely populated coordsI, so am taking a subsection of the datasets
    to use as test data.

    TODO: consider writing the data into BAL format so can compare the results of using ceres-solver on it.
 */
public class SBADataReader {

    public static String SBADIR = "sba";
    protected final static String sep = System.getProperty("file.separator");

    public static class BundleAdjustmentData {

        // TODO: include expected results

        /**
         * using cameras as a term for images, that is, solving each image for intr and extr even if it shares intr
         * with other observations.
         */
        int numCameras;

        /**
         * number of WCS points
         */
        int numPoints;

        int numObs;

        double[][] coordsI;
        double[][] coordsW;

        /** for each camera, add an intrinsic camera matrix.  [(3*nImages) X 3]
         */
        BlockMatrixIsometric intr;

        /**for each camera, a row of radial distortion coefficients [nImages X 2]
         */
        double[][] kRadials;

        //key = camera index, value = set of point indexes
        TIntObjectMap<TIntSet> imageFeaturesMap;

        /**the extrinsic camera parameter rotation vectors
        stacked along the 3 columns, that is the size is [nImages X 3] where
        nImages is coordsI[0].length/coordsW[0].length.  each array is size [1X3]
         */
        double[][] extrRotVecs;
        double[][] extrTrans;

        public BundleAdjustmentData copy() {

            BundleAdjustmentData data = new BundleAdjustmentData();

            if (this.imageFeaturesMap != null) {
                data.imageFeaturesMap = new TIntObjectHashMap<TIntSet>(imageFeaturesMap);
            }
            if (this.coordsI != null) {
                data.coordsI = MatrixUtil.copy(this.coordsI);
            }
            if (this.coordsW != null) {
                data.coordsW = MatrixUtil.copy(this.coordsW);
            }
            if (this.intr != null) {
                data.intr = intr.copy();
            }
            if (this.kRadials != null) {
                data.kRadials = MatrixUtil.copy(this.kRadials);
            }
            if (this.extrRotVecs != null) {
                data.extrRotVecs = MatrixUtil.copy(this.extrRotVecs);
            }
            if (this.extrTrans != null) {
                data.extrTrans = MatrixUtil.copy(this.extrTrans);
            }
            data.numCameras = numCameras;
            data.numPoints = numPoints;
            data.numObs = numObs;

            return data;
        }
    }

    /**
     * note that this assumes data.coordsI is densley written
     * @param data
     * @param fileName
     * @return
     * @throws IOException
     */
    public static void writeForCeresSolver(BundleAdjustmentData data, String fileName) throws IOException, NotConvergedException {

        BufferedWriter writer = getFileWriter(fileName);
        double[][] intr = new double[3][3];
        double[] v3;
        int cIdx, pIdx;
        int idx;
        int nPts = data.numPoints;

        writer.write(String.format("%d %d %d\n", data.numCameras, data.numPoints, data.numObs));

        boolean useR2R4 = true;

        double[][] x = new double[3][nPts];
        Arrays.fill(x[2], 1);

        int i = 0;
        //for numObs:  cam_idx  pt_idx xC yC
        for (cIdx = 0; cIdx < data.numCameras; ++cIdx) {
            data.intr.getBlock(intr, cIdx, 0);
            v3 = data.kRadials[cIdx];

            for (pIdx = 0; pIdx < nPts; ++pIdx) {
                //nCameras * nPoints.  so index is cIdxNew*nPoints + pIdx
                idx = cIdx * nPts + pIdx;

                x[0][pIdx] = data.coordsI[0][idx];
                x[1][pIdx] = data.coordsI[1][idx];
            }

            double[][] xC = Camera.pixelToCameraCoordinates(x, intr, v3, useR2R4);
            for (pIdx = 0; pIdx < nPts; ++pIdx) {
                for (int k = 0; k < 3; ++k) {
                    xC[k][pIdx] /= xC[2][pIdx];
                }
            }

            // for this camera write all points:  cam_idx  pt_idx xC yC
            for (pIdx = 0; pIdx < nPts; ++pIdx) {
                writer.write(String.format("%d %d %.6e %.6e\n", cIdx, pIdx, xC[0][pIdx], xC[1][pIdx]));
                ++i;
            }
        }
        assert(data.numObs == i);

        // then singly listed for these:
        //     first are the nCameras, 9 camera vars each (RotVec,trans,f,k1 and k2)
        //     then the remaining are the 3 WCS points for each numObs
        for (cIdx = 0; cIdx < data.numCameras; ++cIdx) {
            v3 = data.extrRotVecs[cIdx];
            for (int k = 0; k < 3; ++k) {
                writer.write(String.format("%.11e\n", v3[k]));
            }

            v3 = data.extrTrans[cIdx];
            for (int k = 0; k < 3; ++k) {
                writer.write(String.format("%.11e\n", v3[k]));
            }

            data.intr.getBlock(intr, cIdx, 0);
            double f =  (intr[0][0] + intr[1][1])/2.;
            writer.write(String.format("%.11e\n", f));

            v3 = data.kRadials[cIdx];
            for (int k = 0; k < 2; ++k) {
                writer.write(String.format("%.11e\n", v3[k]));
            }
        }

        //     then the remaining are the 3 WCS points for each numObs
        for (pIdx = 0; pIdx < nPts; ++pIdx) {
            for (int k = 0; k < 2; ++k) {
                writer.write(String.format("%.11e\n", data.coordsW[k][pIdx]));
            }
            // for BAL dataset, we need -Z
            writer.write(String.format("%.11e\n", -data.coordsW[2][pIdx]));
        }
        writer.flush();
        writer.close();

        int nLinesExpected = (data.numCameras*data.numPoints) + 9*(data.numCameras) + 3*data.numPoints + 1;
        int nLines = numberOfLinesW(fileName);
        assert(nLinesExpected == nLines);
    }

    protected static BufferedReader getFileReader(String fileName) throws IOException {
        String path = ResourceFinder.findTestResourcesDirectory() + sep + SBADIR
                + sep + fileName;

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        return new BufferedReader(new FileReader(new File(path)));
    }

    protected static BufferedWriter getFileWriter(String fileName) throws IOException {
        String path = ResourceFinder.findOutputTestDirectory() + sep + fileName;

        File f = new File(path);
        if (!f.exists()) {
            f.createNewFile();
        }

        return new BufferedWriter(new FileWriter(new File(path)));
    }

    public static class CamerasAndPoints {
        Set<Integer> cameraIndexes;
        Set<Integer> pointIndexes;
    }

    protected static int _numberOfLines(String filePath) throws IOException {

        File f = new File(filePath);

        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }

        try (Stream<String> stream = Files.lines(Path.of(f.getPath()), StandardCharsets.UTF_8)) {
            return (int) stream.count();
        }
    }

    protected static int numberOfLines(String fileName) throws IOException {

        String filePath = ResourceFinder.findTestResourcesDirectory() + sep + SBADIR
                + sep + fileName;

        return _numberOfLines(filePath);
    }
    protected static int numberOfLinesW(String fileName) throws IOException {

        String filePath = ResourceFinder.findOutputTestDirectory() + sep + fileName;

        return _numberOfLines(filePath);
    }

    /**
     * finds the at least imLimit number of cameras having at least ptLimit intersection in greedy manner, preferring larger
     * number of cameras over larger point intersection.
     * The default fixed limits are imLimit=5 and ptLimit=7.
     * To make more than one such intersection result CamerasAndPoints, one can overload this method to accept an ignore camera set,
     * and keep appending the results of invocation of this method to the ignore camera set until there are no new
     * blocks (while obviously modifying the method to check for whether a camera is in the ignore set).
     *
     * @param varFileName
     * @param ptsFileName
     * @return
     * @throws IOException
     */
    protected static CamerasAndPoints readSBAFileForIntersections(String varFileName, String ptsFileName) throws IOException {

        /*
        my code uses a densely populated coordsI so I'm reading the pts file in and finding an intersection of images
        and points where > 5 images have >= 7 points in common.

        a maximum number of intersection of subsets with min number of subsets and min number of points
        is an exponential algorithm.

        So instead of optimal solution, trying a greedy non-optimal solution.
        This approach is  (nCamera^3) and fails.
         */

        int ptLimit = 7;
        int imLimit = 5;

        int nImages = numberOfLines(varFileName) - 1;
        int nPoints = numberOfLines(ptsFileName) - 1;

        BufferedReader in = getFileReader(ptsFileName);

        int nRead = 0;

        String[] s;
        String line = in.readLine();
        assert(line.startsWith("#"));
        ++nRead;

        //key = camera index, value = set of point indexes
        TIntObjectMap<TIntSet> imagePointsMap = new TIntObjectHashMap<TIntSet>();
        TIntObjectMap<TIntSet> pointsImageMap = new TIntObjectHashMap<TIntSet>();

        int[] imagePointCounts = new int[nImages];
        //int[] pointImageCounts = new int[nPoints];

        int nP, i, cIdx;
        int pIdx = 0;

        TIntSet pSet, cSet;

        line = in.readLine();
        while (line != null) {

            s = line.split("\\s+");

            pSet = pointsImageMap.get(pIdx);
            if (pSet == null) {
                pSet = new TIntHashSet();
                pointsImageMap.put(pIdx, pSet);
            }

            nP = Integer.parseInt(s[3]);
            for (i = 4; i < 4 + nP*3; i+=3) {
                // extract camera numbers
                cIdx = Integer.parseInt(s[i]);

                pSet.add(cIdx);

                cSet = imagePointsMap.get(cIdx);
                if (cSet == null) {
                    cSet = new TIntHashSet();
                    imagePointsMap.put(cIdx, cSet);
                }
                cSet.add(pIdx);

                ++imagePointCounts[cIdx];
            }

            //pointImageCounts[pIdx] = nP;

            ++pIdx;
            line = in.readLine();
            ++nRead;
        }
        in.close();

        assert(imagePointsMap.size() == nImages);
        assert(pointsImageMap.size() == nPoints);

        // looking for the intersection of points in images to get at least 7 points common to at least 5 images.
        // need complexity to be less than exponential, so will make a rough greedy version that is not optimal.

        // 54 images, descending sort by point count.  max is 593 points
        int[] sortedImageIdxs
                = IntStream.range(0, imagePointCounts.length).boxed()
                .sorted((ii, jj)-> Integer.compare(imagePointCounts[jj], imagePointCounts[ii]))
                .mapToInt(ele->ele).toArray();

        // 5207 points.   max is 26 cameras
       /* int[] sortedPointIdxs
                = IntStream.range(0, pointImageCounts.length).boxed()
                .sorted((ii, jj)-> Integer.compare(pointImageCounts[jj], pointImageCounts[ii]))
                .mapToInt(ele->ele).toArray();*/

        TIntSet intersectingPoints = null;

        // O(nCamera^3) with factor of (avg nPoints/image)**3
        for (i = 0; i < sortedImageIdxs.length; ++i) {
            cIdx = sortedImageIdxs[i];
            TIntSet setI = imagePointsMap.get(cIdx);
            if (setI.size() < ptLimit) {
                break;
            }
            List<TIntSet> intersectionSets = new ArrayList<>();

            // O(nCamera^2) with factor of (avg nPoints/image)**2
            for (int j = i+1; j < sortedImageIdxs.length; ++j) {
                int cIdxJ = sortedImageIdxs[j];
                if (imagePointsMap.get(cIdxJ).size() < ptLimit) {
                    break;
                }
                TIntSet setJ = imagePointsMap.get(cIdxJ);
                TIntIterator iter = setJ.iterator();
                TIntSet intersection = new TIntHashSet();
                while (iter.hasNext()) {
                    pIdx = iter.next();
                    if (setI.contains(pIdx)) {
                        intersection.add(pIdx);
                    }
                }
                if (intersection.size() >= ptLimit) {
                    intersectionSets.add(intersection);
                }
            }
            if (intersectionSets.size() < imLimit) {
                continue;
            }
            // any point that appears < 7 times in intersectionSets gets removed.  any setI w/ < 7 points is removed
            trim(intersectionSets, ptLimit);

            // among intersectionSets, are there 5 with intersection size at least 7
            // sort the intersectionSets by size, desc

            if (intersectionSets.size() < imLimit) {
                continue;
            }

            // intersectionSets got sorted by descending order in trim
            // start with set0 and keep subtracting members not in next setI until we have 7
            for (int j = 0; j < intersectionSets.size(); ++j) {
                TIntSet setJ = intersectionSets.get(j);
                TIntSet intersection = new TIntHashSet(setJ);
                for (int k = j+1; k < intersectionSets.size(); ++k) {
                    TIntSet setK = intersectionSets.get(k);
                    TIntIterator iter = intersection.iterator();
                    TIntSet rm = new TIntHashSet();
                    while (iter.hasNext()) {
                        pIdx = iter.next();
                        if (!setK.contains(pIdx)) {
                            rm.add(pIdx);
                        }
                    }
                    int sz = intersection.size();
                    int nS = (k-j+1);
                    if (nS > imLimit && (sz - rm.size()) < ptLimit) {
                        // if we remove these points, the intersection is too small, so break now
                        intersectingPoints = intersection;
                        break;
                    }
                    intersection.removeAll(rm);
                    if (intersection.size() < ptLimit) {
                        // intersection too small, so break to start again with next j
                        break;
                    }
                }
                if (intersection.size() < ptLimit) {
                    // intersection too small, so break to start again with next j
                    break;
                }
                if (intersectingPoints != null) {
                    // we found a greedy intersection for starting with and including current j
                    break;
                }
            }

            if (intersectingPoints != null) {
                // we found a greedy intersection
                break;
            }
        }

        if (intersectingPoints == null) {
            // we did not find a greedy intersection between 5 or more images having the same 7 or more points
            throw new IllegalStateException("cannot use this dataset with my dense format coordsI");
        }

        // extract the images having all of the points in them
        Set<Integer> cameras = new HashSet<>();

        for (i = 0; i < sortedImageIdxs.length; ++i) {
            cIdx = sortedImageIdxs[i];
            TIntSet points = imagePointsMap.get(cIdx);
            if (points.containsAll(intersectingPoints)) {
                cameras.add(cIdx);
            }
        }

        CamerasAndPoints cap = new CamerasAndPoints();
        cap.cameraIndexes = cameras;
        cap.pointIndexes = new HashSet<>();
        TIntIterator iter = intersectingPoints.iterator();
        while (iter.hasNext()) {
            cap.pointIndexes.add(iter.next());
        }
        return cap;
    }

    private static void trim(List<TIntSet> sets, int limit) {
        TIntIntMap countMap = new TIntIntHashMap();
        for (TIntSet set : sets) {
            TIntIterator iter = set.iterator();
            while (iter.hasNext()) {
                int s = iter.next();
                if (countMap.containsKey(s)) {
                    countMap.put(s, countMap.get(s) + 1);
                } else {
                    countMap.put(s,  1);
                }
            }
        }
        TIntSet rm = new TIntHashSet();
        for (int key : countMap.keys()) {
            if (countMap.get(key) < limit) {
                rm.add(key);
            }
        }
        for (int i = 0; i < sets.size(); ++i) {
            TIntSet set = sets.get(i);
            set.removeAll(rm);
        }
        Collections.sort(sets, new Comparator<TIntSet>() {
            @Override
            public int compare(TIntSet o1, TIntSet o2) {
                return Integer.compare(o2.size(), o1.size());
            }
        });
        for (int i = sets.size() - 1; i >= 0; --i) {
            if (sets.get(i).size() < limit) {
                sets.remove(i);
            }
        }
    }


    protected static BundleAdjustmentData readSBAFileIntersection(String varFileName, String ptsFileName, int nCameraParams) throws IOException {

        if (nCameraParams != 17) {
            throw new UnsupportedOperationException("nCameraParams=17 is only file read supported for now");
        }

        // for my dense input coordsI, I would like to reduce the dataset to the intersections of images and points
        // and re-number the data for that purpose.
        // and write the output for this project and in BAL format to check my results with ceres-solver.
        //     to get best results from ceres-solver, I might need to write C++ code for the jacobians etc.
        CamerasAndPoints cap = readSBAFileForIntersections(varFileName, ptsFileName);

        // assign numbers to cap sets
        TIntIntMap cIdxToNewCIdx = new TIntIntHashMap();
        TIntIntMap pIdxToNewPIdx = new TIntIntHashMap();
        int c = 0;
        for (int cIdx : cap.cameraIndexes) {
            cIdxToNewCIdx.put(cIdx, c++);
        }
        c = 0;
        for (int pIdx : cap.pointIndexes) {
            pIdxToNewPIdx.put(pIdx, c++);
        }

        int nCameras = cIdxToNewCIdx.size();
        int nPoints = pIdxToNewPIdx.size();

        BundleAdjustmentData data = new BundleAdjustmentData();
        data.imageFeaturesMap = new TIntObjectHashMap<TIntSet>();
        data.coordsI = new double[3][nCameras * nPoints];
        Arrays.fill(data.coordsI[2], 1);
        data.coordsW = new double[3][nPoints];
        //[(3*nImages) X 3]
        data.intr = new BlockMatrixIsometric(MatrixUtil.zeros(nCameras * 3, 3), 3, 3);
        data.kRadials = new double[nCameras][];
        data.extrRotVecs = new double[nCameras][];
        data.extrTrans = new double[nCameras][];
        data.numCameras = nCameras;
        data.numPoints = nPoints;
        data.numObs = nCameras * nPoints;// this is calculated because we are filtering coordsI to be densely populated

        readVarFile(nCameraParams, data, varFileName, cap, cIdxToNewCIdx);

        readPtsFile(data, ptsFileName, cap, cIdxToNewCIdx, pIdxToNewPIdx);

        return data;
    }

    private static void readPtsFile(BundleAdjustmentData data, String ptsFileName,
                                    CamerasAndPoints cap, TIntIntMap cIdxToNewCIdx, TIntIntMap pIdxToNewPIdx) throws IOException {

        BufferedReader in = getFileReader(ptsFileName);

        String[] s;
        String line = in.readLine();
        assert(line.startsWith("#"));

        s = line.split("\\s+");

        int pIdx = 0;
        int pIdxNew;
        int i;
        int cIdx, cIdxNew;

        int nPts = data.numPoints;
        int nCameras = data.numCameras;
        int nObs = data.numObs;

        int idx;

        line = in.readLine();
        while (line != null) {

            if (!pIdxToNewPIdx.containsKey(pIdx)) {
                line = in.readLine();
                ++pIdx;
                continue;
            }
            pIdxNew = pIdxToNewPIdx.get(pIdx);

            s = line.split("\\s+");

            //# X Y Z  nframes  frame0 x0 y0  frame1 x1 y1 ...
            data.coordsW[0][pIdxNew] = Double.parseDouble(s[0]);
            data.coordsW[1][pIdxNew] = Double.parseDouble(s[1]);
            data.coordsW[2][pIdxNew] = Double.parseDouble(s[2]);

            // we filter the projected points to include only those in cIdxToNewPix
            for (i = 4; i < s.length; i+=3) {
                cIdx = Integer.parseInt(s[i]);
                if (!cIdxToNewCIdx.containsKey(cIdx)) continue;

                cIdxNew = cIdxToNewCIdx.get(cIdx);

                idx = cIdxNew * nPts + pIdxNew;
                //nCameras * nPoints.  so index is cIdxNew*nPoints + pIdx
                data.coordsI[0][idx] = Double.parseDouble(s[i+1]);
                data.coordsI[1][idx] = Double.parseDouble(s[i+2]);

                if (!data.imageFeaturesMap.containsKey(cIdxNew)) {
                    data.imageFeaturesMap.put(cIdxNew, new TIntHashSet());
                }
                data.imageFeaturesMap.get(cIdxNew).add(pIdxNew);
            }

            line = in.readLine();
            ++pIdx;
        }
        in.close();
    }

    private static void readVarFile(int nCameraParams, BundleAdjustmentData data, String varFileName, CamerasAndPoints cap,
                                    TIntIntMap cIdxToNewCIdx) throws IOException {

        BufferedReader in = getFileReader(varFileName);

        String[] s;
        String line = in.readLine();
        assert(line.startsWith("#"));

        s = line.split("\\s+");

        double[][] intr;

        int cIdx = 0;
        int cIdxNew;
        double fx;
        double[] q = new double[4];

        line = in.readLine();
        while (line != null) {

            if (!cIdxToNewCIdx.containsKey(cIdx)) {
                line = in.readLine();
                ++cIdx;
                continue;
            }
            cIdxNew = cIdxToNewCIdx.get(cIdx);

            s = line.split("\\s+");
            assert(s.length == nCameraParams);

            if (nCameraParams==17) {
                //# fu, u0, v0, ar, s   kc(1:5)   quaternion   translation
                //   0   1   2   3   4    56789   10 11 12 13  14 15 16
                fx = Double.parseDouble(s[0]);
                intr = new double[][] {
                        {fx, Double.parseDouble(s[4]), Double.parseDouble(s[1])},
                        {0, fx*Double.parseDouble(s[3]), Double.parseDouble(s[2])},
                        {0, 0, 1}
                };
                data.intr.setBlock(intr, cIdxNew, 0);

                data.kRadials[cIdxNew] = new double[]{Double.parseDouble(s[5]), Double.parseDouble(s[6])};

                // they provide hamilton quaternion which is scalar first.
                // convert to barfoot by placing scalar last, then can use rotationvector method
                q[0] = Double.parseDouble(s[11]);
                q[1] = Double.parseDouble(s[12]);
                q[2] = Double.parseDouble(s[13]);
                q[3] = Double.parseDouble(s[10]);
                data.extrRotVecs[cIdxNew] = Rotation._extractRotationVectorFromQuaternionBarfoot(q);

                data.extrTrans[cIdxNew] = new double[]{Double.parseDouble(s[14]), Double.parseDouble(s[15]), Double.parseDouble(s[16])};
            }
            //TODO" write clauses for the other parameter files as needed

            line = in.readLine();
            ++cIdx;
        }
        in.close();

    }

}
