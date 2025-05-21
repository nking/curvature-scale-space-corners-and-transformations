package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

/**
 * class to read in BAL datasets and format it into input for
 * BundledAdjustment to use for testing.
 *
 * Have paused editing this and will not commit it.
 *
 * To format their data to use as tests for by BA would require some assumptions because their
 * datasets are missing intrinsic camera information.
 * I wrote some notes below about how to roughly estimate the principal points,
 * but this is not a good use of time right now...
 *
 */
public class BALDataReader {

    public static String BALDIR = "bal";
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
    }

    public static BundleAdjustmentData readLadyBug() {
        throw new UnsupportedOperationException("not yet implemented");
        /*
        BundleAdjustmentData data = new BundleAdjustmentData();
        readVarKD("54camsvarKD.txt", data);
        readPoints("54pts.txt", data);
        return data;*/
    }

    public static String writeForCeresSolver() {
        /*
        needs to follow data format at
        https://grail.cs.washington.edu/projects/bal/

        Each problem is provided as a bzip2 compressed text file in the following format.
        <num_cameras> <num_points> <num_observations>
        <camera_index_1> <point_index_1> <x_1> <y_1>
        ...
        <camera_index_num_observations> <point_index_num_observations> <x_num_observations> <y_num_observations>
        <camera_1>
        ...
        <camera_num_cameras>
        <point_1>
        ...
        <point_num_points>

        Where, there camera and point indices start from 0.
        Each camera is a set of 9 parameters - R,t,f,k1 and k2.
        The rotation R is specified as a Rodrigues' vector.
         */
        /* **** 1st line read 3 params:
        FscanfOrDie(fptr, "%d", &num_cameras_);
        FscanfOrDie(fptr, "%d", &num_points_);
        FscanfOrDie(fptr, "%d", &num_observations_);

        point_index_ = new int[num_observations_];
        camera_index_ = new int[num_observations_];
        observations_ = new double[2 * num_observations_];

        num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
        parameters_ = new double[num_parameters_];

        *** next num_observations_ lines
        read  cam_idx  pt_idx x y

        *** next num_parameters_ lines
        read single float for params
        */
        throw new UnsupportedOperationException("not yet implemented");
    }

    protected static BufferedReader getFileReader(String fileName) throws IOException {
        String path = ResourceFinder.findTestResourcesDirectory() + sep + BALDIR
                + sep + fileName;

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        return new BufferedReader(new FileReader(new File(path)));
    }

    public static class CamerasAndPoints {
        Set<Integer> cameraIndexes;
        Set<Integer> pointIndexes;
    }

    protected static CamerasAndPoints readBALFileForIntersections(String fileName) throws IOException {

        BufferedReader in = getFileReader(fileName);

        int nRead = 0;

        String[] s;
        String line = in.readLine();
        ++nRead;

        s = line.split("\\s+");

        final int numCameras = Integer.parseInt(s[0]);
        // number of WCS points:
        final int numPoints = Integer.parseInt(s[1]);
        final int numObs = Integer.parseInt(s[2]);

        //key = camera index, value = set of point indexes
        TIntObjectMap<TIntSet> imagePointsMap = new TIntObjectHashMap<TIntSet>();
        TIntObjectMap<TIntSet> pointsImageMap = new TIntObjectHashMap<TIntSet>();

        int[] imagePointCounts = new int[numCameras];
        int[] pointImageCounts = new int[numPoints];

        int cIdx, pIdx, i;
        int nObsCount = 0;

        // 1 point should be in 5 or so views
        // store camera, pts
        // pts, camera ids

        // for my dense input coordsI, I would like to reduce the dataset to the intersections of images and points
        // and re-number the data for that purpose.
        // and write the output for this project and in BAL format to check my results with ceres-solver.
        //     to get best results from ceres-solver, I might need to write C++ code for the jacobians etc.

        line = in.readLine();
        while (line != null) {

            s = line.split("\\s+");

            if (nObsCount < numObs) {
                //read:  cam_idx  pt_idx x y
                cIdx = Integer.parseInt(s[0]);
                pIdx = Integer.parseInt(s[1]);

                TIntSet set = imagePointsMap.get(cIdx);
                if (set == null) {
                    set = new TIntHashSet();
                    imagePointsMap.put(cIdx, set);
                }
                set.add(pIdx);

                set = pointsImageMap.get(pIdx);
                if (set == null) {
                    set = new TIntHashSet();
                    pointsImageMap.put(pIdx, set);
                }
                set.add(cIdx);

                ++imagePointCounts[cIdx];
                ++pointImageCounts[pIdx];

                ++nObsCount;
            } else {
                break;
            }
            line = in.readLine();
            ++nRead;
        }
        in.close();

        // looking for the intersection of points in images to get at least 7 points common to at least 5 images.
        // need complexity to be less than exponential, so will make a rough greedy version that is not optimal.

        // 49 images, descending sort by point count
        int[] sortedImageIdxs
                = IntStream.range(0, imagePointCounts.length).boxed()
                .sorted((ii, jj)-> Integer.compare(imagePointCounts[jj], imagePointCounts[ii]))
                .mapToInt(ele->ele).toArray();

        // 7776 points
        /*int[] sortedPointIdxs
                = IntStream.range(0, pointImageCounts.length).boxed()
                .sorted((ii, jj)-> Integer.compare(pointImageCounts[jj], pointImageCounts[ii]))
                .mapToInt(ele->ele).toArray();*/

        TIntSet intersectingPoints = null;

        // O(nCamera^3) with factor of (avg nPoints/image)**3
        for (i = 0; i < sortedImageIdxs.length; ++i) {
            cIdx = sortedImageIdxs[i];
            TIntSet set = imagePointsMap.get(cIdx);
            List<TIntSet> sets = new ArrayList<>();

            // O(nCamera^2) with factor of (avg nPoints/image)**2

            for (int j = i+1; j < sortedImageIdxs.length; ++j) {
                int cIdxJ = sortedImageIdxs[j];
                if (imagePointsMap.get(cIdxJ).size() < 7) {
                    break;
                }
                TIntSet setJ = imagePointsMap.get(cIdxJ);
                TIntIterator iter = setJ.iterator();
                TIntSet intersection = new TIntHashSet();
                while (iter.hasNext()) {
                    pIdx = iter.next();
                    if (set.contains(pIdx)) {
                        intersection.add(pIdx);
                    }
                }
                sets.add(intersection);
            }
            if (sets.size() < 5) {
                continue;
            }
            // any point that appears < 7 times in sets gets removed
            trim(sets, 7);

            // among sets, are there 5 with intersection size at least 7
            // sort the sets by size, desc

            if (sets.size() < 5) {
                continue;
            }

            // sets got sorted by descending order in trim
            // start with set0 and keep subtracting members not in next set until we have 7
            for (int j = 0; j < sets.size(); ++j) {
                TIntSet setJ = sets.get(j);
                TIntSet intersection = new TIntHashSet(setJ);
                for (int k = j+1; k < sets.size(); ++k) {
                    TIntSet setK = sets.get(k);
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
                    if (nS >= 5 && (sz - rm.size()) < 7) {
                        intersectingPoints = intersection;
                        break;
                    }
                    intersection.removeAll(rm);
                    if (sz < 7) {
                        break;
                    }
                }
                if (intersection.size() < 7) {
                    break;
                }
                if (intersectingPoints != null) {
                    break;
                }
            }

            if (intersectingPoints != null) {
                break;
            }
        }

        if (intersectingPoints == null) {
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


    protected static BundleAdjustmentData readBALFileIntersection(String fileName) throws IOException {

        // for my dense input coordsI, I would like to reduce the dataset to the intersections of images and points
        // and re-number the data for that purpose.
        // and write the output for this project and in BAL format to check my results with ceres-solver.
        //     to get best results from ceres-solver, I might need to write C++ code for the jacobians etc.
        CamerasAndPoints cap = readBALFileForIntersections(fileName);

        // we don't have the principal point, image center for the LadyBug dataset.
        // browsing the BAL dataset info, this is the same camera, a FLIR LadyBug camera, but no information about the
        // model or number of pixels.
        // we can roughly estimate the sensor size by the range of x and y values which are in camera coordinate frame.
        // but we need the intrinsic camera.
        // for every cameras observations, we could transform the camera coordinates to image coordinates,
        // and then find the maximum imX, imY in all of this dataset, and then estimate that the principal
        // point coordinates is half of thos image max values.
        //      P  =  R * X + t       (conversion from world to camera coordinates)
        //            p  = -P / P.z         (perspective division)
        //            p' =  f * r(p) * p    (conversion to pixel coordinates)

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

        BufferedReader in = getFileReader(fileName);

        int nRead = 0;

        String[] s;
        String line = in.readLine();
        ++nRead;

        s = line.split("\\s+");

        final int numCameras = Integer.parseInt(s[0]);
        // number of WCS points:
        final int numPoints = Integer.parseInt(s[1]);
        final int numObs = Integer.parseInt(s[2]);
        final int numParameters = 9 * numCameras + 3 * numPoints;

        BundleAdjustmentData data = new BundleAdjustmentData();
        data.imageFeaturesMap = new TIntObjectHashMap<TIntSet>();
        data.coordsI = new double[3][cIdxToNewCIdx.size()*pIdxToNewPIdx.size()];
        data.coordsW = new double[3][pIdxToNewPIdx.size()];
        //data.intr = no principal point constraints see notes above

        int cIdx, pIdx, i;
        double xC, yC;
        int nObsCount = 0;
        int nParamCount = 0;

        /*
        They use a negative sign convention for z-axis:

            P  =  R * X + t       (conversion from world to camera coordinates)
            p  = -P / P.z         (perspective division)
            p' =  f * r(p) * p    (conversion to pixel coordinates)
            where P.z is the third (z) coordinate of P.
            In the last equation, r(p) is a function that computes a scaling factor to undo the radial distortion:
            r(p) = 1.0 + k1 * ||p||^2 + k2 * ||p||^4.
         */

        double[] rVec = new double[3];
        double[] trans = new double[3];
        double[] wcs = new double[3];
        double v;
        double f = 0;
        double[] kRadial24 = new double[2];
        int nC = 0;
        int nP = 0;

        line = in.readLine();
        while (line != null) {

            s = line.split("\\s+");

            if (nObsCount < numObs) {
                //read:  cam_idx  pt_idx x y
                cIdx = Integer.parseInt(s[0]);
                pIdx = Integer.parseInt(s[1]);

                if (!cap.cameraIndexes.contains(cIdx) || !cap.pointIndexes.contains(pIdx)) {
                    line = in.readLine();
                    ++nRead;
                    continue;
                }

                xC = Double.parseDouble(s[2]);
                yC = Double.parseDouble(s[3]);

                TIntSet set = data.imageFeaturesMap.get(cIdx);
                if (set == null) {
                    set = new TIntHashSet();
                    data.imageFeaturesMap.put(cIdx, set);
                }
                set.add(pIdx);

                ++nObsCount;

            } else if (nParamCount < numParameters) {
                // singly listed for these:
                // first are the nCameras, 9 camera vars each
                // then the remaining are the 3 WCS points for each numObs
                // read each parameter (there is only 1 per line, but there are 9 per camera frame
                if (nC < numCameras) {
                    for (i = 0; i < 9; ++i) {
                        //R,t,f,k1 and k2
                        v = Double.parseDouble(s[0]);
                        switch (i) {
                            case 0:
                                rVec[0] = v;
                                break;
                            case 1:
                                rVec[1] = v;
                                break;
                            case 2:
                                rVec[2] = v;
                                break;
                            case 3:
                                trans[0] = v;
                                break;
                            case 4:
                                trans[1] = v;
                                break;
                            case 5:
                                trans[2] = v;
                                break;
                            case 6:
                                f = v;
                                break;
                            case 7:
                                kRadial24[0] = v;
                                break;
                            default:
                                kRadial24[1] = v;
                                break;
                        }
                    }
                    // store into data using new idx mappings
                    nParamCount += 9;
                    ++nC;
                } else {
                    // WCS 3D points
                    for (i = 0; i < 3; ++i) {
                        v = Double.parseDouble(s[0]);
                        switch (i) {
                            case 0:
                                wcs[0] = v;
                                break;
                            case 1:
                                wcs[1] = v;
                                break;
                            default:
                                wcs[2] = v;
                                break;
                        }
                    }
                    // store into data using new idx mappings
                    ++nP;
                }
            } else {
                throw new IllegalStateException("expected no more lines in file");
            }
            line = in.readLine();
            ++nRead;
        }
        in.close();

        return data;
    }

}
