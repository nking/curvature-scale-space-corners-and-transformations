package algorithms.imageProcessing.matching;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.signalProcessing.CurveResampler;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import com.spotify.voyager.jni.Index;

import java.util.*;

/**
 * A class to Query a closed curve's shape against a database of closed curve shapes
 * and find the best <em>articulated</em> match for the Query
 * based upon algorithm in paper
 *  "Efficient Partial Shape Matching
 *     of Outer Contours: by Donoser et al.
 *  implemented in PartialShapeMatcher.java
 *
 * <pre>

 *  The databases construction have r.t.c. O(L*(NS^2)) where NS is the common number of points
 *  all curves are resampled to and L is the number of curves.
 *
 *  The Query construction has r.t.c. O(NS^2) and the search of Indexed embeddings is
 *  O(topK*log(topK)*log(NS)*log(L*NS*(log(NS^2))) where topK is the top results.
 *
 *  Details of Index creation:
 *    - the indexes are instances of an approximate nearest neighbor library that indexes embeddings.
 *    The java enabled Spotify Voyager API is used for this purpose.
 *    - descriptor images for each closed curve are constructed following the algorithm by Donoser et al.
 *    in PartialShapeMatcher where the descriptor is part of the data structures built for
 *    for pairwise curve matching.  The descriptors are chord angles formed between relative orientations
 *    of a chord formed from the static reference point i, and a point j and a point j-1 as j marches
 *    around the curve.  N descriptors are formed with each point taking the role of the static reference
 *    point.  The final descriptor is then N X N where N is the number of points in the curve.
 *    - embeddings are made from each row of the descriptor, but offset to start at the reference point
 *    (which is the same as the row number).
 *    An embedding of length N is made for each row.
 *    And then log(N)-1 more embeddings are extracted for the curve starting at its reference point.
 *    e.g. if N=16 and row=2, we have the following embedding indexes (the embeddings would be vectors
 *    of the values for those descriptor chords at those indices):
 *    [2,3,4,5,6,7,8,9,10,11,12,13,14,15,0,1], [2,3,4,5,6,7,8,9], [2,3,4,5], [2,3]
 *    though the later would probably be dropped as being a smaller embedding than a minimum limit.
 *    And so the number of embeddings created is L*N*log(N).
 *    A Voyager Index is made for each unique length of the embeddings (i.e. log(N) or so Indexes are made)
 *    and the embeddings are stored in the Index for their length.
 *    => the runtime complexity for building the database is then O(L*N^2)
 *
 *  Details of query:
 *    - a descriptor image is made for the Query and it is chopped up into embeddings just as the Indexed
 *    target shapes were.
 *    The runtime complexity for that is O(N^2).
 *    - each of the N*log(N) query embeddings are searched against the Index of same length.
 *    Each ANN search for the topK is log(N_I) where N_I is the number of items stored in the Index.
 *    There are log(N) different Indexes so the runtime complexity is something like O(N * log(N) * log(N_I)).
 *    - The topK results are kept and ordered over all searches using the Salukwdze comparison from PartialShapeMatcher.
 *    The runtime complexity is O(N * log(N) * log(topK))
 *    ==> The overall runtime complexity of just the Query is then O(N^2).
 *
 * This algorithm could be implemented to use vectorization (intrinsics in SIMD or ISPC etc) and
 * vector shifts and total embeddings being composed of 8-wide or 16-wide vectors.
 *
 * The Voyager Index API uses vectorization.
 *
 *</pre>
 */
public class MultiPartialShapeMatcher {

    /**
     * the number of points that each closed curve will be scaled to.
     */
    protected final int curveDimension;

    /**
     * the smallest embedding length that will be used in creating internal indexes
     */
    protected final int minBlockSize;

    protected final Indexer indexer;

    protected final List<PairFloatArray> targetCurves;

    /**
     * class to build and hold Voyager indexes of target curves
     */
    protected static class Indexer {

        // see https://spotify.github.io/voyager/java/com/spotify/voyager/package-summary.html
        // NOTE that a SIGSEGV can occur for cases of adding to an Index in 1 batch a very large number
        // of embeddings of large length .  might be due to large use of the native stack memory
        // (near 1GB memory needed for example on OTC laptop where no
        // changes for java stack size have been made nor any for heap)

        /**
         * map of different size indexes to search.  key = embedding length, value = the Voyager indexer
         */
        public final Map<Integer, Index> indexMap;
        public final int nEmbeddings;
        public Indexer(int nEmbeddings) {
            this.nEmbeddings = nEmbeddings;
            this.indexMap = new HashMap<>(nEmbeddings);
        }

        public Index.QueryResults query(float[] embedding, int topK) {
            int len = embedding.length;
            Index index = indexMap.get(len);
            if (index == null) {
                throw new IllegalArgumentException("embedding length not found in Indexes: " + len);
            }
            if (index.getNumElements() < Integer.MAX_VALUE) {
                topK = Math.min(topK, (int) index.getNumElements());
            }
            return index.query(embedding, topK);
        }

        /**
         * search the index for topK approximate nearest neighbors for the given embeddings
         * @param embeddings an array of embeddings of same length
         * @param topK number of best results to return
         * @return topK nearest neighbors, approximately
         */
        public Index.QueryResults[] query(float[][] embeddings, int topK) {
            int len = embeddings[0].length;
            Index index = indexMap.get(len);
            if (index == null) {
                throw new IllegalArgumentException("embedding length not found in Indexes: " + len);
            }
            if (index.getNumElements() < Integer.MAX_VALUE) {
                topK = Math.min(topK, (int) index.getNumElements());
            }
            return index.query(embeddings, topK, -1);
        }

        public void addEmbeddingsIndex(float[][] embeddings, long[] ids) {
            int len = embeddings[0].length;
            if (!indexMap.containsKey(len)) {
                indexMap.put(len, new Index(Index.SpaceType.Euclidean, len));
            }

            // have to batch the embeddings to avoid a SIGSEGV that appears to be due to using too much
            //   memory in native methods
            //
            // default stack size can be obtained from command line:
            // java -XX:+PrintFlagsFinal -version | grep StackSize

            // DEBUG trying a subset.  succeeds with only 10 embeddings of length 30276.
            // succeeds w/ 100 embeddings of length 30276.
            // and also w/ 1000 embeddings of length 30276 but takes quite a long time.
            // fails for 30176 embeddings of length 30276.

            //TODO: batching could be improved w/ a look at memory and stack properties or allow user to
            // more configuration options

            Index index = indexMap.get(len);

            long nE = embeddings.length;
            if (nE * len < 10_000_000) {
                // -1 is numThreads.  If -1 (the default), the number of CPUs available on the current machine will be used.
                index.addItems(embeddings, ids, -1);
            } else {
                double nB = Math.ceil(nE * len / 1E7);
                int batchSize = (int)Math.ceil(nE / nB);
                int b = batchSize;
                for (int i = 0; i < nE; i+= batchSize) {
                    if (nE - i < batchSize) {
                        // batch is smaller than batchsize
                        b = (int)nE - i;
                    }
                    // copy from [i, i+b-1] into arrays to load
                    float[][] embeddings2 = copy(embeddings, i, i+b-1);
                    long[] ids2 = Arrays.copyOfRange(ids, i, i+b);
                    index.addItems(embeddings2, ids2, -1);
                }
            }
        }

        /**
         * copy rows i0 through i1, inclusive of array a
         * @param a
         * @param i0 first index to copy
         * @param i1 last index to copy, inclusive
         * @return copied rows i0 through i1, inclusive, of array a
         */
        private float[][] copy(float[][] a, int i0, int i1) {
            int len = a[i0].length;
            float[][] out = new float[i1 - i0 + 1][];
            for (int i = i0, j = 0; i <= i1; ++i, ++j) {
                out[j] = Arrays.copyOf(a[i0], len);
            }
            return out;
        }
    }

    /**
     * constructor.  r.t.c. is O(L*(NS^2)) where L is closedCurves.size() and NS = curveDimension.
     * @param curveDimension all curves are scaled to this length
     * @param minBlockSize the smallest vector length that will be created during this search.  e.g. should be at least 3
     * though practically you should consider minBlockSize >= 0.1*curveDimension.
     * @param closedCurves the shapes to scale to curveDimension, create embeddings for, and store in internal
     *                     Indexes for future queries.
     */
    public MultiPartialShapeMatcher(int curveDimension, int minBlockSize, List<PairFloatArray> closedCurves) {

        int L = closedCurves.size();

        List<PairFloatArray> targetCurves = new ArrayList<>();
        for (PairFloatArray curve : closedCurves) {
            targetCurves.add(curve.copy());
        }

        //[closedCurves.size][this.curveDimension-1][this.curveDimension-1]
        float[][][] descriptors = createDescriptors(closedCurves, curveDimension);
        assert(descriptors.length == L);
        //assert(descriptors[0].length == this.curveDimension);
        //assert(descriptors[0][0].length == this.curveDimension);

        Indexer indexer = new Indexer(L);
        addEmbeddingsToIndexes(descriptors, indexer, minBlockSize);
        this.indexer = indexer;
        this.targetCurves = targetCurves;
        this.curveDimension = curveDimension;
        this.minBlockSize = minBlockSize;

    }

    protected static void addEmbeddingsToIndexes(float[][][] descriptors, Indexer indexer, int minBlockSize) {
        int L = descriptors.length;

        int n = descriptors[0].length;

        int nr = (int)(Math.log(n)/Math.log(2));

        for (int r = 0; r < nr; ++r){
            int len = n/(1<<r);
            if (len < minBlockSize) {
                break;
            }
            // ids are a composite of iCurve and iDiag.
            // idx = iCurve * n + iDiag.
            // iCurve = idx / n;
            // iDiag = idx % n
            long[] ids = new long[L*n];
            float[][] embeddings = new float[L*n][len];
            int j = 0;
            for (int iCurve = 0; iCurve < L; ++iCurve) {
                for (int iDiag = 0; iDiag < n; ++iDiag) {
                    if ((iDiag + len) > n) {
                        int len1 = n - iDiag;
                        System.arraycopy(descriptors[iCurve][iDiag], iDiag, embeddings[j], 0, len1);
                        // wrap around
                        int len2 = iDiag + len - n;
                        System.arraycopy(descriptors[iCurve][iDiag],0, embeddings[j], len1, len2);
                    } else {
                        System.arraycopy(descriptors[iCurve][iDiag], iDiag, embeddings[j], 0, len);
                    }
                    ids[j] = iCurve * n + iDiag;
                    ++j;
                }
            }
            // store in indexer
            indexer.addEmbeddingsIndex(embeddings, ids);
        }
    }

    // return list of curve index, and offset to start match
    public static class Results {
        final List<PairFloatArray> dbCurves;
        final List<Integer> dbCurveIndexes;
        final List<Integer> offsetsQuery;
        final List<Integer> offsetsTargets;
        final List<Integer> matchingLengths;
        final List<Float> distances;
        public Results(int n) {
            dbCurves = new ArrayList<>();
            dbCurveIndexes = new ArrayList<>();
            offsetsTargets = new ArrayList<>();
            offsetsQuery = new ArrayList<>();
            matchingLengths = new ArrayList<>();
            distances = new ArrayList<>();
        }
        public void add(int dbCurveIndex, PairFloatArray dbCurve, int targetOffset, int queryOffset, int length, float distance) {
            dbCurves.add(dbCurve);
            dbCurveIndexes.add(dbCurveIndex);
            offsetsTargets.add(targetOffset);
            offsetsQuery.add(queryOffset);
            matchingLengths.add(length);
            distances.add(distance);
        }

        public List<PairFloatArray> getDBCurves() {
            return dbCurves;
        }

        public List<Integer> getDBCurveIndexes() {
            return dbCurveIndexes;
        }

        public List<Integer> getOffsetsQuery() {
            return offsetsQuery;
        }

        public List<Integer> getOffsetsTargets() {
            return offsetsTargets;
        }

        public List<Integer> getMatchingLengths() {
            return matchingLengths;
        }

        public List<Float> getDistances() {
            return distances;
        }
    }

    /**
     * search for topK nearest neighbors to queryCurve in the shapes indexes.
     * @param queryCurve the query closed curve
     * @param topK the number of best results to return
     * @return the nearest matched curves as the matching curves in their natural scale, the offset point
     * that the match begins, the length of the match as the number of matching points, and the distance
     * of the match.
     */
    public Results query(PairFloatArray queryCurve, final int topK) {

        PairFloatArray q1 = createScaledCurve(queryCurve, this.curveDimension);

        float[][] descriptor = PartialShapeMatcher.createDescriptorMatrix(q1);

        final float maxLength = descriptor.length;

        // build a TreeSet to sort topK results from queries:
        // calculate maxDist as euclidean distance of 2 vectors of chord differences, where a
        // single chord diff maximum possible value is up to 2*pi
        // so maxDiff = this.curveDimension*2*pi
        final float maxDiff = (float)(1.1 * maxLength * Math.PI * 2.);

        // storing search results as: distance, embedding length, id.
        // using a tree to keep to reduce the sort from O(n*log(n)) to O(n*(log(k)) by removing when tree size > k
        TreeSet<float[]> results = new TreeSet<>(new Comparator<float[]>() {
            @Override
            public int compare(float[] o1, float[] o2) {
                float diff1 = o1[0];
                float len1 = o1[1];
                float diff2 = o2[0];
                float len2 = o2[1];
                float d1 = calcSalukDist(diff1/len1, maxDiff, len1, maxLength);
                float d2 = calcSalukDist(diff2/len2, maxDiff, len2, maxLength);
                return Double.compare(d1, d2);
            }
        });

        /*
        create query embeddings:
        avoiding redundant search of whole curves by using row 0 of descriptor for length N embedding.
        The other embeddings:
        (1) row 0 of descriptor
        (2) the log(N)-1 embeddings of row 0-N-1
         */

        float[] embedding = Arrays.copyOf(descriptor[0], descriptor[0].length);
        Index.QueryResults res = indexer.query(embedding, topK);
        long[] labels = res.getLabels();
        float[] dists = res.getDistances();
        int queryOffset = 0;
        for (int i = 0; i < dists.length; ++i) {
            float[] result = new float[]{dists[i], embedding.length, labels[i], queryOffset};
            results.add(result);
            if (results.size() > topK) {
                results.removeLast();
            }
        }

        int L = this.targetCurves.size();

        int n = descriptor.length;

        int nr = (int)(Math.log(n)/Math.log(2));

        for (int r = 1; r < nr; ++r){
            int len = n/(1<<r);
            if (len < this.minBlockSize) {
                break;
            }
            // idsQ are iDiag, that is, the reference point, the offset from 0
            long[] idsQ = new long[n];
            float[][] embeddingsQ = new float[n][len];
            int j = 0;
            for (int iDiag = 0; iDiag < n; ++iDiag) {
                if ((iDiag + len) > n) {
                    int len1 = n - iDiag;
                    System.arraycopy(descriptor[iDiag], iDiag, embeddingsQ[j], 0, len1);
                    // wrap around
                    int len2 = iDiag + len - n;
                    System.arraycopy(descriptor[iDiag],0, embeddingsQ[j], len1, len2);
                } else {
                    System.arraycopy(descriptor[iDiag], iDiag, embeddingsQ[j], 0, len);
                }
                idsQ[j] = iDiag;
                ++j;
            }
            // there is a Index.QueryResults for every embeddingsQ, which has same index as idsQ
            Index.QueryResults[] indexResults = indexer.query(embeddingsQ, topK);

            for (int ii = 0; ii < indexResults.length; ++ii) {
                queryOffset = (int)idsQ[ii];
                Index.QueryResults indexRes = indexResults[ii];
                labels = indexRes.getLabels();
                dists = indexRes.getDistances();
                for (int i = 0; i < dists.length; ++i) {
                    // store in result, the dist, length, codedIndexLabel, query curve offset index
                    float[] result = new float[]{dists[i], len, labels[i], queryOffset};
                    results.add(result);
                    if (results.size() > topK) {
                        results.removeLast();
                    }
                }
            }
        }

        Results out = new Results(results.size());

        // resolve the curve index, shift index, and iDiag ref index to get the point offset index
        // and then use the scale to estimate the offset index within the original curve
        int i = 0;
        long id;
        for (float[] result : results) {
            //result: the dist, length, codedIndexLabel, query curve offset index
            id = (long)result[2];
            int iCurve = (int)(id / n);
            int iDiag = (int)(id % n);
            PairFloatArray curve = this.targetCurves.get(iCurve);
            // undo the scaling used in curve resampler:
            double factor = (this.curveDimension - 1.)/(curve.getN() - 1.);
            int targetOffset = (int) Math.round(iDiag * factor);
            int length = (int) Math.round(result[1] * factor);
            double qFactor = (this.curveDimension - 1.)/(queryCurve.getN() - 1.);
            queryOffset = (int) Math.round(result[3] * qFactor);
            /*{
                System.out.printf("iCurve=%d, off=%d, off_query=%d, len=%d, dist=%.3e, saluk=%.3e\n",
                        iCurve, targetOffset, queryOffset, length, result[0],
                        calcSalukDist(result[0]/result[1], maxDiff, length, maxLength));
            }*/
            out.add(iCurve, curve, targetOffset, queryOffset, length, result[0]);
            ++i;
        }

        return out;
    }

    protected static float calcSalukDist(float compChord, float maxChord,
                                          float length, float maxMatchable) {
        float d;
        if (maxChord == 0) {
            d = 0;
        } else {
            d = compChord / maxChord;
        }
        float f = 1.f - (length/maxMatchable);
        return f*f + d*d;
    }

    protected static PairFloatArray createScaledCurve(PairFloatArray p, int curveDimension) {
        //TODO: edit to store original number of points
        float[][] pxy = new float[2][p.getN()];
        pxy[0] = Arrays.copyOf(p.getX(), p.getN());
        pxy[1] = Arrays.copyOf(p.getY(), p.getN());
        float[][] xyOut = CurveResampler.resample(pxy, curveDimension);
        /*try {
            xyOut = CurveResampler.resample(pxy, this.curveDimension);
        } catch(Throwable t) {
            int t2 = 2;
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("float[][] pxy = new float[2][%d];\n", p.getN()));
            sb.append(String.format("pxy[0] = new float[]{"));
            for (int i = 0; i < pxy[0].length; ++i) {
                if (i > 0) {
                    sb.append(",");
                }
                sb.append(String.format("%.3f", pxy[0][i])).append("f");
            }
            sb.append("};\n");
            sb.append(String.format("pxy[1] = new float[]{"));
            for (int i = 0; i < pxy[1].length; ++i) {
                if (i > 0) {
                    sb.append(",");
                }
                sb.append(String.format("%.3f", pxy[1][i])).append("f");
            }
            sb.append("};\n");
            System.out.println(sb.toString());
        }*/
        PairFloatArray p2 = new PairFloatArray();
        for (int i = 0; i < xyOut[0].length; ++i) {
            p2.add(xyOut[0][i], xyOut[1][i]);
        }
        return p2;
    }

    public static List<PairFloatArray> convert(List<PairIntArray> p) {
        List<PairFloatArray> out = new ArrayList<>(p.size());
        for (PairIntArray q : p) {
            out.add(convert(q));
        }
        return out;
    }
    public static PairFloatArray convert(PairIntArray p) {
        PairFloatArray f = new PairFloatArray(p.getN());
        for (int i = 0; i < p.getN(); ++i) {
            f.add(p.getX(i), p.getY(i));
        }
        return f;
    }

    /**
     * create a descriptor image for each curve
     * @param closedCurves
     * @return descriptors array with dimensions: [closedCurves.size][this.curveDimension-1][this.curveDimension-1]
     */
    protected static float[][][] createDescriptors(List<PairFloatArray> closedCurves, int curveDimension) {
        float[][][] descriptors = new float[closedCurves.size()][][];
        // building the descriptors:
        // r.t.c. is O(n * curveDimension^2)
        for (int i = 0; i < closedCurves.size(); ++i) {
            PairFloatArray p = closedCurves.get(i);
            PairFloatArray p2 = createScaledCurve(p, curveDimension);
            /*if (true) {
                try {
                    plot(p2, System.currentTimeMillis());
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }*/
            float[][] a2 = PartialShapeMatcher.createDescriptorMatrix(p2);
            descriptors[i] = a2;
        }
        return descriptors;
    }

    private static String plot(PairFloatArray p, long fn) throws Exception {

        float[] x = Arrays.copyOf(p.getX(), p.getN());
        float[] y = Arrays.copyOf(p.getY(), p.getN());
        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
                x, y, x, y, "");

        return plot.writeFile(fn);
    }

}
