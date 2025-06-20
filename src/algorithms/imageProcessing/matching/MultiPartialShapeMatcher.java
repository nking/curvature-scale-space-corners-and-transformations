package algorithms.imageProcessing.matching;

import algorithms.matrix.MatrixUtil;
import algorithms.signalProcessing.CurveResampler;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
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

}
