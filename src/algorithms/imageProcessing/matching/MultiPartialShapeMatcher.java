package algorithms.imageProcessing.matching;

/**
 * A class to Query a closed curve's shape against a database of closed curve shapes
 * and find the best <em>articulated</em> match for the Query
 * based upon algorithm in paper
 *  "Efficient Partial Shape Matching
 *     of Outer Contours: by Donoser
 *  implemented in PartialShapeMatcher.
 *
 * <pre>

 *  The databases construction have r.t.c. O(L*NS^3) where NS is the common number of points
 *  all curves are resampled to and L is the number of curves.
 *
 *  The Query construction has r.t.c. O(NS^2).
 *
 *  The Query search is roughly O(NS*topK*log(L)) where topK is the top results selected for each size of
 *  embedding.
 *
 *  Details of db creation:
 *    - the databases are instances of a nearest neighbor library indexer such as Google's ScANN or Meta's FAISS
 *      or java version of Voyager, etc.
 *    - shifted descriptor images (no subtraction or sum to integral image) are created.
 *       O(L*N_S^3) where L is the number of curves to store in databases.
 *    - stepping along the diagonal of a shifted descriptor image until reach minLength,
 *       -- each block starting at that diagonal index gets flattened into a vector and stored in a
 *          separate db for its length.
 *  Details of query:
 *    - descriptor image:
 *        O(NS^2)
 *    - search of Query against each db: < NS*log(L)
 *    - Salukwdze comparison: O(NS*topK*log(NS*topK))
 *
 * The above estimates do not include the computations for writing blocks to flattened vectors.
 * This algorithm could be implemented to use vectorization (intrinsics in SIMD or ISPC etc) and
 * vector shifts and total embeddings being composed of 8-wide or 16-wide vectors.
 *
 *</pre>
 */
public class MultiPartialShapeMatcher {

}
