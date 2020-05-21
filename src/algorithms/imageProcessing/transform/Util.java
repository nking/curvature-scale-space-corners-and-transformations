package algorithms.imageProcessing.transform;

import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class Util {
    
    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public static DenseMatrix rewriteInto3ColumnMatrix(PairIntArray xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.getN());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.getN(); i++) {
            xy.set(0, i, xyPairs.getX(i));
            xy.set(1, i, xyPairs.getY(i));
            xy.set(2, i, 1);
        }

        return xy;
    }
    
     /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public static DenseMatrix rewriteInto3ColumnMatrix(PairFloatArray xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.getN());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.getN(); i++) {
            xy.set(0, i, xyPairs.getX(i));
            xy.set(1, i, xyPairs.getY(i));
            xy.set(2, i, 1);
        }

        return xy;
    }

     /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public static DenseMatrix rewriteInto3ColumnMatrix(List<PairInt> xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.size());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.size(); i++) {
            PairInt p = xyPairs.get(i);
            xy.set(0, i, p.getX());
            xy.set(1, i, p.getY());
            xy.set(2, i, 1);
        }

        return xy;
    }

    /**
     * write a matrix of size mRows = 3, nCols = xyPairs.getN()
     * @param xyPairs
     * @return
     */
    public static DenseMatrix rewriteFirstItemInto3ColumnMatrix(List<List<PairInt>> xyPairs) {

        DenseMatrix xy = new DenseMatrix(3, xyPairs.size());

        // rewrite xyPairs into a matrix of size 3 X xy.getN();
        // row 0 is x
        // row 1 is y
        // row 2 is all 1's
        for (int i = 0; i < xyPairs.size(); i++) {
            List<PairInt> points = xyPairs.get(i);
            PairInt p = points.get(0);
            xy.set(0, i, p.getX());
            xy.set(1, i, p.getY());
            xy.set(2, i, 1);
        }

        return xy;
    }

    public static String toString(DenseMatrix m) {
        StringBuilder sb = new StringBuilder();
        for (int r = 0; r < m.numRows(); ++r) {
            StringBuilder line = new StringBuilder();
            for (int c = 0; c < m.numColumns(); ++c) {
                line.append(String.format("%.3e  ", m.get(r, c)));
            }
            sb.append(line).append("\n");
        }
        return sb.toString();
    }
}
