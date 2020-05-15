package algorithms.imageProcessing.util;

import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
    
/**
 *
 * @author nichole
 */
public class MatrixUtilTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MatrixUtilTest() {
    }
   
    public void testSparseMatrixMultiply() {
        
        /*
            m     *    n        result
        
         1  0  1    1  1  0    2  2  0
         0  0  0    1  1  0    0  0  0
         1  0  2    1  1  0    3  3  0
        */
        
        FlexCompRowMatrix m = new FlexCompRowMatrix(3, 3);
        m.set(0, 0, 1); m.set(0, 2, 1);
        m.set(2, 0, 1); m.set(2, 2, 2);
        
        FlexCompRowMatrix n = new FlexCompRowMatrix(3, 3);
        n.set(0, 0, 1); n.set(0, 1, 1);
        n.set(1, 0, 1); n.set(1, 1, 1);
        n.set(2, 0, 1); n.set(2, 1, 1);
        
        FlexCompRowMatrix expected = new FlexCompRowMatrix(3, 3);
        expected.set(0, 0, 2); expected.set(0, 1, 2);
        expected.set(2, 0, 3); expected.set(2, 1, 3);
        
        FlexCompRowMatrix result = MatrixUtil.sparseMatrixMultiply(m, n);
        
        Iterator<MatrixEntry> iter = result.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(expected.get(entry.row(), entry.column()), entry.get());
        }
        
        iter = expected.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(result.get(entry.row(), entry.column()), entry.get());
        }
    }
    
    public void testSparseMatrixMultiply2() {
        
        /*
            m     *    n        result
        
         1  0  1    1  1  0    2  2  0
         0  0  0    1  1  0    0  0  0
         1  0  2    1  1  0    3  3  0
        */
        
        FlexCompColMatrix m = new FlexCompColMatrix(3, 3);
        m.set(0, 0, 1); m.set(0, 2, 1);
        m.set(2, 0, 1); m.set(2, 2, 2);
        
        FlexCompColMatrix n = new FlexCompColMatrix(3, 3);
        n.set(0, 0, 1); n.set(0, 1, 1);
        n.set(1, 0, 1); n.set(1, 1, 1);
        n.set(2, 0, 1); n.set(2, 1, 1);
        
        FlexCompColMatrix expected = new FlexCompColMatrix(3, 3);
        expected.set(0, 0, 2); expected.set(0, 1, 2);
        expected.set(2, 0, 3); expected.set(2, 1, 3);
        
        FlexCompColMatrix result = MatrixUtil.sparseMatrixMultiply(m, n);
        
        Iterator<MatrixEntry> iter = result.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(expected.get(entry.row(), entry.column()), entry.get());
        }
        
        iter = expected.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(result.get(entry.row(), entry.column()), entry.get());
        }
    }
    
    public void testSparseMatrixSubtract() {
        
        /*
            m     *    n        result
        
         1  0  1    1  1  0    0  -1  1
         0  0  0    1  1  0    -1 -1  0
         1  0  2    1  1  0    0  -1  2
        */
        
        FlexCompRowMatrix m = new FlexCompRowMatrix(3, 3);
        m.set(0, 0, 1); m.set(0, 2, 1);
        m.set(2, 0, 1); m.set(2, 2, 2);
        
        FlexCompRowMatrix n = new FlexCompRowMatrix(3, 3);
        n.set(0, 0, 1); n.set(0, 1, 1);
        n.set(1, 0, 1); n.set(1, 1, 1);
        n.set(2, 0, 1); n.set(2, 1, 1);
        
        FlexCompRowMatrix expected = new FlexCompRowMatrix(3, 3);
        expected.set(0, 1, -1); expected.set(0, 2, 1);
        expected.set(1, 0, -1); expected.set(1, 1, -1);
        expected.set(2, 1, -1); expected.set(2, 2, 2);
        
        FlexCompRowMatrix result = MatrixUtil.sparseMatrixSubtract(m, n);
        
        Iterator<MatrixEntry> iter = result.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(expected.get(entry.row(), entry.column()), entry.get());
        }
        
        iter = expected.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(result.get(entry.row(), entry.column()), entry.get());
        }
    }
    
    public void testSparseMatrixSubtract2() {
        
        /*
            m     *    n        result
        
         1  0  1    1  1  0    0  -1  1
         0  0  0    1  1  0    -1 -1  0
         1  0  2    1  1  0    0  -1  2
        */
        
        FlexCompColMatrix m = new FlexCompColMatrix(3, 3);
        m.set(0, 0, 1); m.set(0, 2, 1);
        m.set(2, 0, 1); m.set(2, 2, 2);
        
        FlexCompColMatrix n = new FlexCompColMatrix(3, 3);
        n.set(0, 0, 1); n.set(0, 1, 1);
        n.set(1, 0, 1); n.set(1, 1, 1);
        n.set(2, 0, 1); n.set(2, 1, 1);
        
        FlexCompColMatrix expected = new FlexCompColMatrix(3, 3);
        expected.set(0, 1, -1); expected.set(0, 2, 1);
        expected.set(1, 0, -1); expected.set(1, 1, -1);
        expected.set(2, 1, -1); expected.set(2, 2, 2);
        
        FlexCompColMatrix result = MatrixUtil.sparseMatrixSubtract(m, n);
        
        Iterator<MatrixEntry> iter = result.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(expected.get(entry.row(), entry.column()), entry.get());
        }
        
        iter = expected.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            assertEquals(result.get(entry.row(), entry.column()), entry.get());
        }
    }
    
    private DenseMatrix[] readHalfMoonDataset() throws Exception {
        
        DenseMatrix data = new DenseMatrix(2, 100);
        DenseMatrix classes = new DenseMatrix(1, 100);
        
        data.set( 0 , 0 , 0.871318704123 );
        data.set( 1 , 0 , 0.490717552004 );
        data.set( 0 , 1 , 0.715472413369 );
        data.set( 1 , 1 , -0.458667853037 );
        data.set( 0 , 2 , 1.46253829024 );
        data.set( 1 , 2 , -0.386599306373 );
        data.set( 0 , 3 , -0.222520933956 );
        data.set( 1 , 3 , 0.974927912182 );
        data.set( 0 , 4 , 0.327699109739 );
        data.set( 1 , 4 , -0.240277997075 );
        data.set( 0 , 5 , 1.0 );
        data.set( 1 , 5 , 0.0 );
        data.set( 0 , 6 , 0.949055747011 );
        data.set( 1 , 6 , 0.315108218024 );
        data.set( 0 , 7 , 0.0 );
        data.set( 1 , 7 , 0.5 );
        data.set( 0 , 8 , 1.40478334312 );
        data.set( 1 , 8 , -0.414412623016 );
        data.set( 0 , 9 , 0.967294863039 );
        data.set( 1 , 9 , 0.25365458391 );
        data.set( 0 , 10 , 0.0960230259077 );
        data.set( 1 , 10 , 0.995379112949 );
        data.set( 0 , 11 , 0.427883339878 );
        data.set( 1 , 11 , -0.320172254597 );
        data.set( 0 , 12 , 1.09602302591 );
        data.set( 1 , 12 , -0.495379112949 );
        data.set( 0 , 13 , 0.198586378132 );
        data.set( 1 , 13 , -0.0981105304912 );
        data.set( 0 , 14 , 0.0320515775717 );
        data.set( 1 , 14 , 0.999486216201 );
        data.set( 0 , 15 , -0.900968867902 );
        data.set( 1 , 15 , 0.433883739118 );
        data.set( 0 , 16 , 1.15959989503 );
        data.set( 1 , 16 , -0.487181783414 );
        data.set( 0 , 17 , -0.761445958369 );
        data.set( 1 , 17 , 0.648228395308 );
        data.set( 0 , 18 , 0.073083242654 );
        data.set( 1 , 18 , 0.124732995121 );
        data.set( 0 , 19 , 1.03205157757 );
        data.set( 1 , 19 , -0.499486216201 );
        data.set( 0 , 20 , -0.623489801859 );
        data.set( 1 , 20 , 0.781831482468 );
        data.set( 0 , 21 , 1.76144595837 );
        data.set( 1 , 21 , -0.148228395308 );
        data.set( 0 , 22 , 0.345365054421 );
        data.set( 1 , 22 , 0.93846842205 );
        data.set( 0 , 23 , -0.284527586631 );
        data.set( 1 , 23 , 0.958667853037 );
        data.set( 0 , 24 , -0.404783343122 );
        data.set( 1 , 24 , 0.914412623016 );
        data.set( 0 , 25 , 1.87131870412 );
        data.set( 1 , 25 , 0.00928244799606 );
        data.set( 0 , 26 , 1.62348980186 );
        data.set( 1 , 26 , -0.281831482468 );
        data.set( 0 , 27 , 0.838088104892 );
        data.set( 1 , 27 , 0.545534901211 );
        data.set( 0 , 28 , 0.0184408430089 );
        data.set( 1 , 28 , 0.308841371299 );
        data.set( 0 , 29 , -0.871318704123 );
        data.set( 1 , 29 , 0.490717552004 );
        data.set( 0 , 30 , 0.222520933956 );
        data.set( 1 , 30 , 0.974927912182 );
        data.set( 0 , 31 , 1.83808810489 );
        data.set( 1 , 31 , -0.0455349012105 );
        data.set( 0 , 32 , -0.518392568311 );
        data.set( 1 , 32 , 0.855142763005 );
        data.set( 0 , 33 , 0.654634945579 );
        data.set( 1 , 33 , -0.43846842205 );
        data.set( 0 , 34 , 1.57211666012 );
        data.set( 1 , 34 , -0.320172254597 );
        data.set( 0 , 35 , 1.7183493501 );
        data.set( 1 , 35 , -0.195682550603 );
        data.set( 0 , 36 , 1.96729486304 );
        data.set( 1 , 36 , 0.24634541609 );
        data.set( 0 , 37 , 1.99179001382 );
        data.set( 1 , 37 , 0.372122838315 );
        data.set( 0 , 38 , 0.281650649902 );
        data.set( 1 , 38 , -0.195682550603 );
        data.set( 0 , 39 , 0.718349350098 );
        data.set( 1 , 39 , 0.695682550603 );
        data.set( 0 , 40 , 0.284527586631 );
        data.set( 1 , 40 , 0.958667853037 );
        data.set( 0 , 41 , 1.80141362187 );
        data.set( 1 , 41 , -0.0981105304912 );
        data.set( 0 , 42 , -0.718349350098 );
        data.set( 1 , 42 , 0.695682550603 );
        data.set( 0 , 43 , 0.161911895108 );
        data.set( 1 , 43 , -0.0455349012105 );
        data.set( 0 , 44 , 0.99794539275 );
        data.set( 1 , 44 , 0.0640702199807 );
        data.set( 0 , 45 , 0.967948422428 );
        data.set( 1 , 45 , -0.499486216201 );
        data.set( 0 , 46 , 0.761445958369 );
        data.set( 1 , 46 , 0.648228395308 );
        data.set( 0 , 47 , 1.28452758663 );
        data.set( 1 , 47 , -0.458667853037 );
        data.set( 0 , 48 , 0.623489801859 );
        data.set( 1 , 48 , 0.781831482468 );
        data.set( 0 , 49 , 0.032705136961 );
        data.set( 1 , 49 , 0.24634541609 );
        data.set( 0 , 50 , 0.518392568311 );
        data.set( 1 , 50 , 0.855142763005 );
        data.set( 0 , 51 , -0.0960230259077 );
        data.set( 1 , 51 , 0.995379112949 );
        data.set( 0 , 52 , 0.00205460724966 );
        data.set( 1 , 52 , 0.435929780019 );
        data.set( 0 , 53 , -0.967294863039 );
        data.set( 1 , 53 , 0.25365458391 );
        data.set( 0 , 54 , 0.926916757346 );
        data.set( 1 , 54 , 0.375267004879 );
        data.set( 0 , 55 , 1.99794539275 );
        data.set( 1 , 55 , 0.435929780019 );
        data.set( 0 , 56 , -0.345365054421 );
        data.set( 1 , 56 , 0.93846842205 );
        data.set( 0 , 57 , -0.949055747011 );
        data.set( 1 , 57 , 0.315108218024 );
        data.set( 0 , 58 , 0.840400104967 );
        data.set( 1 , 58 , -0.487181783414 );
        data.set( 0 , 59 , -0.926916757346 );
        data.set( 1 , 59 , 0.375267004879 );
        data.set( 0 , 60 , 0.572116660122 );
        data.set( 1 , 60 , 0.820172254597 );
        data.set( 0 , 61 , 1.94905574701 );
        data.set( 1 , 61 , 0.184891781976 );
        data.set( 0 , 62 , 0.404783343122 );
        data.set( 1 , 62 , 0.914412623016 );
        data.set( 0 , 63 , 0.672300890261 );
        data.set( 1 , 63 , 0.740277997075 );
        data.set( 0 , 64 , 0.159599895033 );
        data.set( 1 , 64 , 0.987181783414 );
        data.set( 0 , 65 , 0.801413621868 );
        data.set( 1 , 65 , 0.598110530491 );
        data.set( 0 , 66 , 0.128681295877 );
        data.set( 1 , 66 , 0.00928244799606 );
        data.set( 0 , 67 , 0.777479066044 );
        data.set( 1 , 67 , -0.474927912182 );
        data.set( 0 , 68 , 0.376510198141 );
        data.set( 1 , 68 , -0.281831482468 );
        data.set( 0 , 69 , 0.981559156991 );
        data.set( 1 , 69 , 0.191158628701 );
        data.set( 0 , 70 , -0.838088104892 );
        data.set( 1 , 70 , 0.545534901211 );
        data.set( 0 , 71 , -0.572116660122 );
        data.set( 1 , 71 , 0.820172254597 );
        data.set( 0 , 72 , -0.159599895033 );
        data.set( 1 , 72 , 0.987181783414 );
        data.set( 0 , 73 , 0.00820998617675 );
        data.set( 1 , 73 , 0.372122838315 );
        data.set( 0 , 74 , 0.900968867902 );
        data.set( 1 , 74 , 0.433883739118 );
        data.set( 0 , 75 , -0.99794539275 );
        data.set( 1 , 75 , 0.0640702199807 );
        data.set( 0 , 76 , 0.238554041631 );
        data.set( 1 , 76 , -0.148228395308 );
        data.set( 0 , 77 , 1.92691675735 );
        data.set( 1 , 77 , 0.124732995121 );
        data.set( 0 , 78 , 2.0 );
        data.set( 1 , 78 , 0.5 );
        data.set( 0 , 79 , -0.801413621868 );
        data.set( 1 , 79 , 0.598110530491 );
        data.set( 0 , 80 , 0.991790013823 );
        data.set( 1 , 80 , 0.127877161685 );
        data.set( 0 , 81 , 0.537461709759 );
        data.set( 1 , 81 , -0.386599306373 );
        data.set( 0 , 82 , 0.0509442529893 );
        data.set( 1 , 82 , 0.184891781976 );
        data.set( 0 , 83 , -1.0 );
        data.set( 1 , 83 , 1.22464679915e-16 );
        data.set( 0 , 84 , 0.595216656878 );
        data.set( 1 , 84 , -0.414412623016 );
        data.set( 0 , 85 , 1.34536505442 );
        data.set( 1 , 85 , -0.43846842205 );
        data.set( 0 , 86 , -0.672300890261 );
        data.set( 1 , 86 , 0.740277997075 );
        data.set( 0 , 87 , 1.22252093396 );
        data.set( 1 , 87 , -0.474927912182 );
        data.set( 0 , 88 , 1.98155915699 );
        data.set( 1 , 88 , 0.308841371299 );
        data.set( 0 , 89 , -0.0320515775717 );
        data.set( 1 , 89 , 0.999486216201 );
        data.set( 0 , 90 , -0.981559156991 );
        data.set( 1 , 90 , 0.191158628701 );
        data.set( 0 , 91 , -0.462538290241 );
        data.set( 1 , 91 , 0.886599306373 );
        data.set( 0 , 92 , 0.903976974092 );
        data.set( 1 , 92 , -0.495379112949 );
        data.set( 0 , 93 , -0.991790013823 );
        data.set( 1 , 93 , 0.127877161685 );
        data.set( 0 , 94 , 1.67230089026 );
        data.set( 1 , 94 , -0.240277997075 );
        data.set( 0 , 95 , 0.0990311320976 );
        data.set( 1 , 95 , 0.0661162608824 );
        data.set( 0 , 96 , 1.51839256831 );
        data.set( 1 , 96 , -0.355142763005 );
        data.set( 0 , 97 , 0.462538290241 );
        data.set( 1 , 97 , 0.886599306373 );
        data.set( 0 , 98 , 1.9009688679 );
        data.set( 1 , 98 , 0.0661162608824 );
        data.set( 0 , 99 , 0.481607431689 );
        data.set( 1 , 99 , -0.355142763005 );
            
        classes.set( 0 , 0 , 0 );
        classes.set( 0 , 1 , 1 );
        classes.set( 0 , 2 , 1 );
        classes.set( 0 , 3 , 0 );
        classes.set( 0 , 4 , 1 );
        classes.set( 0 , 5 , 0 );
        classes.set( 0 , 6 , 0 );
        classes.set( 0 , 7 , 1 );
        classes.set( 0 , 8 , 1 );
        classes.set( 0 , 9 , 0 );
        classes.set( 0 , 10 , 0 );
        classes.set( 0 , 11 , 1 );
        classes.set( 0 , 12 , 1 );
        classes.set( 0 , 13 , 1 );
        classes.set( 0 , 14 , 0 );
        classes.set( 0 , 15 , 0 );
        classes.set( 0 , 16 , 1 );
        classes.set( 0 , 17 , 0 );
        classes.set( 0 , 18 , 1 );
        classes.set( 0 , 19 , 1 );
        classes.set( 0 , 20 , 0 );
        classes.set( 0 , 21 , 1 );
        classes.set( 0 , 22 , 0 );
        classes.set( 0 , 23 , 0 );
        classes.set( 0 , 24 , 0 );
        classes.set( 0 , 25 , 1 );
        classes.set( 0 , 26 , 1 );
        classes.set( 0 , 27 , 0 );
        classes.set( 0 , 28 , 1 );
        classes.set( 0 , 29 , 0 );
        classes.set( 0 , 30 , 0 );
        classes.set( 0 , 31 , 1 );
        classes.set( 0 , 32 , 0 );
        classes.set( 0 , 33 , 1 );
        classes.set( 0 , 34 , 1 );
        classes.set( 0 , 35 , 1 );
        classes.set( 0 , 36 , 1 );
        classes.set( 0 , 37 , 1 );
        classes.set( 0 , 38 , 1 );
        classes.set( 0 , 39 , 0 );
        classes.set( 0 , 40 , 0 );
        classes.set( 0 , 41 , 1 );
        classes.set( 0 , 42 , 0 );
        classes.set( 0 , 43 , 1 );
        classes.set( 0 , 44 , 0 );
        classes.set( 0 , 45 , 1 );
        classes.set( 0 , 46 , 0 );
        classes.set( 0 , 47 , 1 );
        classes.set( 0 , 48 , 0 );
        classes.set( 0 , 49 , 1 );
        classes.set( 0 , 50 , 0 );
        classes.set( 0 , 51 , 0 );
        classes.set( 0 , 52 , 1 );
        classes.set( 0 , 53 , 0 );
        classes.set( 0 , 54 , 0 );
        classes.set( 0 , 55 , 1 );
        classes.set( 0 , 56 , 0 );
        classes.set( 0 , 57 , 0 );
        classes.set( 0 , 58 , 1 );
        classes.set( 0 , 59 , 0 );
        classes.set( 0 , 60 , 0 );
        classes.set( 0 , 61 , 1 );
        classes.set( 0 , 62 , 0 );
        classes.set( 0 , 63 , 0 );
        classes.set( 0 , 64 , 0 );
        classes.set( 0 , 65 , 0 );
        classes.set( 0 , 66 , 1 );
        classes.set( 0 , 67 , 1 );
        classes.set( 0 , 68 , 1 );
        classes.set( 0 , 69 , 0 );
        classes.set( 0 , 70 , 0 );
        classes.set( 0 , 71 , 0 );
        classes.set( 0 , 72 , 0 );
        classes.set( 0 , 73 , 1 );
        classes.set( 0 , 74 , 0 );
        classes.set( 0 , 75 , 0 );
        classes.set( 0 , 76 , 1 );
        classes.set( 0 , 77 , 1 );
        classes.set( 0 , 78 , 1 );
        classes.set( 0 , 79 , 0 );
        classes.set( 0 , 80 , 0 );
        classes.set( 0 , 81 , 1 );
        classes.set( 0 , 82 , 1 );
        classes.set( 0 , 83 , 0 );
        classes.set( 0 , 84 , 1 );
        classes.set( 0 , 85 , 1 );
        classes.set( 0 , 86 , 0 );
        classes.set( 0 , 87 , 1 );
        classes.set( 0 , 88 , 1 );
        classes.set( 0 , 89 , 0 );
        classes.set( 0 , 90 , 0 );
        classes.set( 0 , 91 , 0 );
        classes.set( 0 , 92 , 1 );
        classes.set( 0 , 93 , 0 );
        classes.set( 0 , 94 , 1 );
        classes.set( 0 , 95 , 1 );
        classes.set( 0 , 96 , 1 );
        classes.set( 0 , 97 , 0 );
        classes.set( 0 , 98 , 1 );
        classes.set( 0 , 99 , 1 );
        return new DenseMatrix[]{data, classes};
    }
    
    public void testDeterminant() {
        
        /**
         | 1  -5  2 | 
         | 7   3  4 |  
         | 2   1  5 |
         
                  | 3 4 |         | 7 4 |         | 7 3 |
           =  1 * | 1 5 |  +  5 * | 2 5 |  +  2 * | 2 1 |  = 11 + 135 + 2 = 148
         
         */
         
        double[][] m = new double[3][3];
        m[0] = new double[]{1, -5, 2};
        m[1] = new double[]{7, 3, 4};
        m[2] = new double[]{2, 1, 5};
        
        double det = MatrixUtil.determinant(m);
        
        assertEquals(148., det);
        
    }
    
     public void testInverse() throws Exception {
        /**
         * e.g.    | 1  0  -1 |
         *         |-2  1   0 |
         *         | 1 -1   2 |
         * 
         * det = 1
         * 
         * cofactor = 
         *            | (2)   -(-4)   (2-1) |
         *            | 1(-1)  (2+1) -(-1)  |
         *            |  1    -(0-2)    1   |
         * 
         *  cofactor^T = | 2  1  1 |
         *               | 4  3  2 |
         *               | 1  1  1 |
         * 
         *   inv = (1/det) * cofactor^T
         */
        System.out.println("testInverse");

        double[][] m = new double[3][3];
        m[0] = new double[]{1, -2, 1};
        m[1] = new double[]{0, 1, -1};
        m[2] = new double[]{-1, 0, 2};
        
        double[][] result = MatrixUtil.inverse(m);

        double[][] e = new double[3][3];
        e[0] = new double[]{2, 4, 1};
        e[1] = new double[]{1, 3, 1};
        e[2] = new double[]{1, 2, 1};

        assertTrue(result.length == e.length);
        assertTrue(result[0].length == e[0].length);
        
        for (int i = 0; i < e.length; i++) {
            double[] a = result[i];
            double[] b = e[i];
            assertTrue(Arrays.equals(a, b));
        }
    }
}
