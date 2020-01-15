package thirdparty.dlib.optimization;

import algorithms.imageProcessing.util.MatrixUtil;
import gnu.trove.list.array.TDoubleArrayList;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * adapted from dlib class 
 *   optimization/optimization_search_strategies.h

   Limited-memory BFGS is an optimization algorithm in the family of
   quasi-Newton methods (finds zeroes or local maxima and minima of functions)
   that approximates the
   Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm using a limited
   amount of computer memory. It is a popular algorithm for parameter
   estimation in machine learning
 * 
 * (add license here)
 * 
 */
public class LBFGSSearchStrategy {
    
    private int maxSize;
    private boolean beenUsed = false;
   
    //std::vector<matrix<double,0,1> >
    private double[] prev_derivative = null;
    private double[] prev_direction = null;
    private double[] prev_x = null;
    
    // alpha.length is same as data.size
    private double[] alpha = new double[0];
   
    private DataHelper dh_temp;
    
    //sequence<data_helper>::kernel_2a data;
    //  this could be replaced with equiv of dlib kernel_2a
    private final LinkedList<DataHelper> data;
    
    /*
    private sequence<data_helper>::kernel_2a data;
    
    typedef sequence_kernel_2<T,mem_manager>
            kernel_2a;
    
    #include "memory_manager_stateless/memory_manager_stateless_kernel_1.h" 
    
    for data, might need to include methods from algs.h and the memory manager
         header as needed.
    */
    
    public LBFGSSearchStrategy(int maxSize) {
        if (maxSize < 1) {
            throw new IllegalArgumentException("maxSize has to be > 0");
        }
        this.maxSize = maxSize;    
        
        //NOTE: if change to an extended LinkedHashSet, can set the capacity to maxSize
        data = new LinkedList<DataHelper>();
    }
    
    public double get_wolfe_rho() { return 0.01; }

    public double get_wolfe_sigma() { return 0.9; }

    public long get_max_line_search_iterations() { return 100; }

    double[] get_next_direction (
        double[] x, double fValue, 
        double[] funct_derivative) {
        
        prev_direction = Arrays.copyOf(funct_derivative, funct_derivative.length);
        MatrixUtil.multiply(prev_direction, -1.);

        if (!beenUsed) {
        
            beenUsed = true;
        
        } else {
        
            if (dh_temp == null) {
                dh_temp = new DataHelper();
                dh_temp.s = new double[x.length];
                Arrays.fill(dh_temp.s, Double.MAX_VALUE);
                dh_temp.y = new double[funct_derivative.length];
                Arrays.fill(dh_temp.y, Double.MAX_VALUE);
            }
        
            // add an element into the stored data sequence
            dh_temp.s = MatrixUtil.subtract(x, prev_x);
            dh_temp.y = MatrixUtil.subtract(funct_derivative, prev_derivative);

            double temp = MatrixUtil.multiplyByTranspose(dh_temp.s, dh_temp.y);
        
            // only accept this bit of data if temp isn't zero
            if (Math.abs(temp) > 1.e-7) {
        
                dh_temp.rho = 1./temp;
                
                dh_temp = dh_temp.copy();
                data.add(data.size(), dh_temp);                
            } else {
                    
                data.clear();                
            }

            if (data.size() > 0) {
                // This block of code is from algorithm 7.4 in the Nocedal book.
                            
                // makes total size(n) and erases all items after it
                alpha = resize(alpha, data.size());
               
                for (int i = data.size()-1; i > -1; --i) {    
                    
                    alpha[i] = 
                        data.get(i).rho * 
                        MatrixUtil.multiplyByTranspose(
                            data.get(i).s, prev_direction);
                    
                    //prev_direction -= alpha[i]*data[i].y;
                    double[] t = Arrays.copyOf(data.get(i).y, data.get(i).y.length);
                    MatrixUtil.multiply(t, alpha[i]);
                                        
                    for (int j = 0; j < prev_direction.length; ++j) {
                        prev_direction[j] -= t[j];
                    }
                }
                
                // Take a guess at what the first H matrix should be.  
                // This formula below is what is suggested
                // in the book Numerical Optimization by Nocedal and 
                // Wright in the chapter on Large Scale 
                // Unconstrained Optimization (in the L-BFGS section).
                double H_0 = 
                    1.0/data.get(data.size()-1).rho
                    / MatrixUtil.multiplyByTranspose(
                        data.get(data.size()-1).y, 
                        data.get(data.size()-1).y);

                H_0 = putInRange(0.001, 1000.0, H_0);

                MatrixUtil.multiply(prev_direction, H_0);
                    
                for (int i = 0; i < data.size(); ++i) {
                    
                    double beta = 
                        data.get(i).rho * 
                        MatrixUtil.multiplyByTranspose(
                        data.get(i).y, prev_direction);
                    
                    //prev_direction += data[i].s * (alpha[i] - beta);
                    
                    double[] t = Arrays.copyOf(data.get(i).s, data.get(i).s.length);
                    MatrixUtil.multiply(t, alpha[i] - beta);
                
                    for (int j = 0; j < prev_direction.length; ++j) {
                        prev_direction[j] += t[j];
                    }
                }                
            }
        }
        
        if (data.size() > maxSize) {
                    
            // remove the oldest element in the data sequence
            // defined in sequence/sequence_kernel_c.h
            remove(data, 0, dh_temp);
            
            //NOTE: remove is not invoked often so have decided to keep linkedlist.  
            // TODO: in future, extend LinkedHashSet and add an instance 
            //       variable in it to keep track of the last item
            //       in the list.  then change the data type of data to
            //       the extended LinkedHashSet.
            
        }

        prev_x = Arrays.copyOf(x, x.length);
        prev_derivative = Arrays.copyOf(funct_derivative, funct_derivative.length);
        
        if (prev_direction == null) {
            prev_direction = new double[x.length];
        }
                
        return prev_direction;
    }

    private void resize(TDoubleArrayList a, int size) {
        int n = a.size();
        if (n < size) {
            for (int i = n; i < size; ++i) {
                a.add(Double.NEGATIVE_INFINITY);
            }
        } else {
            for (int i = n; i < size; ++i) {
                a.set(i, Double.NEGATIVE_INFINITY);
            }
        }
    }
    
    private double[] resize(double[] a, int size) {
        int n = a.length;
        if (n < size) {
            double[] tmp = Arrays.copyOf(a, size);
            for (int i = n; i < size; ++i) {
                tmp[i] = Double.NEGATIVE_INFINITY;
            }
            return tmp;
        } else {
            for (int i = n; i < size; ++i) {
                a[i] = Double.NEGATIVE_INFINITY;
            }
        }
        return a;
    }
    
    private double putInRange(double a, double b, double val) {
        if (a < b) {
            if (val < a) {
                return a;
            } else if (val > b) {
                return b;
            }
        } else {
            if (val < b) {
                return b;
            } else if (val > a) {
                return a;
            }
        }
        return val;
    }

    static class DataHelper {
        double[] s = null;
        double[] y = null;
        double rho = Double.NEGATIVE_INFINITY;
        
        /*
        public DataHelper multiply(double f) {
            DataHelper tmp = copy();
            MatrixUtil.multiply(tmp.s, f);
            MatrixUtil.multiply(tmp.y, f);
            tmp.rho *= f;
            return tmp;
        }*/
        
        public DataHelper copy() {
            DataHelper tmp = new DataHelper();
            tmp.s = Arrays.copyOf(s, s.length);
            tmp.y = Arrays.copyOf(y, y.length);
            tmp.rho = rho;
            return tmp;
        }
        
        public void swap(DataHelper a, DataHelper b) {
            
            double[] tmp = Arrays.copyOf(a.s, a.s.length);
            a.s = b.s;
            b.s = tmp;
        
            tmp = Arrays.copyOf(a.y, a.y.length);
            a.y = b.y;
            b.y = tmp;
        
            double tmp2 = a.rho;
            a.rho = b.rho;
            b.rho = tmp2;
        }
    }

    //memory_manager_stateless_kernel_2<T,memory_manager<char>::kernel_2a>
    
    // from sequence_kernel_2.h
    private void remove (LinkedList<DataHelper> data,
        int pos, DataHelper item) {
        
        data.removeFirst();
        data.addFirst(item);
        data.removeLast();
        
        //NOTE, using svm requires additional logic here
    }
}
