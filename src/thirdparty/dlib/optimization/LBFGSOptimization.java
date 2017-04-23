package thirdparty.dlib.optimization;

import algorithms.imageProcessing.util.MatrixUtil;
import java.util.Arrays;

/**
 * a port to java of dlib optimization method find_min.
 * The only search strategy ported is LBFGS, so the 
 * search strategy argument is specific, but could be 
 * generalized if more than one is implemented,
 * and then find_min abstracted...
 * 
 * the dlib project has copyright:
 * Copyright (C) 2008  Davis E. King (davis@dlib.net)
   License: Boost Software License   See LICENSE.txt for the full license.
 * @author nichole
 */
public class LBFGSOptimization {
  
    // f is a function, that when given the coefficients,
    //   returns the generated function
    
    // df is a function, that when given the coefficients and any
    //   other arguments, returns the coefficient gradient
    
    public static interface IFunction {
        
        // evaluates the objective at this set of coefficients.
        // e.g. could return chiSq
        double f (double[] coeffs);
        
        // estimates the coeficient gradient by comparing the
        //   data with the model
        double[] der(double[] coeffs);
    }
  
    public double findMin(LBFGSSearchStrategy searchStrategy,
        ObjectiveDeltaStopStrategy stopStrategy, 
        IFunction f,
        double[] x /* e.g. coeffs if a polynomial*/,
        double minF) {
      
        double[] s;
        
        double fValue = f.f(x);
        double[] g = f.der(x);

        if (!Double.isFinite(fValue)) {
            throw new IllegalStateException(
                "The objective function generated non-finite outputs");
        }
        for (double gv : g) {
            if (!Double.isFinite(gv)) {
                throw new IllegalStateException(
                    "The objective function generated non-finite outputs");
            }
        }
        
        long count = 0;
        
        while (stopStrategy.shouldContinueSearch(x, fValue, g) &&
            fValue > minF) {
        
            s = searchStrategy.get_next_direction(x, fValue, g);
            s = Arrays.copyOf(s, s.length);
            
            LineSearchFunction fls = new LineSearchFunction(
                f, x, s, fValue);
            
            LineSearchFunction flsDer = new LineSearchFunction(
                f, x, s, g);
            
            double alpha = line_search(
                fls, fValue, flsDer,
                MatrixUtil.multiplyByTranspose(g, s), // <== gradient * delX
                searchStrategy.get_wolfe_rho(), 
                searchStrategy.get_wolfe_sigma(), 
                minF, 100);

            
            // Take the search step indicated by the above line search
            //x += alpha*s;
            for (int i = 0; i < s.length; ++i) {
                x[i] += (alpha * s[i]);
            }
             
            //NLK: adding this for stop criteria
            fValue = f.f(x);
            g = f.der(x);
            
            count++;
            
            if (!Double.isFinite(fValue)) {
                throw new IllegalStateException(
                    "The objective function generated non-finite outputs");
            }
            for (double gv : g) {
                if (!Double.isFinite(gv)) {
                    throw new IllegalStateException(
                        "The objective function generated non-finite outputs");
                }
            }
        }

        return fValue;
    }
    
    //from optimization_line_search.h
    private double line_search(
        LineSearchFunction f,
        double f0, LineSearchFunction der,
        double d0, double rho, double sigma, 
        double minF, int maxIter) {
        
        assert(0 < rho && rho < sigma && sigma < 1 && maxIter > 0);
        /*
            "\tdouble line_search()" +
             "\n\tYou have given invalid arguments to this function"
            + "\n\t sigma:    " + sigma
            + "\n\t rho:      " + rho 
            + "\n\t max_iter: " + max_iter 
        */

        // The bracketing phase of this function is implemented according to block 2.6.2 from
        // the book Practical Methods of Optimization by R. Fletcher.   The sectioning 
        // phase is an implementation of 2.6.4 from the same book.

        // 1 <= tau1a < tau1b. Controls the alpha jump size during the 
        // bracketing phase of
        // the search.
        final double tau1a = 1.4;
        final double tau1b = 9;

        // it must be the case that 0 < tau2 < tau3 <= 1/2 for the algorithm 
        // to function correctly but the specific values of tau2 and tau3 
        // aren't super important.
        final double tau2 = 1.0/10.0;
        final double tau3 = 1.0/2.0;

        // Stop right away and return a step size of 0 if the gradient is 0 at the starting point
        if (Math.abs(d0) <= Math.abs(f0) * 2.22e-16)
            return 0;

        // Stop right away if the current value is good enough according to min_f
        if (f0 <= minF)
            return 0;

        // Figure out a reasonable upper bound on how large alpha can get.
        final double mu = (minF - f0)/(rho * d0);


        double alpha = 1;
        if (mu < 0)
            alpha = -alpha;
        alpha = putInRange(0, 0.65*mu, alpha);
        
        double last_alpha = 0;
        double last_val = f0;
        double last_val_der = d0;

        // The bracketing stage will find a range of points [a,b]
        // that contains a reasonable solution to the line search
        double a, b;

        // These variables will hold the values and derivatives of f(a) and f(b)
        double a_val, b_val, a_val_der, b_val_der;

        // This thresh value represents the Wolfe curvature condition
        final double thresh = Math.abs(sigma*d0);

        int itr = 0;
        // do the bracketing stage to find the bracket range [a,b]
        while (true) {
                        
            ++itr;
            final double val = f.operator(alpha, false);
            final double val_der = der.operator(alpha, true);
            
            // we are done with the line search since we found a value smaller
            // than the minimum f value
            if (val <= minF) {
                System.out.println("L1 alpha=" + alpha);
                return alpha;
            }

            if (val > f0 + rho*alpha*d0 || val >= last_val) {
            
                a_val = last_val;
                a_val_der = last_val_der;
                b_val = val;
                b_val_der = val_der;

                a = last_alpha;
                b = alpha;
                                
                break;
            }

            if (Math.abs(val_der) <= thresh) {
                return alpha;
            }

            // if we are stuck not making progress then quit with the current alpha
            if (last_alpha == alpha || itr >= maxIter) {
                return alpha;
            }

            if (val_der >= 0) {
                a_val = val;
                a_val_der = val_der;
                b_val = last_val;
                b_val_der = last_val_der;

                a = alpha;
                b = last_alpha;
                                
                break;
            }



            final double temp = alpha;
            // Pick a larger range [first, last].  We will pick the next alpha in that
            // range.
            double first, last;
            if (mu > 0) {
                first = Math.min(mu, alpha + tau1a*(alpha - last_alpha));
                last  = Math.min(mu, alpha + tau1b*(alpha - last_alpha));
            
            } else {
                
                first = Math.max(mu, alpha + tau1a*(alpha - last_alpha));
                last  = Math.max(mu, alpha + tau1b*(alpha - last_alpha));
            
            }
            

            // pick a point between first and last by doing some kind of interpolation
            if (last_alpha < alpha) {
                alpha = last_alpha + (alpha-last_alpha)
                    * poly_min_extrap(last_val, last_val_der, 
                    val, val_der, 1e10);
                            
            } else {
                alpha = alpha + (last_alpha-alpha)
                    *poly_min_extrap(val, val_der, 
                    last_val, last_val_der, 1e10);
                
            }
            
            alpha = putInRange(first, last, alpha);

            last_alpha = temp;

            last_val = val;
            last_val_der = val_der;
            
        }

        // Now do the sectioning phase from 2.6.4
        while (true) {
            
            ++itr;
            double first = a + tau2*(b-a);
            double last = b - tau3*(b-a);

            // use interpolation to pick alpha between first and last
            alpha = a + (b-a)
                *poly_min_extrap(a_val, a_val_der, b_val, b_val_der);
            alpha = putInRange(first,last,alpha);
            
            final double val = f.operator(alpha, false);
            final double val_der = der.operator(alpha, true);

            // we are done with the line search since we found a value smaller
            // than the minimum f value or we ran out of iterations.
            if (val <= minF || itr >= maxIter) {
                return alpha;
            }

            // stop if the interval gets so small that it isn't shrinking any more due to rounding error 
            if (a == first || b == last) {
                return b;
            }

            // If alpha has basically become zero then just stop.  Think of it like this,
            // if we take the largest possible alpha step will the objective function
            // change at all?  If not then there isn't any point looking for a better
            // alpha.
            final double max_possible_alpha = Math.max(Math.abs(a), Math.abs(b));
            if (Math.abs(max_possible_alpha*d0) <= Math.abs(f0) * 2.2e-16) {
                return alpha;
            }


            if (val > f0 + rho*alpha*d0 || val >= a_val) {
                b = alpha;
                b_val = val;
                b_val_der = val_der;
            } else {
                if (Math.abs(val_der) <= thresh) {
                    return alpha;
                }

                if ( (b-a)*val_der >= 0) {
                    
                    b = a;
                    b_val = a_val;
                    b_val_der = a_val_der;                    
                }

                a = alpha;
                a_val = val;
                a_val_der = val_der;                
            }
        }
    }
    
    public static class LineSearchFunction {
                
        private double scalarR = 0;
        private double[] start;
        private double[] direction;
        private IFunction funct;
        private double[] matrixR = null;
        
        public IFunction getFunction() {
            return funct;
        }
        
        public LineSearchFunction(
            IFunction f, double[] start_, double[] direction_) { 
            this.funct = f;
            this.start = start_;
            this.direction = direction_;
        }

        public LineSearchFunction(
            IFunction f, double[] start_, double[] direction_, double[] r) { 
            this.funct = f;
            this.start = start_;
            this.direction = direction_;
            this.matrixR = r;
        }
        
        //make_line_search_function(f, x, s, fValue),
        public LineSearchFunction(
            IFunction f, double[] start_, double[] direction_, double fValue) { 
            this.funct = f;
            this.start = start_;
            this.direction = direction_;
            this.scalarR = fValue;
        }
        
        public double operator(double x, boolean isGradient) {
            
            //return get_value(f(start + x*direction));
            
            double[] v0 = Arrays.copyOf(direction, direction.length);
            MatrixUtil.multiply(v0, x);
            
            for (int i = 0; i < start.length; ++i) {
                v0[i] += start[i];
            }
            
            
            if (isGradient) {
                double[] gValue = funct.der(v0);
                return get_value(gValue);
            } else {
                double fValue = funct.f(v0);
                return get_value(fValue);
            }            
        }
        
        public double get_value(double[] x) {
            
            if (matrixR != null) {
                //NOTE: not used
                matrixR = Arrays.copyOf(x, x.length);
            }
            
            double result = MatrixUtil.multiplyByTranspose(x, direction);
            
            return result;
        }
        
        private double get_value(double r) {
            // save a copy of this value for later
            if (scalarR > 0) {
                scalarR = r;
            }

            return r;
        }
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

    private double poly_min_extrap (
        double f0, double d0,
        double f1, double d1) {
        return poly_min_extrap(f0, d0, f1, d1, 1.);
    }
    
    private double poly_min_extrap (
        double f0, double d0,
        double f1, double d1, double limit) {
                
        final double n = 3*(f1 - f0) - 2*d0 - d1;
        final double e = d0 + d1 - 2*(f1 - f0);


        // find the minimum of the derivative of the polynomial

        double temp = Math.max(n*n - 3*e*d0,0.0);

        if (temp < 0)
            return 0.5;

        temp = Math.sqrt(temp);

        if (Math.abs(e) <= 2.2e-16)
            return 0.5;

        // figure out the two possible min values
        double x1 = (temp - n)/(3*e);
        double x2 = -(temp + n)/(3*e);

        // compute the value of the interpolating polynomial at these two points
        double y1 = f0 + d0*x1 + n*x1*x1 + e*x1*x1*x1;
        double y2 = f0 + d0*x2 + n*x2*x2 + e*x2*x2*x2;

        // pick the best point
        double x;
        if (y1 < y2)
            x = x1;
        else
            x = x2;

        // now make sure the minimum is within the allowed range of [0,limit] 
        return putInRange(0,limit,x);
    }

}
