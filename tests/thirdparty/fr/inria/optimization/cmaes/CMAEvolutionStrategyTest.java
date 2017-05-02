package thirdparty.fr.inria.optimization.cmaes;

import java.util.Arrays;
import thirdparty.fr.inria.optimization.cmaes.fitness.IObjectiveFunction;
import junit.framework.TestCase;

/**
 *
 */
public class CMAEvolutionStrategyTest extends TestCase {
    
    public CMAEvolutionStrategyTest() {
    }
   
    /*
    some of the unit tests are ported from gsl project file
     * multimin/test.c
     * downloaded from 
     * http://mirror.team-cymru.org/gnu/gsl/
     * the code has license 
     * GNU General Public License.
     * https://www.gnu.org/software/gsl/
     * http://www.gnu.org/copyleft/gpl.html
    */
    
    public void test2() {
        
        rosenbrock(new double[]{-1.2, 1.0});
        rosenbrock(new double[] {2., 2.0});
        
        roth();
        wood();
        simpleAbs();
        parabola();
    }
   
    public void rosenbrock(double[] init) {
        
        IObjectiveFunction fitfun = new Rosenbrock();

		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
        cma.setDimension(2); // overwrite some loaded properties
		cma.setInitialX(init);
        cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		cma.options.stopFitness = 1e-9;       // optional setting

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// core iteration step 
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				while (!fitfun.isFeasible(pop[i]))   //    test whether solution is feasible,  
					pop[i] = cma.resampleSingle(i);  //       re-sample solution until it is feasible  
				fitness[i] = fitfun.valueOf(pop[i]); //    compute fitness value, where fitfun
			}	                                     //    is the function to be minimized
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution

			// output to console and files
			cma.writeToDefaultFiles();
			int outmod = 150;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1)
				cma.println(); 
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 
        
		// final output
		cma.writeToDefaultFiles(1);
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
				+ " at evaluation " + cma.getBestEvaluationNumber());
			
		// we might return cma.getBestSolution() or cma.getBestX()

        
        CMASolution soln = cma.getBestSolution();
        
        System.out.println("final=" + Arrays.toString(soln.getX()));
        
        assertTrue(Math.abs(soln.getX()[0] - 1.) < 0.01);
        assertTrue(Math.abs(soln.getX()[1] - 1.) < 0.01);

    }
    
    public void roth() {
        
        double[] init = new double[] {4.5, 3.5};
        
        IObjectiveFunction fitfun = new Roth();

		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
        cma.setDimension(2); // overwrite some loaded properties
		cma.setInitialX(init);
        cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		cma.options.stopFitness = 1e-9;       // optional setting

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// core iteration step 
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				while (!fitfun.isFeasible(pop[i]))   //    test whether solution is feasible,  
					pop[i] = cma.resampleSingle(i);  //       re-sample solution until it is feasible  
				fitness[i] = fitfun.valueOf(pop[i]); //    compute fitness value, where fitfun
			}	                                     //    is the function to be minimized
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution

			// output to console and files
			cma.writeToDefaultFiles();
			int outmod = 150;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1)
				cma.println(); 
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 
        
		// final output
		cma.writeToDefaultFiles(1);
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
				+ " at evaluation " + cma.getBestEvaluationNumber());
			
		// we might return cma.getBestSolution() or cma.getBestX()

        
        CMASolution soln = cma.getBestSolution();
        
        System.out.println("final=" + Arrays.toString(soln.getX()));
        
        assertTrue(Math.abs(soln.getX()[0] - 5.) < 0.01);
        assertTrue(Math.abs(soln.getX()[1] - 4.) < 0.01);

    }
    
    public void wood() {
        
        double[] init = new double[] {-3.0, -1.0, -3.0, -1.0};
        
        IObjectiveFunction fitfun = new Wood();

		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
        cma.setDimension(4); // overwrite some loaded properties
		cma.setInitialX(init);
        cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		cma.options.stopFitness = 1e-9;       // optional setting

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// core iteration step 
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				while (!fitfun.isFeasible(pop[i]))   //    test whether solution is feasible,  
					pop[i] = cma.resampleSingle(i);  //       re-sample solution until it is feasible  
				fitness[i] = fitfun.valueOf(pop[i]); //    compute fitness value, where fitfun
			}	                                     //    is the function to be minimized
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution

			// output to console and files
			cma.writeToDefaultFiles();
			int outmod = 150;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1)
				cma.println(); 
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 
        
		// final output
		cma.writeToDefaultFiles(1);
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
				+ " at evaluation " + cma.getBestEvaluationNumber());
			
		// we might return cma.getBestSolution() or cma.getBestX()

        
        CMASolution soln = cma.getBestSolution();
        
        System.out.println("final=" + Arrays.toString(soln.getX()));
        
        assertTrue(Math.abs(soln.getX()[0] - 1.) < 0.01);
        assertTrue(Math.abs(soln.getX()[1] - 1.) < 0.01);
        assertTrue(Math.abs(soln.getX()[2] - 1.) < 0.01);
        assertTrue(Math.abs(soln.getX()[3] - 1.) < 0.01);

    }
    
    public void simpleAbs() {
        
        double[] init = new double[] {1., 2.0};
        
        IObjectiveFunction fitfun = new SimpleAbs();

		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
        cma.setDimension(2); // overwrite some loaded properties
		cma.setInitialX(init);
        cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		cma.options.stopFitness = 1e-9;       // optional setting

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// core iteration step 
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				while (!fitfun.isFeasible(pop[i]))   //    test whether solution is feasible,  
					pop[i] = cma.resampleSingle(i);  //       re-sample solution until it is feasible  
				fitness[i] = fitfun.valueOf(pop[i]); //    compute fitness value, where fitfun
			}	                                     //    is the function to be minimized
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution

			// output to console and files
			cma.writeToDefaultFiles();
			int outmod = 150;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1)
				cma.println(); 
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 
        
		// final output
		cma.writeToDefaultFiles(1);
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
				+ " at evaluation " + cma.getBestEvaluationNumber());
			
		// we might return cma.getBestSolution() or cma.getBestX()

        
        CMASolution soln = cma.getBestSolution();
        
        System.out.println("final=" + Arrays.toString(soln.getX()));
        
        assertTrue(Math.abs(soln.getX()[0] - 1.) < 0.01);
        assertTrue(Math.abs(soln.getX()[1] - 2.) < 0.01);

    }
    
    public void parabola() {
        
        double[] init = new double[] {1., 1., 1.};
        
        IObjectiveFunction fitfun = new Parabola();

		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
        cma.setDimension(3); // overwrite some loaded properties
		cma.setInitialX(init);
        cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		cma.options.stopFitness = 1e-9;       // optional setting

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// core iteration step 
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				while (!fitfun.isFeasible(pop[i]))   //    test whether solution is feasible,  
					pop[i] = cma.resampleSingle(i);  //       re-sample solution until it is feasible  
				fitness[i] = fitfun.valueOf(pop[i]); //    compute fitness value, where fitfun
			}	                                     //    is the function to be minimized
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution

			// output to console and files
			cma.writeToDefaultFiles();
			int outmod = 150;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1)
				cma.println(); 
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 
        
		// final output
		cma.writeToDefaultFiles(1);
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
				+ " at evaluation " + cma.getBestEvaluationNumber());
			
		// we might return cma.getBestSolution() or cma.getBestX()

        
        CMASolution soln = cma.getBestSolution();
        
        System.out.println("final=" + Arrays.toString(soln.getX()));
        
        assertTrue(Math.abs(soln.getX()[0] - 3.) < 0.01);
        assertTrue(Math.abs(soln.getX()[1] - 2.) < 0.01);
        assertTrue(Math.abs(soln.getX()[2] - 1.) < 0.01);

    }
    
    class Rosenbrock implements IObjectiveFunction { 

        public double valueOf(double[] x) {
            double res = 0;
            for (int i = 0; i < x.length - 1; ++i) {
                res += 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1])
                    + (x[i] - 1.) * (x[i] - 1.);
            }
            return res;
        }

        public boolean isFeasible(double[] x) {
            return true;
        } // entire R^n is feasible
    }
    
    private static class Roth implements IObjectiveFunction {
        
        public Roth() {}

        public double valueOf(double[] coeffs) {
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
            double b = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
            
            double c = -2 + v * (10 - 3 * v);
            double d = -14 + v * (2 + 3 * v);
            
            double m = a * a + b * b;
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }
        
        public boolean isFeasible(double[] x) {
            return true;
        }
    }

    private static class Wood implements IObjectiveFunction {
        
        public Wood() {}

        public double valueOf(double[] coeffs) {
            
            double u1 = coeffs[0];
            double u2 = coeffs[1];
            double u3 = coeffs[2];
            double u4 = coeffs[3];

            double t1 = u1 * u1 - u2;
            double t2 = u3 * u3 - u4;
        
            double m = 100 * t1 * t1 + (1 - u1) * (1 - u1)
                + 90 * t2 * t2 + (1 - u3) * (1 - u3)
                + 10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4))
                + 19.8 * (1 - u2) * (1 - u4);
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }
        
        public boolean isFeasible(double[] x) {
            return true;
        }
    }

    private static class SimpleAbs implements IObjectiveFunction {
        
        public SimpleAbs() {}

        public double valueOf(double[] coeffs) {
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = u - 1;
            double b = v - 2;
  
            double sign0 = u - 1;
            if (sign0 >= 0) {
                sign0 = 1;
            } else {
                sign0 = -1;
            }
            
            double sign1 = v - 2;
            if (sign1 >= 0) {
                sign1 = 1;
            } else {
                sign1 = -1;
            }
            
            double m = Math.abs(a) + Math.abs(b);
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }
        
        public boolean isFeasible(double[] x) {
            return true;
        }
    }
    
    private static class Parabola implements IObjectiveFunction {
        
        //test data from: http://rosettacode.org/wiki/Polynomial_Fitting
        double[] xp = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] yp = new double[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
    
        public Parabola() {}

        public double valueOf(double[] coeffs) {
            
            double[] gen = new double[xp.length];
            double[] gradient = new double[coeffs.length];
            double[] diffY = new double[xp.length];
       
            double sumDiff = calcGradient(coeffs, gen, gradient, diffY);
            
            return sumDiff;
        }
        
        public boolean isFeasible(double[] x) {
            return true;
        }
        
        void generatePolynomial(double[] coeffs, 
            double[] xPoly, double[] gen) {

            for (int i = 0; i < 11; ++i) {
                gen[i] = 0.;
            }        

            for (int i = 0; i < 11; ++i) {
                double x2 = 1.;
                for (int j = 2; j > -1; j--) {
                    double c = coeffs[j];
                    gen[i] += (c * x2);
                    x2 *= xPoly[i];
                }
            }
        }
        
        double calcGradient(double[] coeffs, 
            double[] gen, double[] outputCoeffGrad,
            double[] outputDiffY) {

            generatePolynomial(coeffs, xp, gen);

            double sumDiff = 0;
            for (int i = 0; i < 11; ++i) {
                outputDiffY[i] = (gen[i] - yp[i]);
                sumDiff += (outputDiffY[i] * outputDiffY[i]);
            }
            sumDiff = Math.sqrt(sumDiff);

            for (int j = 0; j < 3; ++j) {
                outputCoeffGrad[j] = 0.;
            }
            
            for (int i = 0; i < 11; ++i) {
                double x2 = 1;
                double dyAtX = outputDiffY[i];
                for (int j = 2; j > -1; j--) {
                    int varIdx = 3 - j - 1;

                    //dy * (dc/dy)
                    outputCoeffGrad[varIdx] += (dyAtX/x2);
                                        
                    x2 *= xp[i];

                    if (x2 == 0.0) {
                        break;
                    }
                }
            }

            for (int j = 0; j < 3; ++j) {
                outputCoeffGrad[j] /= xp.length;
            }
            
            return sumDiff;
        }
    }
}
