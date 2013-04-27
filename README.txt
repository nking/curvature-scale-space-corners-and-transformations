two-point-correlation 
================================================================
finds clusters in data by looking for regions 
whose density is 2.5 times the background density (that is 2.5 sigma 
above 'noise').

Given x, y, xerror, and yerror data points, the code calculates 
the two-point density function of rectangular voids, fits that 
with a Generalized Extreme Value (GEV) distribution to determine 
the background density, and uses that as a threshold for finding 
clusters. The results are available as data and visualized through 
generated html plots.

The algorithm assumes that the data contains background points 
and clustered points whose 2 distributions are different from one 
another, but roughly homogeneous within their own. It assumes that 
the distributions do not need to be fit well for the rough range of 
the background density to be learned.

More specifically:

The location of the 'background' points in two dimensional space are 
likely Poisson, that is their locations in a fixed interval of space 
are independent of one another and occurred randomly. The areas between 
the smallest voids in such a distribution are well fit by 
Generalized Extreme Value (GEV) distributions. Extreme value 
distributions are used to describe the maximum or minimum of values 
drawn from a sample distribution that is essentially exponential.

The fits improve as N, the number of data points, increases.

The GEV curve contains 3 independent fitting parameters and the curve 
is an exponential combined with a polynomial, so it's resulting fitted 
parameters are not unique, but the curve is useful for characterizing 
the background point distribution and then integrating under the curve. 
The points within a cluster are distributed radially, and those are not 
modeled separately.

-----
Build
-----
The project uses ant to build. You'll need to have installed on your 
computer java and ant. The other libraries are contained in the project.

To list the targets:
  ant

To build:
  ant compile

To run the tests:
  ant runTests

-----------
Use the API
-----------
    TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, n);
    clusterFinder.calculateBackground();
    clusterFinder.findClusters();
    clusterFinder.calculateHullsOfClusters();
    String plotFilePath = clusterFinder.plotClusters();

This adjusts the the background calculation to be a smaller value to 
compensate for the sparsity of background points that were needed in 
the two point void calculation.

To set the background value yourself:
    TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, n);
    clusterFinder.setBackground(0.03f, 0.003f);
    clusterFinder.findClusters();
    clusterFinder.calculateHullsOfClusters();
    String plotFilePath = clusterFinder.plotClusters();
    

If the centers of the cluster hulls are needed for something else, 
seeds for a Voronoi diagram, for instance, one can use:
    float[] xSeeds = clusterFinder.getXHullCentroids();
    float[] ySeeds = clusterFinder.getYHullCentroids();

-------------------------
Use from the Command Line
-------------------------
Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.
Uses the default parameters at this time.
    java -jar bin/two-point-correlation.jar --file /path/to/file/fileName.txt

    optional flags:
       --twosigma
       --threesigma
       --background <value> (requires backgrounderror guesstimate at least)
       --backgrounderror <value>

-------------------
Performance Metrics  
-------------------

Roughly determined by estimation and measured with a very small number of iterations
on a computer with characteristics:
   64-bit Intel Core 2 Duo processor w/ speed 2 GHz w/ 4 MB L2 Cache
   3GB RAM and a bus speed of 800 MHz

JVM characteristics:
    J2SE 1.6, 64-bit HotSpot server w/ initial heap size 512 MB and max heap size 1024 MB.

Measurements are from the April 27, 2013 code version tagged as v20130427.

Note that the runtime (RT) dependency of voidFits on M is not well constrained due to the iterations needed
to construct a histogram and perform a fit.

    N   |  voidDens  big-oh              mem  | voidFits    M       mem   | Find Clusters  | Sys load |  Total RT
 points |  RT[sec]                       [kB] |  RT[sec]  density   [kB]  |    RT[sec]     | at start |    [sec]
        |                                     |           points          |                |          |
--------------------------------------------------------------------------------------------------------------------------
     99 |      22    O(N^4)              31.2 |   23        589     38.6  |     0          |   1.09   |     60
        |                                     |                           |                |          |

     99 |      25    O(N^4)              26.1 |   25        494     37.1  |     0          |   1.09   |     60
        |                                     |                           |                |          |

   1089 |       0  O((1.5*N/2)^1.8)      72.7 |    1       1380     87.5  |     0          |   1.32   |      1
        |                                     |                           |                |          |

   1089 |       0  O((1.5*N/2)^1.8)      72.2 |    1       1351     87.0  |     0          |   1.14   |      1
        |                                     |                           |                |          |

  12625 |      39  O(N^2.2)             909.3 |   40      17442    909.3  |     5          |   1.14   |     84
        |                                     |                           |                |          |

  13635 |      46  O(N^2.2)            1008.7 |   47      19317    826.1  |     6          |   1.32   |    103
        |                                     |                           |                |          |

  Extrapolating to N=100,000 points, expect roughly RT ~ 110 * 100 sec ~ 183 min ~ 3hr and will 
  need 10 MB of memory using similar architecture, system load, and jvm properties.
  (Note, the N relationships are *very* roughly approximated and will be improved at a later date.)
