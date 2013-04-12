two-point-correlation 
================================================================
finds clusters in data by looking for regions 
whose density is 3 times the background density (that is 3 sigma 
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

If you have data that has already been reduced so that there are few to 
no outliers between clusters, use:
    TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, n);
    clusterFinder.calculateBackground(true); <=====
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


-------------------
Performance Metrics
-------------------
Roughly determined by estimation and measured with a very small number of iterations
on computer with characteristics:
   64-bit Intel Core 2 Duo processor w/ speed 2 GHz w/ 4 MB L2 Cache
   3GB RAM and a bus speed of 800 MHz

JVM characteristics:
    J2SE 1.6, 64-bit HotSpot Server w/ initial heap size 512 MB and max heap size 1024 MB.

Measurements are from build from April 11, 2013

Note that the runtime (RT) dependency of voidFits is not well constrained due to the iterations needed
to construct a histogram and perform a fit.

    N   |  voidDens  big-oh        mem  | voidFits    M       mem   | Find Clusters  | Sys load |  Total RT
 points |  RT[sec]                 [kB] |  RT[sec]  density   [kB]  |    RT[sec]     | at start |    [sec]
        |                               |           points          |                |          |
--------------------------------------------------------------------------------------------------------------------------
     99 |      30    O(N^4)        30.4 |   30        539     38.9  |     0          |   0.27   |     60
        |                               |                           |                |          |

   1089 |       4  O(N*lg_2(N))   114.5 |    5       2199     39.5  |     0          |   0.56   |      9
        |                               |                           |                |          |
   1089 |       4  O(N^2 - N)     113.9 |    4       2159     39.5  |     0          |   1.45   |      8
        |                               |                           |                |          |

  14544 |    5502  O(N*lg_2(N))  1913.1 | 5507      36761    523.9  |    14          |   0.59   |  11023  = 184 min
        |                               |                           |                |          |
  13433 |    4513  O(N*lg_2(N))  1752.3 | 4515      33685    483.9  |    12          |   1.04   |   9040  = 151 min
        |                               |                           |                |          |
  13130 |    4163  O(N*lg_2(N))  1705.7 | 4165      32797    473.0  |    11          |   1.14   |   8339  = 139 min
        |                               |                           |                |          |
  13130 |    4208  O(N*lg_2(N))  1715.0 | 4210      32929    473.0  |    11          |   1.10   |   8429  = 140 min
        |                               |                           |                |          |
  12524 |    3739  O(N*lg_2(N))  1632.4 | 3741      31365    451.2  |    10          |   1.03   |   7490  = 125 min
        |                               |                           |                |          |
  12019 |    3271  O(N*lg_2(N))  1554.7 | 3273      29888    433.0  |    10          |   1.10   |   6554  = 109 min
        |                               |                           |                |          |
  10807 |    2457  O(N*lg_2(N))  1392.4 | 2458      26717    389.3  |     8          |   1.06   |   4923  = 82 min
        |                               |                           |                |          |

  Extrapolating to N=100,000 points, expect RT ~ 300 * 180 min using similar architecture and system load.
