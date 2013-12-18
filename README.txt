two-point-correlation 
================================================================

Find clusters in data in an unsupervised manner by calculating 
the two-point density function of rectangular voids, by fitting 
a GEV (generalized extreme value) curve to the density histogram 
to estimate the background, and then using that times a 
threshold factor (default is 2.5) to find groups of data points 
closer than the critical density.  The results are available as 
data and visualized through generated html plots.

More about the density distribution:

The location of the 'background' points in two dimensional 
space are likely Poisson, that is their locations in a fixed 
interval of space are independent of one another and occur 
randomly.
 
The 2-point void densities are formed into a histogram that 
is well fit by Generalized Extreme Value (GEV) curves. 
Extreme value distributions are used to describe the maximum 
or minimum of values drawn from a sample distribution that 
is essentially exponential.

There are 2 methods for determining clusters in this code:
  
  (1) For datasets in which there are background points:
  The peak of the GEV fit should represent the background density.  
  The clusters are then defined statistically as being 2 to 3 
  times 'above the background', that is having separations 2 to 
  3 times more dense than the background density. The code by 
  default uses a factor of 2.5, but methods are supplied to 
  allow the user to set the background to 2 or 3 instead, and 
  there's also a method to set the background manually.  The 
  later manual setting is useful for a case where perhaps one 
  determined the background density in one dataset and need to 
  apply that to a 2nd dataset which has the same background, 
  but is 'saturated' with foreground points.  
  
  (2) For datasets in which there are no background points:
  Datasets which are only points which should be in groups, and 
  essentially have no background points are referred to as 
  sparse background datasets.  For these datasets, the 
  background density is zero, so we define the level above the 
  background by the edges of the densities of the group.  This 
  edge density is already  much larger than the background so it 
  is the threshold density for membership already.  This 
  threshold density is the first x bin in a well formed 
  histogram of 2-point densities.

The code automatically determines which of method (1) and (2) to 
use.
  
  If the user has better knowledge of which should be applied, 
  can set that with:
     useFindMethodForDataWithoutBackgroundPoints() 
  or useFindMethodForDataWithBackgroundPoints()


The GEV curve contains 3 independent fitting parameters and the 
curve is an exponential combined with a polynomial, so it's 
resulting fitted parameters are not unique, but the curve is 
useful for characterizing the background point distribution by 
then integrating under the curve.
 
The points within a cluster are may have interesting 
distributions that can be better modeled after they've been 
found by these means.

Usage as an API:
{{{
To use the code with default settings:
  
       TwoPointCorrelation clusterFinder = new 
           TwoPointCorrelation(x, y, xErrors, yErrors, 
           totalNumberOfPoints);
  
       clusterFinder.calculateBackground();
       
       clusterFinder.findClusters();
  
  The results are available as group points or as convex hulls 
  surrounding the groups:
      int n = clusterFinder.getNumberOfGroups()
      
      int groupNumber = 0;

      To get the hull for groupId 0:
          ArrayPair hull0 = 
              clusterFinder.getGroupHull(groupNumber);

      To get the points in groupId 0:
          ArrayPair group0 = clusterFinder.getGroup(groupNumber)
      
      To plot the results:
          String plotFilePath = clusterFinder.plotClusters();

 If debugging is turned on, plots are generated and those file 
 paths are printed to standard out, and statements are printed 
 to standard out.
 
  To set the background density manually:
      TwoPointCorrelation clusterFinder = 
          new TwoPointCorrelation(x, y, xErrors, yErrors, 
          totalNumberOfPoints);
      clusterFinder.setBackground(0.03f, 0.003f);
      clusterFinder.findClusters();
      String plotFilePath = clusterFinder.plotClusters();
}}}


If the centers of the cluster hulls are needed for something 
else, seeds for a Voronoi diagram, for instance, one can use:
{{{
    ArrayPair seeds = clusterFinder.getHullCentroids();
}}}

The scatter plots and histograms below use d3.js (http://d3js.org d3)

Note that improvements in the histogram code is *in progress*.  
Currently datasets with a small number of points may have less 
than ideal solutions.

Note also that the code has the ability to refine a solution:  
that is to determine groups and then subtract them from the data 
and then re-determine the background density from the remaining 
points.  The ability is not enabled by default, but can be with 
the method setAllowRefinement().
More information is in docs/clustering_and_refinement.pdf

-------

The citation for use of this code in a publication is:
    `http://code.google.com/p/two-point-correlation/`, 
    Nichole King,  "Unsupervised Clustering Based Upon Voids in 
    Two-Point Correlation". March 15, 2013. <date accessed>


The license for using the source code is the MIT open source
license.  See LICENSE.txt


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

Measurements are from the April 27, 2013 code version tagged as v20131213.

    N    |  voidDens  RT complexity      mem  | Find Clusters   | Sys load |  Total RT
 points  |  RT[sec]                      [MB] |     RT[sec]     | at start |    [sec]
         |                                    |                 |          |
--------------------------------------------------------------------------------------------------------------------------
      99 |      0    O(N^4)              0.2  |     0           |   1.4    |      0

    1089 |      0  O((1.5*N/2)^1.8)        2  |     0           |   1.4    |      1

   13332 |      8  O(N^2.2)                1  |     1           |   1.4    |     12

  105100 |     10  O(N^2.2)              0.6  |    204          |   1.8    |    222
