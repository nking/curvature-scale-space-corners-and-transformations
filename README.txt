density based clustering code using a distance transform
================================================================

Find clusters in data in an unsupervised manner by calculating 
the distance transform to nearest points and the critical
two-point density from that.
The critical density combined with a factor above background
is used to define the critical separation of two points,
below which they are members of the same cluster.

Usage as an API:

    DTClusterFinder clusterFinder = new DTClusterFinder(points,
        imageWidth, imageHeight);
                
    clusterFinder.setToDebug();

    clusterFinder.calculateCriticalDensity();

    // alternatively, if the density is known, set instead of calculate:
    //clusterFinder.setCriticalDensity(dens);

    clusterFinder.findClusters();

    int nGroups = clusterFinder.getNumberOfClusters();

    List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
    for (int k = 0; k < nGroups; ++k) {
        Set<PairInt> set = clusterFinder.getCluster(k);
        groupList.add(set);
    }

-----
Build
-----
Requires java to run the code and ant to build it.
The other libraries are contained in the project.

To list the targets:
  ant

To build and package just the latest clustering package that uses a 
distance transform:
  ant package2

To compile the main source code:
  ant compile

To run all tests:
  ant runTests

To run a specific test:
  ant runTest -Dtest=package.TestName

-------------------------
Performance Metrics
-------------------------

The runtime complexity is roughly O(Npixels X log2(Npixels)) where 
Npixels is the width times height of the image. 
The space complexity is roughly platform word size X width X height, 
so for a width of 5000 and height of 5000, the code must be run with 
java arguments to increase the stack size or the data must be 
reduced in size... like knapsack, the code is using dynamic programmining 
using arrays that are as long as needed for capacity.

----------------------
Numerical Resolution
----------------------

floating point data can be converted to integers and scaled to keep the desired numerical
resolution before using this code.

For example, for data where one has an interest in exponential similarity functions,
one would want to apply the exponential operations to the data before use here and
scale the data for numerical resolution such that an integer holds the
significant difference between points
    exp(a - b) is exp(a)/exp(b) 

---------------------
Miscellaneous Advice
---------------------
To keep track of which point in Set<PairInt> points belongs to which pixel, one can extend
PairInt and add a field for the pixel index.  The final groups will have the same points
given to the code so will retain the specialization.
