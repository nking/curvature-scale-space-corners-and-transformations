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
The project uses ant to build. You'll need to have installed on your 
computer java and ant. The other libraries are contained in the project.

To list the targets:
  ant

To build and package just the clustering package that uses a distance transform:
  ant package2

-------------------------
Performance Metrics
-------------------------

  The runtime complexity is ~ O(N_pixels) + ~O(N_points * lg2(N_points))

  rough space complexity estimate:
     creates arrays to hold data and the limiting value is roughly word size * width * height,
     that is platform word size * N_pixels.
     for a width of 5000 and height of 5000, code must be run with java arguments to
     increase the stack size or the data must be reduced in size... like knapsack,
     the code is using dynamic programmining using arrays that are as long as needed
     for capacity.

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
given to the code so will retain he specialization.
