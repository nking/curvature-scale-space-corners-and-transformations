# Requirements: #
  * Running the code requires **[java](https://java.com/en/download/index.jsp)**
  * Building the code additionally requires **[ant](http://ant.apache.org/bindownload.cgi)**

# Running existing build: #
There is a built, packaged form of the current project in the dist directory in jar archive file format.   That jar file has a few command line runners in it to run some of the code's features through the command line.
To see the top level usage that presents more detailed options:
> ` java -jar dist/scalespace.jar`

# Building the project: #
  * `ant clean compile package`

# Running the tests: #
  * `ant runTests`

  * `ant runTest -Dtest=algorithms.imageProcessing.CornersOfLabTest`

Using "runTestWithAspects" runs a specific test in detailed debug mode:
  * `ant runTestWithAspects -Dtest=algorithms.imageProcessing.CornersOfLabTest`

> (The display of debug data won't time out, so you'll have to cntrl-c to end the test.  It's meant to make it easier to examine intermediate files.)

# Misc Other: #

There are code quality scripts and programs that can be run to look at code coverage, dependency, and control flow.  These are currently tailored only for posix systems (unix and linux).   The code coverage script looks for a few files in a maven repository, so it isn't tidied up for general use yet.

The dtrace code is in the tests directory and is called Run\_dtrace\_Test.  It's currently not enabled by default, so edit that.  The user may need to set up ACL for superuser to run it.  The current Run\_dtrace\_Test is tracking class visits for a look at dependencies, but the method level visits flow script can be found in the dtrace directory instead (and edit Run\_dtrace\_Test to use that).