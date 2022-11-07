test data from:
http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
now at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
"Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
by Jean-Yves Bouguet
Camera Calibration Toolbox for Matlab

the web page is saved in pdf format in miscNotes/bouguetj_5th_calibration_example.pdf

There are 14 images from the left and right cameras each, but only the first
pair are in this directory to use in tests.

There is a python load script that can be used to read the .mat files
at https://github.com/malcolmreynolds/calib_bouguet
The output to that is here as left_summary.txt and right_summary.txt

Each checkerboard square is 30mm in WCS metric.
There are 7 rows in the checkerboard and 9 columns.
The checkerboard width is 305 pixels in image0.jpg,
and the height is 240 pixels.
The image width is 640 pixels.
The image height is 480 pixels.
So the image width is about 640 pixels * ((30 mm * 9)/305 pixels) ~ 570 mm.
The image height is about 480 pixels * ((30 mm * 7)/240 pixels) ~ 420 mm.

The checkerboard image coordinates and WCS coordinates are not present in the dataset.
One can extract the image coordinates using code in this project, which is at a rougher
stage (cannot mask out other image features for example),
or use Matlab or OpenCV or python's opencv, etc.

The rough range of WCS coordinates that can be checked with results from Triangulate are
in the pdf spanshot of the calibration toolbox in miscNotes directory.

checkerboard position for image 1: distances (Z coords) look to be 450:500
checkerboard postiion for image 2: distances (Z coords) look to be 300:450

     relative rotation between right and left cameras:
       om = 0.00611, 0.00489, -0.00359
        t = -99.84929, 0.82221, 0.43647
     after optimization:
       om = 0.00669, 0.00452, -0.00350 +- 0.0027, 0.00308, 0.00029
        t = -99.80198, 1.12443, 0.05041 +- 0.142, 0.11352, 0.49773

