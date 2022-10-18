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
