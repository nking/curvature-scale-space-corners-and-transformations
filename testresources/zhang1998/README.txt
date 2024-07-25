downloaded from
https://www.microsoft.com/en-us/research/project/a-flexible-new-technique-for-camera-calibration-2/

Here's a description of the data:
The model plane contains a pattern of 8×8 squares, so there are 256 corners. 
The size of the pattern is 17cm x 17cm. 
The 2D coordinates (in inches) of these points are in Model.txt.

We have taken five an off-the-shelf PULNiX CCD camera with 6 mm lens. 
The image resolution is 640×480. The five images are available here: 
Image 1, Image 2, Image 3, Image 4 and Image 5. 
The first two are shown below. We can observe a significant lens distortion in the images.

And here is what the calibration tells us about the camera: 
The pixel is square (aspect ratio = 1); 
the focal length = 832.5 pixels; 
the image center is at (303.959, 206.585); 
there is a significant radial distortion: k1 = -0.228601, k2 = 0.190353. 
The complete calibration result is available here. 
(The format of the calibration file [Calib_Results.mat] is: 
   a, c, b, u0, v0, k1, k2, then the rotation matrix and translation vector for the 
first image, the rotation matrix and translation vector for the second image, etc.)
where a and b the scale factors in image u and v axes
c the parameter describing the skewness of the two image axes.
(u0, v0) the coordinates of the principal point,
    for radial distortion, they use the k1*r + k2*r^2 pattern.

Model.txt is 64 lines of 8 real numbers
     they hold the WCS (X,Y) coordinates of the checkerboard corners, w/ assumption that Z=0
data?.txt are each 64 lines of 8 real numbers.
     they hold the image coordinates of the checkerboard corners

