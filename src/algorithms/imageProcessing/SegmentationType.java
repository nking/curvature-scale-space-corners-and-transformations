package algorithms.imageProcessing;

/**
<pre>
    GREYSCALE_KMPP = KMPP for k=2;
    COLOR_POLARCIEXY_ADAPT = color to greyscale using polar theta of CIEXY 
        scale space followed by adaptive mean thresholding (scl=4?);
    COLOR_POLARCIEXY = color to greyscale using polar theta of CIEXY;
    COLOR_POLARCIEXY_LARGE = COLOR_POLARCIEXY with a larger range of blob sizes, 
        5000 to 100000;
    ADAPTIVE_MEAN uses adaptive mean thresholding w h=2
</pre>
@author nichole
 */
public enum SegmentationType {

    GREYSCALE_KMPP,
    COLOR_POLARCIEXY_LARGE,
    COLOR_POLARCIEXY_ADAPT,
    COLOR_POLARCIEXY,
    ADAPTIVE_MEAN;

}
