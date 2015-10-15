package algorithms.imageProcessing;

/**
<pre>
    BINARY = KMPP for k=2 followed by adaptive mean thresholding 
        (scl=20/binFactor);
    GREYSCALE_KMPP = KMPP for k=2;
    COLOR_POLARCIEXY_ADAPT = color to greyscale using polar theta of CIEXY 
        scale space followed by adaptive mean thresholding (scl=4?);
    COLOR_POLARCIEXY = color to greyscale using polar theta of CIEXY.
</pre>
@author nichole
 */
public enum SegmentationType {

    BINARY,
    GREYSCALE_KMPP,
    COLOR_POLARCIEXY_ADAPT,
    COLOR_POLARCIEXY;

}
