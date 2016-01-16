package algorithms.imageProcessing;

/**
<pre>
    GREYSCALE_WAVELET = greyscale segmentation based upon coarsest coefficient
                        of a trous wavelet transform;
    GREYSCALE_CANNY = greyscale segmentation based upon Canny edge detector;
    GREYSCALE_HIST = greyscale segmentation based on histogram with a black and
                     white mask
    GREYSCALE_KMPP = KMPP for k=2;
    DT_CLUSTERING = density based clustering filtered for top k results
    COLOR_POLARCIEXY_ADAPT = color to greyscale using polar theta of CIEXY 
        scale space followed by adaptive mean thresholding (scl=4?);
    COLOR_POLARCIEXY = color to greyscale using polar theta of CIEXY;
    ADAPTIVE_MEAN uses adaptive mean thresholding w h=2
    NONE = no segmentation is performed allowing special handling of the points.
</pre>
@author nichole
 */
public enum SegmentationType {
    GREYSCALE_WAVELET, GREYSCALE_CANNY,
    GREYSCALE_HIST,
    GREYSCALE_KMPP,
    DT_CLUSTERING,
    COLOR_POLARCIEXY_ADAPT,
    COLOR_POLARCIEXY,
    ADAPTIVE_MEAN,
    NONE;
}
