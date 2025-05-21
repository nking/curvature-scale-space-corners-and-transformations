package algorithms.dimensionReduction;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.ResourceFinder;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;

public class KaggleDatasets {

    public static String EYESDIR = "EYES";
    protected final static String sep = System.getProperty("file.separator");

    /**
     * @misc{pavel biza_2021,
     *                 title={Female and male eyes},
     *                 url={https://www.kaggle.com/ds/1438879},
     *             DOI={10.34740/KAGGLE/DS/1438879},
     *                     publisher={Kaggle},
     *                     author={Pavel Biza},
     *                     year={2021}
     *
     *  The dataset has 6324 male eye images and 5202 female eye images,
     *  but this project uses less than 100 of each and reduces them in size.
     * @return
     */
    public double[][] getMaleEyes(int width) throws IOException {

        String dirPath = ResourceFinder.findTestResourcesDirectory() + sep + EYESDIR
                + sep + "male" + sep;
        String path;
        File f;
        ImageExt image;
        int i;
        for (i = 0; i < 100; ++i) {
            path = dirPath + Integer.toString(i) + ".jpg";
            f = new File(path);
            if (!f.exists()) {
                continue;
                //"could not find file at " + path;
            }
            image = ImageIOHelper.readImageExt(path);
        }
        throw new UnsupportedEncodingException("not yet implemented");
    }
}
