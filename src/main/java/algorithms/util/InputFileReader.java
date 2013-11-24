package algorithms.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class InputFileReader implements IInputFileReader {

    protected float[] x = null;
    protected float[] y = null;
    protected float[] xErrors = null;
    protected float[] yErrors = null;
    protected int nPoints = 0;

    protected Logger log = null;

    protected String filePath = null;

    public InputFileReader(){
        log = Logger.getLogger(this.getClass().getName());
    }

    protected void initArrays() {
        int n = 100;
        x = new float[n];
        y = new float[n];
        xErrors = new float[n];
        yErrors = new float[n];
    }

    public void read(String pathToFile) throws IOException {

        this.filePath = pathToFile;

        initArrays();

        FileReader reader = null;
        BufferedReader in = null;

        try {
            File fl = new File(filePath);

            if (!fl.exists()) {
                throw new IOException("file " + filePath + " does not exist");
            }

            reader = new FileReader(fl);
            in = new BufferedReader(reader);

            String line = in.readLine();

            log.info("loading " + filePath);

            while (line != null) {

                String[] values = line.split("\t");
                if (values == null || values.length != 4) {
                    throw new IOException("file data has to be tab separated, and 4 columns of: x y xError yError");
                }

                expandArraysIfNecessary();

                for (int i = 0; i < values.length; i++) {

                    Float value = Float.valueOf(values[i]);

                    switch (i) {
                        case 0:
                            x[nPoints] = value;
                            break;
                        case 1:
                            y[nPoints] = value;
                            break;
                        case 2:
                            xErrors[nPoints] = value;
                            break;
                        case 3:
                            yErrors[nPoints] = value;
                            break;
                        default:
                            break;
                    }
                }

                nPoints++;

                line = in.readLine();
            }

        } finally {
            if (reader != null) {
                reader.close();
            }
            if (in != null) {
                in.close();
            }
            condenseArrays();
            log.info("finished loading file");
        }
    }

    protected void expandArraysIfNecessary() {
        if (nPoints > (x.length - 2)) {
            int expand = (int)(1.3f*x.length);
            x = Arrays.copyOf(x, expand);
            y = Arrays.copyOf(y, expand);
            xErrors = Arrays.copyOf(xErrors, expand);
            yErrors = Arrays.copyOf(yErrors, expand);
        }
    }
    protected void condenseArrays() {
        x = Arrays.copyOf(x, nPoints);
        y = Arrays.copyOf(y, nPoints);
        xErrors = Arrays.copyOf(xErrors, nPoints);
        yErrors = Arrays.copyOf(yErrors, nPoints);
    }

    @Override
    public float[] getX() {
        return x;
    }
    @Override
    public float[] getY() {
        return y;
    }
    @Override
    public float[] getXErrors() {
        return xErrors;
    }
    @Override
    public float[] getYErrors() {
        return yErrors;
    }
}
