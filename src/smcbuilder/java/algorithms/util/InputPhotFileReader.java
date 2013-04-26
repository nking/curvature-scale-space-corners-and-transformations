package algorithms.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author nichole
 */
public class InputPhotFileReader implements IInputFileReader{

    protected float[] x = null;
    protected float[] y = null;
    protected float[] xErrors = null;
    protected float[] yErrors = null;
    protected int nPoints = 0;

    protected Logger log = null;

    protected String filePath = null;

    public InputPhotFileReader(){
    }

    protected void initArrays() {
        int n = 100;
        x = new float[n];
        y = new float[n];
        xErrors = new float[n];
        yErrors = new float[n];
    }

    protected String getFilePath() {
        return filePath;
    }

    public void read(String pathToFile) throws IOException {

        this.filePath = pathToFile;

        log = Logger.getLogger(this.getClass().getName());

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

            Pattern p = Pattern.compile("\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s*");

            while (line != null) {

                Matcher m = p.matcher(line);

                expandArraysIfNecessary();

                if (m.matches()) {

                    float xp = Float.valueOf(m.group(1));
                    float yp = Float.valueOf(m.group(2));
                    float mI = Float.valueOf(m.group(3));

                    boolean t1 = (mI < 20.5f);
                    boolean t2 = true;
                    if (filePath.contains("smc111.4")) {
                        t2 = (xp < 2050);
                    } else if (filePath.contains("smc115.5")) {
                        t2 = (xp < 2100);
                    } else if (filePath.contains("smc110.3")) {
                        t2 = (yp < 4050);
                    } else if (filePath.contains("smc125.4")) {
                        t2 = (xp < 2050);
                    } else if (filePath.contains("smc100.5")) {
                        t2 = (yp < 4050);
                    } else if (filePath.contains("smc100.1")) {
                        t2 = (yp < 4050);
                        t2 = (xp < 2050);
                    }

                    if (t1 && t2){

                        x[nPoints] = xp;
                        xErrors[nPoints] = 0.3f;// assuming centering errors are fraction of a pixel

                        y[nPoints] = yp;
                        xErrors[nPoints] = 0.3f;// assuming centering errors

                        nPoints++;
                    }

                } else {
                    log.info("Did not match=>" + line + "<=");
                }

                line = in.readLine();
            }

            condenseArrays();

        } finally {
            if (reader != null) {
                reader.close();
            }

            log.info("finished loading file of " + nPoints + " points");
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

    public float[] getX() {
        return x;
    }
    public float[] getY() {
        return y;
    }
    public float[] getXErrors() {
        return xErrors;
    }
    public float[] getYErrors() {
        return yErrors;
    }
}
