package algorithms.curves;

/**
 *
 * @author nichole
 */
public class CurveMisc {

    /**
     *
     * @param x
     * @param y
     * @param xyIndex
     * @param isStepFunction
     * @param xScaleFactor
     * @param yScaleFactor
     * @return
     */
    public static float calculateArea(float[] x, float[] y, int xyIndex,
        boolean isStepFunction, float xScaleFactor, float yScaleFactor) {

        if (x.length == 0) {

            return 0.0f;

        } else if (x.length == 1) {

            float w = xScaleFactor * x[0];

            float h = yScaleFactor * y[0];

            return w * h;

        }

        float x0, xmid, x1, y0, ymid, y1;

        /*                              *
         *                          .      .
         *   *         *         *  ........  *
         *     .     .              .      .
         *     ...*...              .      .
         *     .     .
         *     .     .
         *
         */

        xmid = x[xyIndex];
        ymid = y[xyIndex];

        if (xyIndex == 0) {

            float xDelta = ((x[xyIndex + 1] - x[xyIndex])/2.0f);

            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;

            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                float yDelta = (y[1] - y[0])/2.0f;
                y0 = y[0] - yDelta;
                y1 = y[0] + yDelta;
            }

        } else if (xyIndex == (y.length - 1)) {

            float xDelta = ((x[xyIndex] - x[xyIndex - 1])/2.0f);

            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;

            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                float yDelta = (y[xyIndex - 1] - y[xyIndex])/2.0f;
                y0 = ymid + yDelta;
                y1 = ymid - yDelta;
            }

        } else {

            x0 = ((x[xyIndex] + x[xyIndex - 1])/2.0f);
            x1 = ((x[xyIndex + 1] + x[xyIndex])/2.0f);

            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                y0 = (y[xyIndex - 1] + y[xyIndex])/2.0f;
                y1 = (y[xyIndex] + y[xyIndex + 1])/2.0f;
            }
        }

        float areaBase = xScaleFactor*(x1 - x0)*yScaleFactor*(ymid);

        if (isStepFunction) {
            return areaBase;
        }

        float areaTop0 = 0.5f*xScaleFactor*(xmid - x0)*yScaleFactor*(y0 - ymid);
        float areaTop1 = 0.5f*xScaleFactor*(x1 - xmid)*yScaleFactor*(y1 - ymid);

        float area = areaTop0 + areaTop1 + areaBase;

        return area;
    }
}
