package algorithms.imageProcessing;

import java.util.List;

/**
 * Looks through each curve to see if the start and ends meet within an allowed
 * gap, and if so, replaces the PairIntArray with a PairIntArrayWithColor and
 * color='1'.  Note that gaps are also filled in the replaced edge.
 * 
 * The code also looks for T-junctions and returns a list of points for those.
 * 
 * @author nichole
 */
public class ClosedCurveAndJunctionFinder {
    
    private int allowedGapPlusOne = 3;
    
    /*
       @
      @ @
     @   @
    @__*__@  radius of 2, npoints = 12
     @   @
      @ @
       @
    */
    private int minEdgePoints = 12;
    
    public ClosedCurveAndJunctionFinder() {
        
    }
    
    public void findClosedCurves(List<PairIntArray> edges) {
        
        for (int i = 0; i < edges.size(); i++) {
            
            PairIntArray edge = edges.get(i);
            
            if (edge.getN() < minEdgePoints) {
                continue;
            }
            
            int x0 = edge.getX(0);
            int y0 = edge.getY(0);
            
            int xn = edge.getX(edge.getN() - 1);
            int yn = edge.getY(edge.getN() - 1);
            
            int diffX = x0 - xn;
            if (diffX < 0) {
                diffX *= -1;
            }
            if (diffX > allowedGapPlusOne) {
                continue;
            }
            int diffY = y0 - yn;
            if (diffY < 0) {
                diffY *= -1;
            }
            if (diffY > allowedGapPlusOne) {
                continue;
            }
            
            //  0  1  2  3
            
            if ((diffX > 1) || diffY > 1) {
                
                /* fill in the gaps
                    (1)           (2)          (3)
                              |       xn,yn | x0,y0
                x0,y0   xn,yn | x0,y0       |       xn,yn
                --------------------------------------------
                              |       x0,y0 | xn,yn
                xn,yn   x0,y0 | xn,yn       |       x0,y0
                   (4)           (5)           (6)
                --------------------------------------------
                  xn,yn       |  x0,y0
                  x0,y0       |  xn,yn
                   (7)           (8)
                */
                
                //append points starting with those closer to (xn, yn)
                
                if (x0 < xn) {
                    if (y0 == yn) {
                        //(1)
                        for (int x = (xn - 1); x > x0; x--) {
                            edge.add(x, y0);
                        }
                    } else if (y0 < yn) {
                        //(2)
                        for (int x = (xn - 1); x > x0; x--) {
                            for (int y = (yn - 1); y > y0; y--) {
                                edge.add(x, y);
                            }
                        }
                    } else {
                        //(3)
                        for (int x = (xn - 1); x > x0; x--) {
                            for (int y = yn + 1; y < y0; y++) {
                                edge.add(x, y);
                            }
                        }
                    }
                } else if (x0 > xn) {
                    if (y0 == yn) {
                        //(4)
                        for (int x = (xn + 1); x < x0; x++) {
                            edge.add(x, y0);
                        }
                    } else if (yn < y0) {
                        //(5)
                        for (int x = (xn + 1); x < x0; x++) {
                            for (int y = (yn + 1); y < y0; y++) {
                                edge.add(x, y);
                            }
                        }
                    } else {
                        //(6)
                        for (int x = (xn + 1); x < x0; x++) {
                            for (int y = (yn - 1); y > y0; y--) {
                                edge.add(x, y);
                            }
                        }
                    }
                } else {
                    // x0 == xn
                    if (yn > y0) {
                        // (7)
                        for (int y = (yn - 1); y > y0; y--) {
                            edge.add(x0, y);
                        }
                    } else {
                        // (8)
                        for (int y = (yn + 1); y < y0; y++) {
                            edge.add(x0, y);
                        }
                    }
                }
            }
            
            PairIntArrayWithColor closedCurve = new PairIntArrayWithColor(edge);
            closedCurve.setColor(1);
            
            edges.set(i, closedCurve);
        }
    }
    
    /**
     * not yet implemented.  returns an empty point list.
     * 
     * @param edges
     * @return 
     */
    public PairIntArray findTJunctions(List<PairIntArray> edges) {
        
        PairIntArray junctions = new PairIntArray();
        
        return junctions;
    }
}
