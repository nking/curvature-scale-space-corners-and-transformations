package algorithms.compGeometry.voronoi;

import algorithms.util.PairFloat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;

/**
 * A voronoi diagram is a partitioning of points
 * into cells defined by seed centers.  The cells
 * edges are defined as perpendicular bisectors
 * between 2 seeds, which results in the cells 
 * being regions which are closer to the seed
 * point than to any other seed point.
 * 
 * The sweep line algorithm runtime complexity
 * is  O(n log n).
 * 
  adapted from the following codes:
 
   a port to java of the C++ port of
 * Steven Fortune's original c code.
 * The C++ port by Shane O' Sullivan is at 
 * http://skynet.ie/~sos/mapviewer/voronoi.php
 * and uses a liberal AT&T license.
 * Then 2 ports to java.  
 * https://sourceforge.net/projects/simplevoronoi/postdownload?source=dlp
 * All have liberal licenses listed here.
 
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 
 * This code was originally written by Stephan Fortune in C code.  I, Shane O'Sullivan,
 * have since modified it, encapsulating it in a C++ class and, fixing memory leaks and
 * adding accessors to the Voronoi Edges.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
  
 * Java Version by Zhenyu Pan
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 *
 * 
 * in this project:
 * changes to the above java port are essentially 
 * to the comparator and sorting, and an extra check for
 * zero length edges at the end of the method clipLine.
 * Format changes such as camel case are present 
 * here.  Also moved the data classes to inner classes
 * here, and am using specialized classes that already
 * existed in this project such as the PairFloat for 
 * coordinates.
 */
public class VoronoiFortunesSweep {
    
    private static int LE = 0;
    private static int RE = 1;

    private HalfEdge ELleftend = null;
    private HalfEdge ELrightend = null;
    private int ELhashsize;
    private float xmin, xmax, ymin, ymax, deltax, deltay;

    private Site[] sites = null;
    private int nSites = 0;
    private int siteIdx = 0;
    private int sqrtNSites;
    private int nVertices = 0;
    private int nEdges = 0;
    
    private Site bottomsite;
	
    private int	PQcount;
    private int	PQmin;
    private int PQhashsize;
    private HalfEdge[] PQhash;

    private float borderMinX, borderMaxX, borderMinY, borderMaxY;

    private LinkedList<GraphEdge> allEdges = null;

    private float minDistanceBetweenSites = 0;
    
    private HalfEdge[] ELhash = null;
    
    public VoronoiFortunesSweep() {
        allEdges = new LinkedList<GraphEdge>();
    }
    
    public LinkedList<GraphEdge> getAllEdges() {
        return allEdges;
    }
    
    public boolean generateVoronoi(
        float[] xValues, float[] yValues, float minX, float maxX, 
        float minY, float maxY, float minDist) {
        
        if (xValues.length != yValues.length) {
            throw new IllegalArgumentException(
            "xValues and yValues must be same length");
        }

        minDistanceBetweenSites = minDist;
    
        siteIdx = 0;
        nSites = xValues.length;
        nVertices = 0;
        nEdges = 0;
        double sn = (double) nSites + 4;
        sqrtNSites = (int) Math.sqrt(sn);

        sites = new Site[nSites];

        xmin = xValues[0];
        ymin = yValues[0];
        xmax = xValues[0];
        ymax = yValues[0];
        
        int i;
        for(i = 0; i < nSites; i++) {
            sites[i] = new Site();
            sites[i].coord = new PairFloat(xValues[i], yValues[i]);
            sites[i].sitenbr = i;
            sites[i].refcnt = 0;

            if (xValues[i] < xmin) {
                xmin = xValues[i];
            } else if (xValues[i] > xmax) {
                xmax = xValues[i];
            }
            
            if (yValues[i] < ymin) {
                ymin = yValues[i];
            } else if (yValues[i] > ymax) {
                ymax = yValues[i];
            }
            
            //printf("\n%f %f\n",xValues[i],yValues[i]);
        }
        
        // N log_2 N
        qsort(sites);
        
        float temp = 0;
        if (minX > maxX) {
            temp = minX;
            minX = maxX;
            maxX = temp;
        }
        if (minY > maxY) {
            temp = minY;
            minY = maxY;
            maxY = temp;
        }
        borderMinX = minX;
        borderMinY = minY;
        borderMaxX = maxX;
        borderMaxY = maxY;
        
        deltay = ymax - ymin;
        deltax = xmax - xmin;
        
        return voronoi();        
    }
    
    /*
       implicit parameters: nsites, sqrt_nsites, 
       xmin, xmax, ymin, ymax, deltax,
       deltay (can all be estimates). 
       Performance suffers if they are wrong;
       better to make nsites, deltax, 
       and deltay too big than too small. (?)
     */
    private boolean voronoi() {
        
        Site newsite, bot, top, temp, p;
        Site v;
        PairFloat newintstar = null;
        int pm;
        HalfEdge lbnd, rbnd, llbnd, rrbnd, bisector;
        Edge e;

        PQInitialize();
        ELInitialize();

        bottomsite = nextOne();
        newsite = nextOne();
        
        while (true) {
            
            if (!PQEmpty()) {
                newintstar = PQMin();
            }
            
            // if the lowest site has a smaller y value than 
            // the lowest vector intersection,
            // process the site otherwise process the vector 
            // intersection

            if (newsite != null
            && (PQEmpty() || 
                newsite.coord.getY() < newintstar.getY() || 
                (newsite.coord.getY() == newintstar.getY() && 
                newsite.coord.getX() < newintstar.getX()))) {
                
                // new site is smallest -this is a site event
                         
                // get the first HalfEdge to the LEFT of the new site
                lbnd = ELLeftBnd((newsite.coord));
                
                // get the first HalfEdge to the RIGHT of the new site
                rbnd = lbnd.ELright;
                
                // if this halfedge has no edge,bot =bottom site 
                bot = rightReg(lbnd);
                
                // create a new edge that bisects
                e = bisect(bot, newsite);

                // create a new HalfEdge, setting its ELpm field to 0
                bisector = HECreate(e, LE);
                
                // insert this new bisector edge between the left 
                // and right vectors in a linked list
                ELInsert(lbnd, bisector);

                // if the new bisector intersects with the left edge,
                // remove the left edge's vertex, and put in the new one
                if ((p = intersect(lbnd, bisector)) != null) {
                    PQDelete(lbnd);
                    PQInsert(lbnd, p, dist(p, newsite));
                }
                lbnd = bisector;
                
                // create a new HalfEdge, setting its ELpm field to 1
                bisector = HECreate(e, RE);
                
                // insert the new HE to the right of the original bisector
                // earlier in the IF stmt
                ELInsert(lbnd, bisector);

                // if this new bisector intersects with the new HalfEdge
                if ((p = intersect(bisector, rbnd)) != null) {
                    // push the HE into the ordered linked list of vertices
                    PQInsert(bisector, p, dist(p, newsite));
                }
                
                newsite = nextOne();
                
            } else if (!PQEmpty()) {
            
                // intersection is smallest - this is a vector event
            
                // pop the HalfEdge with the lowest vector off the ordered list
                // of vectors
                lbnd = PQExtractMin();
                
                // get the HalfEdge to the left of the above HE
                llbnd = lbnd.ELleft;
                
                // get the HalfEdge to the right of the above HE
                rbnd = lbnd.ELright;
                
                // get the HalfEdge to the right of the HE to the right of the
                // lowest HE
                rrbnd = rbnd.ELright;
                
                // get the Site to the left of the left HE which it bisects
                bot = leftReg(lbnd);
                
                // get the Site to the right of the right HE which it bisects
                top = rightReg(rbnd);

                // get the vertex that caused this event
                v = lbnd.vertex;
                
                // set the vertex number - couldn't do this
                makeVertex(v);
                
                // earlier since we didn't know when it would be processed
                endPoint(lbnd.ELedge, lbnd.ELpm, v);
                
                // set the endpoint of
                // the left HalfEdge to be this vector
                endPoint(rbnd.ELedge, rbnd.ELpm, v);
                
                // set the endpoint of the right HalfEdge to
                // be this vector
                ELDelete(lbnd); // mark the lowest HE for
                
                // deletion - can't delete yet because there might be pointers
                // to it in Hash Map
                PQDelete(rbnd);
                
                // remove all vertex events to do with the right HE
                ELDelete(rbnd); 
                
                // mark the right HE for
                // deletion - can't delete yet because there might be pointers
                // to it in Hash Map
                
                // set the pm variable to zero
                pm = LE;

                if (bot.coord.getY() > top.coord.getY()) {
                    // if the site to the left of the event is 
                    // higher than the Site
                    // to the right of it, then swap them and 
                    // set the 'pm' variable to 1
                    temp = bot;
                    bot = top;
                    top = temp;
                    pm = RE;
                }
                
                // create an Edge (or line)
                // that is between the two Sites. This creates the formula of
                // the line, and assigns a line number to it
                e = bisect(bot, top); 
                
                // create a HE from the Edge 'e',
                // and make it point to that edge
                // with its ELedge field
                bisector = HECreate(e, pm); 
                
                // insert the new bisector to the
                // right of the left HE
                ELInsert(llbnd, bisector);               
                
                // set one endpoint to the new edge
                // to be the vector point 'v'.
                // If the site to the left of this bisector is higher than the
                // right Site, then this endpoint
                // is put in position 0; otherwise in pos 1
                endPoint(e, RE - pm, v);
                
                // if left HE and the new bisector intersect, then delete
                // the left HE, and reinsert it
                if ((p = intersect(llbnd, bisector)) != null) {
                    PQDelete(llbnd);
                    PQInsert(llbnd, p, dist(p, bot));
                }

                // if right HE and the new bisector intersect, then
                // reinsert it
                if ((p = intersect(bisector, rrbnd)) != null) {
                    PQInsert(bisector, p, dist(p, bot));
                }
            } else {
                break;
            }
        }

        for (lbnd = ELleftend.ELright; lbnd != ELrightend; 
            lbnd = lbnd.ELright) {
            
            e = lbnd.ELedge;
            clipLine(e);
        }

        return true;
    }

    private class SiteComp implements Comparator<Site> {
        @Override
        public int compare(Site o1, Site o2) {
            if (o1.coord.getY() < o2.coord.getY()) {
                return -1;
            } else if (o1.coord.getY() > o2.coord.getY()) {
                return 1;
            } else if (o1.coord.getX() < o2.coord.getX()) {
                return -1;
            } else if (o1.coord.getX() > o2.coord.getX()) {
                return 1;
            }
            return 0;
        }
    }

    private void qsort(Site[] sites) {
        SiteComp comp = new SiteComp();
        Arrays.sort(sites, comp);
    }
    
    public class Site {
        PairFloat coord;
        int	sitenbr;
	    int	refcnt;
        public PairFloat getCoord() {
            return coord;
        }
    }

    private class Edge {
        float a = 0;
        float b = 0;
        float c = 0;
	    Site[] ep = new Site[2];
	    Site[] reg = new Site[2];
	    int edgenbr;
    }
    
    public class GraphEdge {
        float x1, y1, x2, y2;
	    public int site1;
        public int site2;
    }
     
    private class HalfEdge {
	    HalfEdge ELleft = null; 
        HalfEdge ELright = null;
	    Edge ELedge = null;
        boolean deleted;
	    int ELpm;
	    Site vertex = null;
	    float ystar;
	    HalfEdge PQnext = null;
    }

    // return a single in-storage site
    private Site nextOne() {
        Site s;
        if (siteIdx < nSites) {
            s = sites[siteIdx];
            siteIdx += 1;
            return (s);
        } else {
            return (null);
        }
    }

    private Edge bisect(Site s1, Site s2) {
        float dx, dy, adx, ady;
        Edge newEdge;

        newEdge = new Edge();

        // store the sites that this edge is bisecting
        newEdge.reg[0] = s1;
        newEdge.reg[1] = s2;
        // to begin with, there are no endpoints on the bisector - it goes to
        // infinity
        newEdge.ep[0] = null;
        newEdge.ep[1] = null;

        // get the difference in x dist between the sites
        dx = s2.coord.getX() - s1.coord.getX();
        dy = s2.coord.getY() - s1.coord.getY();
        
        // make sure that the difference in positive
        adx = dx > 0 ? dx : -dx;
        ady = dy > 0 ? dy : -dy;
        
        // get the slope of the line
        newEdge.c = 
            (s1.coord.getX() * dx + s1.coord.getY() 
                * dy + (dx * dx + dy
                * dy) * 0.5f);

        if (adx > ady) {
            // set formula of line, with x fixed to 1
            newEdge.a = 1.0f;
            newEdge.b = dy / dx;
            newEdge.c /= dx;
        } else {
            // set formula of line, with y fixed to 1
            newEdge.b = 1.0f;
            newEdge.a = dx / dy;
            newEdge.c /= dy;
        }

        newEdge.edgenbr = nEdges;

        nEdges += 1;
        
        return newEdge;
    }
    
    public Site[] getSites() {
        return sites;
    }

    private void makeVertex(Site v) {
        v.sitenbr = nVertices;
        nVertices += 1;
    }

    private void PQInitialize() {
        PQcount = 0;
        PQmin = 0;
        PQhashsize = 4 * sqrtNSites;
        PQhash = new HalfEdge[PQhashsize];

        for (int i = 0; i < PQhashsize; i += 1) {
            PQhash[i] = new HalfEdge();
        }
    }

    private int PQBucket(HalfEdge he) {
        int bucket;

        bucket = (int) ((he.ystar - ymin) / deltay * PQhashsize);
        if (bucket < 0) {
            bucket = 0;
        }
        if (bucket >= PQhashsize) {
            bucket = PQhashsize - 1;
        }
        if (bucket < PQmin) {
            PQmin = bucket;
        }
        return bucket;
    }

    // push the HalfEdge into the ordered linked list of vertices
    private void PQInsert(HalfEdge he, Site v, float offset) {
        
        HalfEdge last;
        HalfEdge next;

        he.vertex = v;
        he.ystar = (v.coord.getY() + offset);
        last = PQhash[PQBucket(he)];
  
        // NOTE: for the cases where the number
        //  of items in the linked list is getting
        //  large enough that this is not approx O(1),
        //  could initialize a specialized comparator
        //  with v.coord.getX()
        //  and store the nodes as an array in this bin
        //  instead of a linked list
        //  and use binarySearch and the comparator with
        //  a small scan to find the last node matching
        // the criteria.
        while ((next = last.PQnext) != null
            && (he.ystar > next.ystar
               || (he.ystar == next.ystar
               && v.coord.getX() > next.vertex.coord.getX()))) {
            last = next;
        }
        he.PQnext = last.PQnext;
        last.PQnext = he;
        PQcount += 1;
    }

    // remove the HalfEdge from the list of vertices
    private void PQDelete(HalfEdge he) {
        
        HalfEdge last;

        if (he.vertex != null) {
            last = PQhash[PQBucket(he)];
            while (last.PQnext != he) {
                last = last.PQnext;
            }

            last.PQnext = he.PQnext;
            PQcount -= 1;
            he.vertex = null;
        }
    }

    private boolean PQEmpty() {
        return (PQcount == 0);
    }

    private PairFloat PQMin() {
        
        while (PQhash[PQmin].PQnext == null) {
            PQmin += 1;
        }
        
        PairFloat answer = new PairFloat(
            PQhash[PQmin].PQnext.vertex.coord.getX(),
            PQhash[PQmin].PQnext.ystar);
        
        return answer;
    }

    private HalfEdge PQExtractMin() {
        HalfEdge curr;

        curr = PQhash[PQmin].PQnext;
        PQhash[PQmin].PQnext = curr.PQnext;
        PQcount -= 1;
        return curr;
    }

    private HalfEdge HECreate(Edge e, int pm) {
        HalfEdge answer;
        answer = new HalfEdge();
        answer.ELedge = e;
        answer.ELpm = pm;
        answer.PQnext = null;
        answer.vertex = null;
        return answer;
    }

    private void ELInitialize() {
        
        ELhashsize = 2 * sqrtNSites;
        ELhash = new HalfEdge[ELhashsize];

        ELleftend = HECreate(null, 0);
        ELrightend = HECreate(null, 0);
        
        ELleftend.ELleft = null;
        ELleftend.ELright = ELrightend;
        ELrightend.ELleft = ELleftend;
        ELrightend.ELright = null;
        
        ELhash[0] = ELleftend;
        ELhash[ELhashsize - 1] = ELrightend;
    }

    private Site leftReg(HalfEdge he) {
        if (he.ELedge == null) {
            return bottomsite;
        }
        return he.ELpm == LE ? he.ELedge.reg[LE] : 
            he.ELedge.reg[RE];
    }

    private void ELInsert(HalfEdge lb, HalfEdge newHe) {
        newHe.ELleft = lb;
        newHe.ELright = lb.ELright;
        (lb.ELright).ELleft = newHe;
        lb.ELright = newHe;
    }

    /*
     * This delete routine can't reclaim node, since pointers from hash table
     * may be present.
     */
    private void ELDelete(HalfEdge he) {
        (he.ELleft).ELright = he.ELright;
        (he.ELright).ELleft = he.ELleft;
        he.deleted = true;
    }

    // Get entry from hash table, pruning any deleted nodes
    private HalfEdge ELGetHash(int b) {
        HalfEdge he;

        if (b < 0 || b >= ELhashsize) {
            return null;
        }
        he = ELhash[b];
        if (he == null || !he.deleted) {
            return he;
        }

        // Hash table points to deleted half edge. Patch 
        // as necessary.
        ELhash[b] = null;
        
        return null;
    }

    private HalfEdge ELLeftBnd(PairFloat p) {
        int i, bucket;
        HalfEdge he;

        // Use hash table to get close to desired halfedge
        // use the hash function to find the place in the hash map that this
        // HalfEdge should be
        bucket = (int) ((p.getX() - xmin) / deltax * ELhashsize);

        // make sure that the bucket position in within the range of the hash
        // array
        if (bucket < 0) {
            bucket = 0;
        }
        if (bucket >= ELhashsize) {
            bucket = ELhashsize - 1;
        }

        he = ELGetHash(bucket);
        if (he == null) {
            
            //TODO: could improve this with a datastructure
            // that has predecessor and successor
            
            // if the HE isn't found, search backwards and forwards in the hash map
            // for the first non-null entry
        
            for (i = 1; i < ELhashsize; i += 1) {
                if ((he = ELGetHash(bucket - i)) != null) {
                    break;
                }
                if ((he = ELGetHash(bucket + i)) != null) {
                    break;
                }
            }
        }
        
        // Now search linear list of halfedges for the correct one
        if (he == ELleftend || 
        (he != ELrightend && rightOf(he, p))) {
            // keep going right on the list until either the end is reached, or
            // you find the 1st edge which the point isn't to the right of
            do {
                he = he.ELright;
            } while (he != ELrightend && rightOf(he, p));
            he = he.ELleft;
        } else {
            // if the point is to the left of the HalfEdge, then search left for
            // the HE just to the left of the point
        
            do {
                he = he.ELleft;
            } while (he != ELleftend && !rightOf(he, p));
        }

        // Update hash table and reference counts
        if (bucket > 0 && bucket < ELhashsize - 1) {
            ELhash[bucket] = he;
        }
        return he;
    }

    private void pushGraphEdge(Site leftSite, Site rightSite, 
        float x1, float y1, float x2, float y2) {
        
        GraphEdge newEdge = new GraphEdge();
        allEdges.add(newEdge);
        newEdge.x1 = x1;
        newEdge.y1 = y1;
        newEdge.x2 = x2;
        newEdge.y2 = y2;

        newEdge.site1 = leftSite.sitenbr;
        newEdge.site2 = rightSite.sitenbr;
    }

    private void clipLine(Edge e) {
        float pxmin, pxmax, pymin, pymax;
        Site s1, s2;
        float x1 = 0, x2 = 0, y1 = 0, y2 = 0;

        x1 = e.reg[0].coord.getX();
        x2 = e.reg[1].coord.getX();
        y1 = e.reg[0].coord.getY();
        y2 = e.reg[1].coord.getY();

        // if the distance between the two points this line 
        // was created from is
        // less than the square root of 2, then ignore it
        if (Math.sqrt(((x2 - x1) * (x2 - x1)) 
            + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites) {
            return;
        }
        pxmin = borderMinX;
        pxmax = borderMaxX;
        pymin = borderMinY;
        pymax = borderMaxY;

        if (e.a == 1.0 && e.b >= 0.0) {
            s1 = e.ep[1];
            s2 = e.ep[0];
        } else {
            s1 = e.ep[0];
            s2 = e.ep[1];
        }

        if (e.a == 1.0) {
            y1 = pymin;
            if (s1 != null && s1.coord.getY() > pymin) {
                y1 = s1.coord.getY();
            }
            if (y1 > pymax) {
                y1 = pymax;
            }
            x1 = e.c - e.b * y1;
            y2 = pymax;
            if (s2 != null && s2.coord.getY() < pymax) {
                y2 = s2.coord.getY();
            }

            if (y2 < pymin) {
                y2 = pymin;
            }
            x2 = (e.c) - (e.b) * y2;
            if (((x1 > pxmax) & (x2 > pxmax)) | ((x1 < pxmin) & (x2 < pxmin))) {
                return;
            }
            if (x1 > pxmax) {
                x1 = pxmax;
                y1 = (e.c - x1) / e.b;
            }
            if (x1 < pxmin) {
                x1 = pxmin;
                y1 = (e.c - x1) / e.b;
            }
            if (x2 > pxmax) {
                x2 = pxmax;
                y2 = (e.c - x2) / e.b;
            }
            if (x2 < pxmin) {
                x2 = pxmin;
                y2 = (e.c - x2) / e.b;
            }
        } else {
            x1 = pxmin;
            if (s1 != null && s1.coord.getX() > pxmin) {
                x1 = s1.coord.getX();
            }
            if (x1 > pxmax) {
                x1 = pxmax;
            }
            y1 = e.c - e.a * x1;
            x2 = pxmax;
            if (s2 != null && s2.coord.getX() < pxmax) {
                x2 = s2.coord.getX();
            }
            if (x2 < pxmin) {
                x2 = pxmin;
            }
            y2 = e.c - e.a * x2;
            if (((y1 > pymax) & (y2 > pymax)) | ((y1 < pymin) & (y2 < pymin))) {
                return;
            }
            if (y1 > pymax) {
                y1 = pymax;
                x1 = (e.c - y1) / e.a;
            }
            if (y1 < pymin) {
                y1 = pymin;
                x1 = (e.c - y1) / e.a;
            }
            if (y2 > pymax) {
                y2 = pymax;
                x2 = (e.c - y2) / e.a;
            }
            if (y2 < pymin) {
                y2 = pymin;
                x2 = (e.c - y2) / e.a;
            }
        }
        
        if (Math.sqrt(((x2 - x1) * (x2 - x1)) 
            + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites) {
            return;
        }
       
        pushGraphEdge(e.reg[0], e.reg[1], x1, y1, x2, y2);
    }

    private void endPoint(Edge e, int lr, Site s) {
        e.ep[lr] = s;
        if (e.ep[RE - lr] == null) {
            return;
        }
        clipLine(e);
    }

    // returns 1 if p is to right of halfedge e
    private boolean rightOf(HalfEdge el, PairFloat p) {
        Edge e;
        Site topsite;
        boolean rightOfSite;
        boolean above, fast;
        float dxp, dyp, dxs, t1, t2, t3, yl;

        e = el.ELedge;
        topsite = e.reg[1];
        if (p.getX() > topsite.coord.getX()) {
            rightOfSite = true;
        } else {
            rightOfSite = false;
        }
        if (rightOfSite && el.ELpm == LE) {
            return true;
        }
        if (!rightOfSite && el.ELpm == RE) {
            return false;
        }

        if (e.a == 1.0) {
            dyp = p.getY() - topsite.coord.getY();
            dxp = p.getX() - topsite.coord.getX();
            fast = false;
            if ((!rightOfSite & (e.b < 0.0)) | (rightOfSite & (e.b >= 0.0))) {
                above = dyp >= e.b * dxp;
                fast = above;
            } else {
                above = p.getX() + p.getY() * e.b > e.c;
                if (e.b < 0.0) {
                    above = !above;
                }
                if (!above) {
                    fast = true;
                }
            }
            if (!fast) {
                dxs = topsite.coord.getX() - (e.reg[0]).coord.getX();
                above = e.b * (dxp * dxp - dyp * dyp) < dxs * dyp
                    * (1.0 + 2.0 * dxp / dxs + e.b * e.b);
                if (e.b < 0.0) {
                    above = !above;
                }
            }
        } else /* e.b==1.0 */ {
            yl = e.c - e.a * p.getX();
            t1 = p.getY() - yl;
            t2 = p.getX() - topsite.coord.getX();
            t3 = yl - topsite.coord.getY();
            above = t1 * t1 > t2 * t2 + t3 * t3;
        }
        return el.ELpm == LE ? above : !above;
    }

    private Site rightReg(HalfEdge he) {
        if (he.ELedge == (Edge) null) {
            // if this halfedge has no edge, return the bottom site (whatever
            // that is)
        
            return bottomsite;
        }

        // if the ELpm field is zero, return the site 0 that this edge bisects,
        // otherwise return site number 1
        return he.ELpm == LE ? he.ELedge.reg[RE] : 
            he.ELedge.reg[LE];
    }

    private float dist(Site s, Site t) {
        float dx, dy;
        dx = s.coord.getX() - t.coord.getX();
        dy = s.coord.getY() - t.coord.getY();
        return (float) (Math.sqrt(dx * dx + dy * dy));
    }

    // create a new site where the HalfEdges el1 and el2 intersect - note that
    // the Point in the argument list is not used, don't know why it's there
    private Site intersect(HalfEdge el1, HalfEdge el2) {
        Edge e1, e2, e;
        HalfEdge el;
        float d, xint, yint;
        boolean rightOfSite;
        Site v;

        e1 = el1.ELedge;
        e2 = el2.ELedge;
        if (e1 == null || e2 == null) {
            return null;
        }

        // if the two edges bisect the same parent, return null
        if (e1.reg[1] == e2.reg[1]) {
            return null;
        }

        d = e1.a * e2.b - e1.b * e2.a;
        if (-1.0e-10 < d && d < 1.0e-10) {
            return null;
        }

        xint = (e1.c * e2.b - e2.c * e1.b) / d;
        yint = (e2.c * e1.a - e1.c * e2.a) / d;

        if ((e1.reg[1].coord.getY() < e2.reg[1].coord.getY())
            || (e1.reg[1].coord.getY() == e2.reg[1].coord.getY() 
            && e1.reg[1].coord.getX() < e2.reg[1].coord.getX())) {
            el = el1;
            e = e1;
        } else {
            el = el2;
            e = e2;
        }

        rightOfSite = xint >= e.reg[1].coord.getX();
        if ((rightOfSite && el.ELpm == LE)
            || (!rightOfSite && el.ELpm == RE)) {
            return null;
        }

        // create a new site at the point of intersection - this is a new vector
        // event waiting to happen
        v = new Site();
        v.coord = new PairFloat(xint, yint);
        return v;
    }
    
}
