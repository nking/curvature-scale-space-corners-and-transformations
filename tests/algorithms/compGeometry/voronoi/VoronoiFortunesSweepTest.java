package algorithms.compGeometry.voronoi;

import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.GraphEdge;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class VoronoiFortunesSweepTest extends TestCase {
    
    public VoronoiFortunesSweepTest() {
    }
    
    public void est0() throws FileNotFoundException, IOException {
        
        int nCities = 129;
        float[] x = new float[nCities];
        float[] y = new float[nCities];
        int[] population = new int[nCities];
        String[] names = new String[nCities];
        
        float xmin = Float.MAX_VALUE;
        float xmax = Float.MIN_VALUE;
        float ymin = Float.MAX_VALUE;
        float ymax = Float.MIN_VALUE;

        String path = ResourceFinder.findFileInTestResources("miles.dat");

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(path)));
            
            String line = null;

            String pattern = "^([a-zA-Z\\s-]+,\\s[a-zA-Z]{2})\\[(\\d+),(\\d+)\\](\\d+)$";
            Pattern p = Pattern.compile(pattern);
            Matcher m = null;
            line = in.readLine();
            do {
                //Youngstown, OH[4110,8065]115436
                m = p.matcher(line);
                if (m.matches()) {
                    String name = m.group(1);
                    names[count] = name;
                    x[count] =  ( 13000 - Integer.valueOf(m.group(3)) );
                    y[count] = ( Integer.valueOf(m.group(2)) );
                    population[count] = ( Integer.valueOf(m.group(4)) );

                    if (x[count] < xmin) {
                        xmin = x[count];
                    }
                    if (x[count] > xmax) {
                        xmax = x[count];
                    }
                    if (y[count] < ymin) {
                        ymin = y[count];
                    }
                    if (y[count] > ymax) {
                        ymax = y[count];
                    }
                    count++;
                }
                line = in.readLine();

                if (count > (nCities - 1) ) break;

            } while (line != null);
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        x = Arrays.copyOf(x, count);
        y = Arrays.copyOf(y, count);
        population = Arrays.copyOf(population, count);
        names = Arrays.copyOf(names, count);
        nCities = count;

        xmin = (xmin > 0) ? xmin : xmin;
        ymin = (ymin > 0) ? ymin : ymin;
        float min = (xmin < ymin) ? xmin : ymin;
        float max = (xmax > ymax) ? xmax : ymax;

        min = 0.7f*min;
        max = 1.05f*max;
        
        //TODO: revise this
        int minDist = 10;

        VoronoiFortunesSweep voronoi = 
            new VoronoiFortunesSweep();
        
        voronoi.generateVoronoi(x, y, min, max, min, max, 
            minDist);
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        assertNotNull(edges);
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(min, max, min, max);
        
        float[] xPolygon = null;
        float[] yPolygon = null;
        
        plotter.addPlot(x, y, xPolygon, yPolygon, "points");
        
        int n = edges.size();
        xPolygon = new float[2*n];
        yPolygon = new float[2*n];
        count = 0;
        for (GraphEdge edge : edges) {
            xPolygon[count] = edge.x1;
            yPolygon[count] = edge.y1;
            xPolygon[count + 1] = edge.x2;
            yPolygon[count + 1] = edge.y2;
            count += 2;
        }
        plotter.addPlot(xPolygon, yPolygon, xPolygon, yPolygon, 
            "edhes");
        String filePath = plotter.writeFile();
        System.out.println("wrote file=" + filePath);
    }
    
    public void test1() throws FileNotFoundException, IOException {
        
        float[] x = new float[]{ 43.f, 85.f,  85.f,  85.f, 127.f, 169.f, 211.f, 253.f, 253.f, 253.f};
        float[] y = new float[]{121.f, 79.f, 163.f, 205.f, 121.f, 205.f, 163.f,  79.f, 121.f, 205.f};

        float xmin = 0.f;
        float xmax = 300.f;
        float ymin = 0.f;
        float ymax = 300.f;

        //TODO: revise this
        int minDist = 10;

        VoronoiFortunesSweep voronoi = 
            new VoronoiFortunesSweep();
        
        voronoi.generateVoronoi(x, y, xmin, xmax, ymin, ymax, 
            minDist);
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        assertNotNull(edges);
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
        
        float[] xPolygon = null;
        float[] yPolygon = null;
        
        plotter.addPlot(x, y, xPolygon, yPolygon, "points");
        
        int n = edges.size();
        xPolygon = new float[2*n];
        yPolygon = new float[2*n];
        int count = 0;
        for (GraphEdge edge : edges) {
            xPolygon[count] = edge.x1;
            yPolygon[count] = edge.y1;
            xPolygon[count + 1] = edge.x2;
            yPolygon[count + 1] = edge.y2;
            count += 2;
        }
        plotter.addPlot(xPolygon, yPolygon, 
            xPolygon, yPolygon, "edhes");
        String filePath = plotter.writeFile();
        System.out.println("wrote file=" + filePath);
    }
}
