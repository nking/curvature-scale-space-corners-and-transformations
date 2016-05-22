package algorithms.imageProcessing.optimization.segmentation;

/**
 Class to estimate the difference between 
 regions in images segmented by humans and by a segmentation
 algorithm.
 
 The class follows the benchmark map comparison section
 of "Learning to Detect Natural Image Boundaries
 Using Local Brightness, Color, and Texture Cues"
 by Martin, Fowlkes, and Malik.
 
 -- test data and human segmentation data (model) are 
    extracted to get the boundaries of the segmented regions.
 -- the weight between a data pixel and model pixel is
    the distance between them.
    -- all boundary pixels matched beyond some threshold 
       d_max are nonhits.
  -- the authors designed a method to determine bipartite
     min-cost pairings faster than O(n^2) using a
     sparse assignment problem. 
     They use Goldberg’s CSA package, 
     which implements the best known algorithms
     for min-cost sparse assignment [42], [43]. 
     The CSA code appears to run in time linear in the 
     size of the graph, O(n), so they designed a sparse graph
     solution adapted to use these structures:
     [42] A. Goldberg and R. Kennedy, 
     “An Efficient Cost Scaling Algorithm
     for the Assignment Problem,” SIAM J. Discrete Math., 1993.
     [43] B.V. Cherkassky and A. V. Goldberg, 
     “On Implementing PushRelabel Method for the 
     Maximum Flow Problem,” Proc. Fourth Integer 
     Programming and Combinatorial Optimization Conf., 
     pp. 157-171, May 1995.
     
     - they include in the graph only those edges 
       with weight w less than or equal to d_max, 
       since an edge with w > d_max marks a missed
       human boundary to a test false positive.
     - then they remove isolated nodes and count 
       them as a "miss" or "false positive"
     perfect matching is needed because the degree
     of the graph cannot be known apriori after 
     sparsification.
     
     - they add outlier nodes on both sides of the match.
       All edges incident on an outlier node have higher 
       weight than any real edge in the graph, ensuring 
       that they are used only when necessary to extend 
       the true min-cost partial matching to a valid
       perfect matching

       "Given a sparse assignment problem with n_L nodes 
       on the left side and n_R nodes on the right, we add 
       n_R outlier nodes to the left and n_L outlier nodes 
       to the right. This squared problem has enough nodes to 
       ensure a perfect matching."
       The outlier connections have identical weights and large
       redundancy.
       Keeping a constant number of outlier connections per 
       node keeps the size of the graph linear.
       The authors use d = 6 connectivity, resulting 
       in d random outlier connections to each real node 
       and d random outlier connections node.
       Another step to ensure perfect matching is 
       perfect matching of high cost connections to match
       each real node with an outlier mode in parallel.
       This false high cost step is before the random
       outlier connections are added without replacement.
       the result is best correspondence between the test
       boundary maps and human boundary maps with a
       maximum localization tolerance of d_max.
       
       the weakness in the method is regarding junctions,
       but because they are so few in number compared to
       points not in junctions, corrections are not needed.

       The code below ports the berkeley implementation in
       matlab to java here.
       http://www.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/
 
       Note that the berkeley segmentation research page
       includes data and code.
       The data have this restriction:
       
       "You are free to download a portion of the dataset for non-commercial research and educational purposes.  In exchange, we request only that you make available to us the results of running your segmentation or boundary detection algorithm on the test set as described below.  Work based on the dataset should cite our ICCV 2001 paper:
        @InProceedings{MartinFTM01,
          author = {D. Martin and C. Fowlkes and D. Tal and J. Malik},
          title = {A Database of Human Segmented Natural Images and its
                   Application to Evaluating Segmentation Algorithms and
                   Measuring Ecological Statistics},
          booktitle = {Proc. 8th Int'l Conf. Computer Vision},
          year = {2001},
          month = {July},
          volume = {2},
          pages = {416--423}
        }
       "
 
 For the test files, need to make soft boundary map, 
 preferably after non-max suppression, stored in a 
 BMP image file.
 Soft boundary map is a one pixel wide boundary, 
 valued from zero to one where high values signify 
 greater confidence in the existence of a boundary.
 
 ---> Either follow the remaining directions, or 
      use the benchNewAlg.m script.

(11) Run boundaryBench('BENCH/PRES/ALG','PRES').  
This will compute precision/recall and put the 
results in text files in the ALG directory.  
This takes 120 minutes on my 3GHz P-IV.

(12) Run boundaryBenchGraphs('BENCH/PRES/ALG').  
This will create precision/recall graphs in the 
ALG directory.  

--> If you have other algorithms to benchmark, 
    repeat steps 6-12 for each algorithm.

(13) Run boundaryBenchGraphsMulti('BENCH').  
This will create graphs comparing pairs of algorithms.

    fwrite(2,'  Calculating precision/recall (fast method) ');
    [thresh,cntR,sumR,cntP,sumP] = boundaryPRfast(pb,segs,nthresh);

    R = cntR ./ (sumR + (sumR==0));
    P = cntP ./ (sumP + (sumP==0));
    F = fmeasure(R,P);
    [bestT,bestR,bestP,bestF] = maxF(thresh,R,P);
  
    scores(i,:) = [iid bestT bestR bestP bestF];
    fname = fullfile(pbDir,'scores.txt');
    fid = fopen(fname,'w');
    if fid==-1, 
      error(sprintf('Could not open file %s for writing.',fname));
    end
    fprintf(fid,'%10d %10g %10g %10g %10g\n',scores(1:i,:)');
    fclose(fid);
  
    cntR_total = cntR_total + cntR;
    sumR_total = sumR_total + sumR;
    cntP_total = cntP_total + cntP;
    sumP_total = sumP_total + sumP;

    R = cntR_total ./ (sumR_total + (sumR_total==0));
    P = cntP_total ./ (sumP_total + (sumP_total==0));
    F = fmeasure(R,P);
    [bestT,bestR,bestP,bestF] = maxF(thresh,R,P);
  
    fname = fullfile(pbDir,'pr.txt');
    fid = fopen(fname,'w');
    if fid==-1, 
      error(sprintf('Could not open file %s for writing.',fname));
    end
  
    % compute f-measure fromm recall and precision
    function [f] = fmeasure(r,p)
    f = 2*p.*r./(p+r+((p+r)==0));

    % interpolate to find best F and coordinates thereof
    function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
    bestT = thresh(1);
    bestR = R(1);
    bestP = P(1);
    bestF = fmeasure(R(1),P(1));
    for i = 2:numel(thresh),
      for d = linspace(0,1),
        t = thresh(i)*d + thresh(i-1)*(1-d);
        r = R(i)*d + R(i-1)*(1-d);
        p = P(i)*d + P(i-1)*(1-d);
        f = fmeasure(r,p);
        if f > bestF,
          bestT = t;
          bestR = r;
          bestP = p;
          bestF = f;
        end
      end
    end

(14) Run boundaryBenchHtml('BENCH').  This will 
generate all the web pages for the benchmark.
 
 @author nichole
 */
public class BenchmarkMeasurer {
    
}
