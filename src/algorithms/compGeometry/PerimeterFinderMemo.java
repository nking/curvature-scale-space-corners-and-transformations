package algorithms.compGeometry;

import algorithms.util.PairInt;
import algorithms.util.PathStep;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

/**
 * 
 * @author nichole
 */
public class PerimeterFinderMemo {
    
    private Map<PairInt, Set<PathStep>> noSolutionMap = new HashMap<PairInt, Set<PathStep>>();
   
    public void storeNoSolutionState(LinkedList<PathStep> path) {
        
        /*
        store entire path that lead to no solution:
        
        last node:      (x,y), set<pairint> possibleStepsAtInstantiation,  integer allVisitedBitvector, set<pairint> availMoves
        last node - 1:  (x,y), set<pairint> possibleStepsAtInstantiation,  integer allVisitedBitvector, set<pairint> availMoves
        ...
        0            :  (x,y), set<pairint> possibleStepsAtInstantiation,  integer allVisitedBitvector, set<pairint> availMoves
        
        store each key and node
        */
        
        for (PathStep step : path) {
            PairInt p = step.getCoords();
            Set<PathStep> nsNodes = noSolutionMap.get(p);
            if (nsNodes == null) {
                nsNodes = new HashSet<PathStep>();
                noSolutionMap.put(p, nsNodes);
            }
            nsNodes.add(step.copy());
        }
    }
    
    public boolean isANoSolutionState(PathStep step) {
                
        Set<PathStep> states = noSolutionMap.get(step.getCoords());
        
        if (states == null) {
            return false;
        }
        
        return states.contains(step);
    }
}
