package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * iterate over every pair to find those with distances less than threshhold*threshholdFactor
 * from one another.
 * 
 * runtime complexity O(N^2)
 * 
 * @author nichole
 */
public class BruteForceGroupFinder extends AbstractGroupFinder {

    public BruteForceGroupFinder(float threshhold, float threshholdFactor) {
        super(threshhold, threshholdFactor);
    }
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }
    
    @Override
    protected void findClusters(DoubleAxisIndexer indexer) {
        
        SimpleLinkedListNode[] tmpGroupMembership = new SimpleLinkedListNode[10];
        int nTmpGroups = 0;
        for (int i = 0; i < tmpGroupMembership.length; i++) {
            tmpGroupMembership[i] = new SimpleLinkedListNode();
        }

        int[] tmpPointToGroupIndex = new int[indexer.getNXY()];
        Arrays.fill(tmpPointToGroupIndex, -1);

        // the indexes stored in the instance vars such as pointToGroupIndex are w.r.t. the arrays sorted by y
        float[] x = indexer.getXSortedByY();
        float[] y = indexer.getYSortedByY();

        float densityThreshold = threshholdFactor*threshhold;
        
        // 2 or 3 times the background density
        if (debug) {
            log.info("For clusters, using densityThreshold=" + densityThreshold);
        }

        SimpleLinkedListNode tmpPointsInMoreThanOneGroup = new SimpleLinkedListNode();

        // store the distances between points when the separation implies it is above density
        for (int i = 0; i < indexer.getNXY(); i++) {

            // groups are points whose separation is an implied density higher than 2 or 3 times the background.
            // see if this point is part of a group already, by looking for it in the pointTo2ndPointIndex entry.
            // if it is a part of a group already, set the groupNumber,
            // else, leave the groupNumber as -1 until another point is close enough to start a new group.

            // this will be -1 if not yet assigned
            int groupNumber = tmpPointToGroupIndex[i];

            for (int j = (i + 1); j < indexer.getNXY(); j++) {

                float d = (float) Math.sqrt( Math.pow((x[i] - x[j]), 2) + Math.pow((y[i] - y[j]), 2) );

                float dens = 2.f/d;

                if (dens > densityThreshold) {

                    if (tmpGroupMembership.length < (nTmpGroups + 1)) {
                        int oldN = tmpGroupMembership.length;
                        int n = (int) (1.5f * oldN);
                        if (n < (oldN + 1)) {
                            n = oldN + 1;
                        }
                        tmpGroupMembership = Arrays.copyOf(tmpGroupMembership, n);
                        for (int k = oldN; k < n; k++) {
                            tmpGroupMembership[k] = new SimpleLinkedListNode();
                        }
                    }

                    if (groupNumber == -1) {

                        if (tmpPointToGroupIndex[j] != -1) {

                            groupNumber = tmpPointToGroupIndex[j];

                        } else {

                            groupNumber = nTmpGroups;
                        }
                        tmpGroupMembership[groupNumber].insert(i);
                        tmpPointToGroupIndex[i] = groupNumber;
                        nTmpGroups++;

                    }// else, is already a member of a group

                    if ( (tmpPointToGroupIndex[j] != -1) && ( tmpPointToGroupIndex[j] != groupNumber ) ) {

                        tmpPointsInMoreThanOneGroup.insertIfDoesNotAlreadyExist(j);
                    }

                    tmpGroupMembership[groupNumber].insertIfDoesNotAlreadyExist(j);

                    tmpPointToGroupIndex[j] = groupNumber;
                }
            }
        }

        // consolidate overlapping groups.  high density backgrounds will have larger number of overlapping clusters
        SimpleLinkedListNode latest = tmpPointsInMoreThanOneGroup;
        while (latest != null && latest.key != -1) {

            int pointIndex = latest.key;

            int gn1 = tmpPointToGroupIndex[pointIndex];

            for (int i = 0; i < nTmpGroups; i++) {

                if (i == gn1) {
                    continue;
                }

                SimpleLinkedListNode grNode = tmpGroupMembership[i].search(pointIndex);

                if (grNode != null) {

                    int n1 = tmpGroupMembership[gn1].getNumberOfKeys();
                    int n2 = tmpGroupMembership[i].getNumberOfKeys();

                    if (n1 > n2) {
                        //put all of group 2 points into group 1
                        SimpleLinkedListNode group2 = tmpGroupMembership[i];
                        while ( (group2 != null) && (group2.key != -1) ) {

                            int aPointInG2 = group2.key;

                            tmpPointToGroupIndex[aPointInG2] = gn1;

                            tmpGroupMembership[gn1].insertIfDoesNotAlreadyExist(aPointInG2);

                            tmpGroupMembership[i].delete(aPointInG2);

                            //group2 = group2.next;
                            //instead of using next, start from top again due to delete while in iteration
                            group2 = tmpGroupMembership[i];
                        }
                    } else {
                        //put all of group 1 points into group 2
                        SimpleLinkedListNode group1 = tmpGroupMembership[gn1];
                        while ( (group1 != null) && (group1.key != -1) ) {

                            int aPointInG1 = group1.key;

                            tmpPointToGroupIndex[aPointInG1] = i;

                            tmpGroupMembership[i].insertIfDoesNotAlreadyExist(aPointInG1);

                            tmpGroupMembership[gn1].delete(aPointInG1);

                            //group2 = group2.next;
                            //instead of using next, start from top again due to delete while in iteration
                            group1 = tmpGroupMembership[gn1];
                        }
                    }

                    i = nTmpGroups;
                }
            }
            tmpPointsInMoreThanOneGroup.delete(pointIndex);
            latest = tmpPointsInMoreThanOneGroup;
        }

        // iterate over group data and store in instance vars when group membership is > 2
        for (int i = 0; i < nTmpGroups; i++) {

            int currentTmpGroupNumber = i;

            SimpleLinkedListNode grpMembersNode = tmpGroupMembership[currentTmpGroupNumber];

            int nGrpMembers = grpMembersNode.getNumberOfKeys();

            if (nGrpMembers >= minimumNumberInCluster) {

                // store these associated with current group number which is nGroups

                int groupNumber = nGroups;
                nGroups++;

                while ( (grpMembersNode != null) && (grpMembersNode.key != -1) ) {

                    int pointIndex = grpMembersNode.key;

                    checkAndExpandGroupMembershipArray();

                    pointToGroupIndex[pointIndex] = groupNumber;

                    groupMembership[groupNumber].insertIfDoesNotAlreadyExist(pointIndex);

                    grpMembersNode = grpMembersNode.next;
                }
            }
        }

    }
}
