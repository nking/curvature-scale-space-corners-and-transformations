package thirdparty;

import algorithms.util.OneDIntArray;
import java.util.*;

/**
 The algorithm below was adapted from code by Gary Baker

 * Edits were made to the code below to allow rectangular matrix with more
 * columns than rows.  Comments were also added.
 * 
 * An implementation of the classic hungarian algorithm for the assignment
 * problem.

   Hungarian method can solve the assignment problem in
    O(mn + n^2(log n)), where n := |X| = |Y | and m := |E|
       (Ramshaw & Tarjan, 2012)
 *
 * Copyright 2007 Gary Baker (GPL v3)
 *
 * @author gbaker
 */
public class HungarianAlgorithm {

    /**
     * NOTE: matrix is modified as a side effect.
     *
     * @param matrix cost matrix w/ first dimension being rows and 2nd being
     * columns.
     * @return minimum cost bipartite matching
     */
    public int[][] computeAssignments(float[][] matrix) {

        int nRows = matrix.length;
        int nCols = matrix[0].length;
        // there can be more columns than rows, but not vice versa
        if (nRows > nCols) {
            throw new IllegalArgumentException(
                "matrix must have same number or fewer rows than columns");
        }

        // subtract minumum value from rows and columns to create lots of zeroes
        reduceMatrix(matrix);

        // non negative values are the index of the starred or primed zero in the row or column
        int[] starsByRow = new int[matrix.length];
        Arrays.fill(starsByRow, -1);
        int[] starsByCol = new int[matrix[0].length];
        Arrays.fill(starsByCol, -1);
        int[] primesByRow = new int[matrix.length];
        Arrays.fill(primesByRow, -1);

        // 1s mean covered, 0s mean not covered
        int[] coveredRows = new int[matrix.length];
        int[] coveredCols = new int[matrix[0].length];

        // star any zero that has no other starred zero in the same row or column
        // populate starsByRow: each row has the assigned column
        // populate starsByCol: each col has the assigned row
        initStars(matrix, starsByRow, starsByCol);
        // populate coveredCols: has a zero where unassigned
        coverColumnsOfStarredZeroes(starsByCol, coveredCols);

        // while still an unassigned column
        while (!allAreCovered(coveredCols)) {

            int[] primedZero = primeSomeUncoveredZero(matrix, primesByRow,
                coveredRows, coveredCols);

            boolean changed = true;

            while (primedZero == null) {
                // keep making more zeroes until we find something that we can 
                // prime (i.e. a zero that is uncovered)

                changed = makeMoreZeroes(matrix, coveredRows, coveredCols);

                if (!changed) {
                    break;
                }

                primedZero = primeSomeUncoveredZero(matrix, primesByRow,
                    coveredRows, coveredCols);
            }

            if (!changed) {
                break;
            }

            // check if there is a starred zero in the primed zero's row
            //(columnIndex is > -1 if it's assignable)
            int columnIndex = starsByRow[primedZero[0]];
            if (-1 == columnIndex) {

                // if not, then we need to increment the zeroes and start over
                incrementSetOfStarredZeroes(primedZero, starsByRow, starsByCol, 
                    primesByRow);
                Arrays.fill(primesByRow, -1);
                Arrays.fill(coveredRows, 0);
                Arrays.fill(coveredCols, 0);
                coverColumnsOfStarredZeroes(starsByCol, coveredCols);
            } else {

                // cover the row of the primed zero and uncover the column of 
                // the starred zero in the same row
                coveredRows[primedZero[0]] = 1;
                coveredCols[columnIndex] = 0;
            }
        }

        // ok now we should have assigned everything
        int[][] retval = new int[matrix.length][];
        
        /*
        // take the starred zeroes in each column as the correct assignments
        for (int i = 0; i < starsByCol.length; i++) {
            retval[i] = new int[]{starsByCol[i], i};
        }
        */
        // choosing by row instead to allow uneven set sizes.
        //   there can be more columns than rows
        for (int i = 0; i < starsByRow.length; i++) {
            retval[i] = new int[]{i, starsByRow[i]};
        }
        
        return retval;
    }

    private boolean allAreCovered(int[] coveredCols) {
        for (int covered : coveredCols) {
            if (0 == covered) {
                return false;
            }
        }
        return true;
    }

    /**
     * the first step of the hungarian algorithm is to find the smallest element
     * in each row and subtract it's values from all elements in that row
     *
     * @return the next step to perform
     */
    private void reduceMatrix(float[][] matrix) {

        for (int i = 0; i < matrix.length; i++) {

            // find the min value in the row
            float minValInRow = Float.MAX_VALUE;
            for (int j = 0; j < matrix[i].length; j++) {
                if (minValInRow > matrix[i][j]) {
                    minValInRow = matrix[i][j];
                }
            }

            // subtract it from all values in the row
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] -= minValInRow;
            }
        }

        for (int i = 0; i < matrix[0].length; i++) {
            float minValInCol = Float.MAX_VALUE;
            for (int j = 0; j < matrix.length; j++) {
                if (minValInCol > matrix[j][i]) {
                    minValInCol = matrix[j][i];
                }
            }

            for (int j = 0; j < matrix.length; j++) {
                matrix[j][i] -= minValInCol;
            }

        }

    }

    /**
     * init starred zeroes
     *
     * for each column find the first zero if there is no other starred zero in
     * that row then star the zero, cover the column and row and go onto the
     * next column
     *
     * @param costMatrix
     * @param starredZeroes
     * @param coveredRows
     * @param coveredCols
     * @return the next step to perform
     */
    private void initStars(float costMatrix[][], int[] starsByRow, int[] starsByCol) {

        int[] rowHasStarredZero = new int[costMatrix.length];
        int[] colHasStarredZero = new int[costMatrix[0].length];

        for (int i = 0; i < costMatrix.length; i++) {
            for (int j = 0; j < costMatrix[i].length; j++) {
                // if cell cost is now a zero (=min) and is not assigned
                if (0 == costMatrix[i][j] && 0 == rowHasStarredZero[i] && 
                    0 == colHasStarredZero[j]) {

                    //row i can be assigned column j:
                    starsByRow[i] = j;

                    // col j can be assigned row i:
                    starsByCol[j] = i;

                    // mark row i as assigned:
                    rowHasStarredZero[i] = 1;

                    // mark col j as assigned:
                    colHasStarredZero[j] = 1;

                    break; // move onto the next row
                }
            }
        }
    }

    /**
     * just marke the columns covered for any column containing a starred zero.
     * coveredCols is populated as 1's where starsByCol is already assigned,
     * else coveredCols has a 0.
     *
     * @param starsByCol
     * @param coveredCols
     */
    private void coverColumnsOfStarredZeroes(int[] starsByCol, int[] coveredCols) {
        for (int j = 0; j < starsByCol.length; j++) {
            coveredCols[j] = (-1 == starsByCol[j]) ? 0 : 1;
        }
    }

    /**
     * finds some uncovered zero and primes it. populates primesByRow with the
     * first possible assignable column for a row and returns int[]{assignable
     * row index, assignable column index}.
     *
     * @param matrix
     * @param primesByRow
     * @param coveredRows
     * @param coveredCols
     * @return
     */
    private int[] primeSomeUncoveredZero(float matrix[][], int[] primesByRow,
        int[] coveredRows, int[] coveredCols) {

        // find an uncovered zero and prime it
        for (int i = 0; i < matrix.length; i++) {
            if (1 == coveredRows[i]) {
                continue;
            }
            for (int j = 0; j < matrix[i].length; j++) {
                // if it's a zero and the column is not covered
                if (0 == matrix[i][j] && 0 == coveredCols[j]) {

                    // ok this is an unstarred (i.e. unassigned) zero
                    // prime it
                    primesByRow[i] = j;
                    return new int[]{i, j};
                }
            }
        }

        //didn't find an assignable row with an assignable column
        return null;

    }

    /**
     *
     * @param unpairedZeroPrime the unassigned {row index, col index}
     * @param starsByRow array holding columns that are assignable
     * @param starsByCol array holding rows that are assignable
     * @param primesByRow array holding the assignable columns
     */
    private void incrementSetOfStarredZeroes(int[] unpairedZeroPrime,
        int[] starsByRow, int[] starsByCol, int[] primesByRow) {

        int j, i = unpairedZeroPrime[1];

        Set<OneDIntArray> zeroSequence = new LinkedHashSet<OneDIntArray>();
        zeroSequence.add(new OneDIntArray(unpairedZeroPrime));
        boolean paired = false;
        do {
            j = starsByCol[i];
            paired = (-1 != j) && zeroSequence.add(new OneDIntArray(new int[]{j, i}));
            if (!paired) {
                break;
            }

            i = primesByRow[j];
            paired = -1 != i && zeroSequence.add(new OneDIntArray(new int[]{j, i}));

        } while (paired);

        // unstar each starred zero of the sequence
        // and star each primed zero of the sequence
        for (OneDIntArray zeroW : zeroSequence) {
            int[] zero = zeroW.a;
            if (starsByCol[zero[1]] == zero[0]) {
                starsByCol[zero[1]] = -1;
                starsByRow[zero[0]] = -1;
            }
            if (primesByRow[zero[0]] == zero[1]) {
                starsByRow[zero[0]] = zero[1];
                starsByCol[zero[1]] = zero[0];
            }
        }
    }

    /**
     * return true if successfully unset a previously assigned row, else returns
     * false.
     *
     * @param matrix
     * @param coveredRows
     * @param coveredCols
     * @return
     */
    private boolean makeMoreZeroes(float[][] matrix, int[] coveredRows,
        int[] coveredCols) {

        // find the minimum uncovered value
        float minUncoveredValue = Float.MAX_VALUE;
        for (int i = 0; i < matrix.length; i++) {
            // row i is assignable
            if (0 == coveredRows[i]) {
                for (int j = 0; j < matrix[i].length; j++) {
                    // col j is assignable
                    if (0 == coveredCols[j] && matrix[i][j] < minUncoveredValue) {
                        minUncoveredValue = matrix[i][j];
                    }
                }
            }
        }

        if (minUncoveredValue == Float.MAX_VALUE) {
            return false;
        }

        boolean didNotChange = true;

        // add the min value to all covered rows
        for (int i = 0; i < coveredRows.length; i++) {
            if (1 == coveredRows[i]) {
                // coveredRows assigned row i is now reset to higher than 0 for all columns
                for (int j = 0; j < matrix[i].length; j++) {
                    matrix[i][j] += minUncoveredValue;
                    didNotChange = false;
                }
            }
        }

        // subtract the min value from all uncovered columns
        for (int j = 0; j < coveredCols.length; j++) {
            if (0 == coveredCols[j]) {
                for (int i = 0; i < matrix.length; i++) {
                    matrix[i][j] -= minUncoveredValue;
                    didNotChange = false;
                }
            }
        }

        return !didNotChange;
    }

}
