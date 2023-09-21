import java.util.List;

public class Said {
    /**
     *
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @param coefficients Array List of doubles,
     *                    each element is the coefficient before the variable that is inside the objective function
     * @return Array List of basis coefficients.
     */
    public static List<Double> basisCoefficients(List<List<Double>> constraintsMatrix, List<Double> coefficients){
        return null;
    }

    /**
     *
     * @param maxDeltaColumnIndex Index of the column inside the constraintsMatrix,
     *                            where the calculated delta is the biggest positive
     * @param rightNumbers ArrayList of doubles, each element is number from the right side of constraint.
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @return Array List of ratios (rightNumbers / column of constraintsMatrix with index = maxDeltaColumnIndex)
     */
    public static List<Double> calculateRatios(int maxDeltaColumnIndex, List<Double> rightNumbers,
                                               List<List<Double>> constraintsMatrix){
        return null;
    }

    /**
     *
     * @param basisVectorIndices Array List of indices of columns that construct the basis in the constraints matrix
     * @param maxDeltaColumnIndex Index of the column inside the constraintsMatrix,
     *                            where the calculated delta is the biggest positive
     * @param minRatioRowIndex Index of the row inside the constraintsMatrix,
     *                            which has the minimum positive ratio
     * @return new basisVectorIndices after replacement
     */
    public static List<Integer> newBasisVectorsIndices(List<Integer> basisVectorIndices,
                                                       int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    /**
     *
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @param maxDeltaColumnIndex Index of the column inside the constraintsMatrix,
     *                            where the calculated delta is the biggest positive
     * @param minRatioRowIndex Index of the row inside the constraintsMatrix,
     *                            which has the minimum positive ratio
     * @return new constraintsMatrix after applying the algorithm of triangles
     */
    public static List<List<Double>> newConstraintsMatrix(List<List<Double>> constraintsMatrix,
                                                          int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }
}
