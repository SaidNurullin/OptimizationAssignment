import java.util.List;

public class Ezekiel {
    /**
     * Example:<br>
     * constraintsMatrix =<br>
     * 1 2 3 0 4 0<br>
     * 0 6 4 0 1 1<br>
     * 0 3 2 1 6 0<br>
     * basis =<br>
     * 1 0 0<br>
     * 0 1 0<br>
     * 0 0 1<br>
     * indices: {0, 5, 3}
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @return Array List of indices of columns that construct the basis in the constraints matrix as in the example above.
     */
    public static List<Integer> basisVectorsIndices(List<List<Double>> constraintsMatrix){
        return null;
    }

    /**
     *
     * Algorithm:<br>
     * For each column of constraintsMatrix, delta of this column is
     * the sum of results of following expression applied to each element of this column:<br>
     * (element of basisCoefficients with index = index of current row of constraintsMatrix) *
     * (current element of constraintsMatrix) - (element of coefficients with index = index of current column of constraintsMatrix)
     *
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @param coefficients Array List of doubles,
     *                    each element is the coefficient before the variable that is inside the objective function
     * @param basisCoefficients Array List of doubles,
     *                          each element is the coefficient before the variable from objective function
     *                          that construct the basis in the constraints matrix
     * @return Array List of deltas calculated using algorithm above.
     */
    public static List<Double> calculateDeltas(List<List<Double>> constraintsMatrix,
                                               List<Double> coefficients, List<Double> basisCoefficients){
        return null;
    }
}
