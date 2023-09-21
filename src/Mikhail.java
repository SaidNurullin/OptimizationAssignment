import java.util.List;

public class Mikhail {
    /**
     * Algorithm:<br>
     * Replace the element of basisCoefficients with index = minRatioRowIndex with element of coefficients
     * with index = maxDeltaColumnIndex
     * @param basisCoefficients Array List of doubles,
     *                          each element is the coefficient before the variable from objective function
     *                          that construct the basis in the constraints matrix
     * @param coefficients Array List of doubles,
     *                    each element is the coefficient before the variable that is inside the objective function
     * @param maxDeltaColumnIndex Index of the column inside the constraintsMatrix,
     *                            where the calculated delta is the biggest positive
     * @param minRatioRowIndex Index of the row inside the constraintsMatrix,
     *                            which has the minimum positive ratio
     * @return new basisCoefficients using algorithm above.
     */
    public static List<Double> newBasisCoefficients(List<Double> basisCoefficients, List<Double> coefficients,
                                                    int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    /**
     * Algorithm:<br>
     * sum all the results of following expression applied to each element of rightNumbers:<br>
     * (value of element of basisCoefficients with index = index of current element of rightNumbers) *
     * (value current element of rightNumbers)
     *
     * @param rightNumbers ArrayList of doubles, each element is number from the right side of constraint.
     * @param basisCoefficients Array List of doubles,
     *                          each element is the coefficient before the variable from objective function
     *                          that construct the basis in the constraints matrix
     * @return the answer using algorithm above.
     */
    public static double calculateAnswer(List<Double> rightNumbers, List<Double> basisCoefficients){
        return 0;
    }
}
