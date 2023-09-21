import java.util.List;

public class Daniil {
    /**
     *
     * @param ratios Array list of elements of type double (rightNumbers / maxDeltaColumn)
     * @return the index of the minimum positive element.
     */
    public static int minRatioRowIndex(List<Double> ratios){
        return 0;
    }

    /**
     *
     * Algorithm:<br>
     * 1. Divide the element of rightNumbers with index that is equal to minRatioRowIndex (let call it "special element")
     * by the element of constraintsMatrix located in the intersection of column with index = maxDeltaColumnIndex
     * and row with index = minRatioRowIndex.<br>
     * 2. For each rest element of rightNumbers:<br>
     * new element = (old element) - (element of the constraintsMatrix with column index = maxDeltaColumnIndex and
     * row index = index of old element == index of new element) * ("special element")
     *
     * @param rightNumbers ArrayList of doubles, each element is number from the right side of constraint.
     * @param constraintsMatrix 2 dimensional matrix of doubles,
     *                          each element is a coefficient before the variable that is inside the constraint
     * @param maxDeltaColumnIndex Index of the column inside the constraintsMatrix,
     *                            where the calculated delta is the biggest positive
     * @param minRatioRowIndex Index of the row inside the constraintsMatrix,
     *                            which has the minimum positive ratio
     * @return New rightNumbers using the algorithm above.
     */
    public static List<Double> newRightNumbers(List<Double> rightNumbers, List<List<Double>> constraintsMatrix,
                                               int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }
}
