import java.util.List;

public class Timur {
    /**
     *
     * @param deltas Array List of doubles, each element is delta calculated for each column from constraints matrix
     * @return the index of maximum positive value among deltas, if no positive deltas, returns -1.
     */
    public static int maxDeltaColumnIndex(List<Double> deltas){
        return -1;
    }

    /**
     * Put values from rightNumbers to elements of output Array List
     * with indices corresponding to the values of basisVectorsIndices Array List
     * Example:<br>
     * basisVectorsIndices = {4, 5, 3}<br>
     * rightNumbers = {11, 8, 36}<br>
     * size = 7<br>
     * output= {0, 0, 0, 36, 11, 8, 0, 0}
     *
     * @param basisVectorsIndices Array List of indices of columns that construct the basis in the constraints matrix
     * @param rightNumbers ArrayList of doubles, each element is number from the right side of constraint.
     * @param size Size of output Array List
     * @return Array List of values of variables from objective function
     */
    public static List<Double> decisionVariables(List<Integer> basisVectorsIndices, List<Double> rightNumbers, int size){
        return null;
    }
}
