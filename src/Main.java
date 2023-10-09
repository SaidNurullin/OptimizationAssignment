import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

public class Main {
    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @return A matrix containing a subset of columns from constraintsMatrix that
     *         construct identity matrix.
     */
    private static Matrix basisMatrix(Matrix constraintsMatrix) {
        int k = 0;
        int idx = 0;
        Matrix currentColumn = new Matrix(constraintsMatrix.getRows(), 1);
        Matrix basMatrix = new Matrix(constraintsMatrix.getRows(), constraintsMatrix.getRows());
        for (int j = 0; j < constraintsMatrix.getColumns(); j++) {
            int numb0 = 0, numb1 = 0;
            for (int i = 0; i < constraintsMatrix.getRows(); i++) {
                currentColumn.setElement(i, 0, constraintsMatrix.getElement(i, j));
                if (constraintsMatrix.getElement(i, j) == 0) {
                    numb0++;
                } else if (constraintsMatrix.getElement(i, j) == 1) {
                    numb1++;
                    idx = i;
                }
            }
            if (numb1 == 1 && (numb0 + numb1) == constraintsMatrix.getRows()) {
                for (int i = 0; i < currentColumn.getRows(); i++){
                    basMatrix.setElement(i, idx, currentColumn.getElement(i, 0));
                }
                k++;
                if(k == basMatrix.getRows()){
                    break;
                }
            }
        }

        return basMatrix;
    }

    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @return A matrix containing the columns from constraintsMatrix that not
     *         construct identity matrix
     */
    private static Matrix nonBasicMatrix(Matrix constraintsMatrix) {
        int k = 0;
        Matrix currentColumn = new Matrix(constraintsMatrix.getRows(), 1);
        Matrix basMatrix = new Matrix(constraintsMatrix.getRows(), constraintsMatrix.getColumns() - constraintsMatrix.getRows());
        for (int j = 0; j < constraintsMatrix.getColumns(); j++) {
            int numb0 = 0, numb1 = 0;
            for (int i = 0; i < constraintsMatrix.getRows(); i++) {
                currentColumn.setElement(i, 0, constraintsMatrix.getElement(i, j));
                if (constraintsMatrix.getElement(i, j) == 0) {
                    numb0++;
                } else if (constraintsMatrix.getElement(i, j) == 1) {
                    numb1++;
                }
            }
            if (!(numb1 == 1 && (numb0 + numb1) == constraintsMatrix.getRows())) {
                for (int i = 0; i < currentColumn.getRows(); i++){
                    basMatrix.setElement(i, k, currentColumn.getElement(i, 0));
                }
                k++;
            }
        }
        return basMatrix;
    }

    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @param coefficients      A matrix representing the coefficients of the
     *                          objective function.
     * @return Coefficients corresponding to the basis columns
     */
    private static Matrix basisCoefficients(Matrix constraintsMatrix, Matrix coefficients) {
        List<Integer> basicIndexes = basisVectorsIndices(constraintsMatrix);
        int rows = 1;
        Matrix matrix = new Matrix(rows, basicIndexes.size());
        int k = 0;
        for (Integer element : basicIndexes){
            matrix.setElement(rows - 1, k, coefficients.getElement(element, 0));
            ++k;
        }
        return matrix;
    }

    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @param coefficients      A matrix representing the coefficients of the
     *                          objective function.
     * @return Coefficients corresponding to the non-basic columns.
     */
    private static Matrix nonBasicCoefficients(Matrix constraintsMatrix, Matrix coefficients) {
        List<Integer> nonbasicIndexes = nonBasisVectorsIndices(constraintsMatrix);
        int rows = 1;
        Matrix matrix = new Matrix(rows, nonbasicIndexes.size());
        int k = 0;
        for (Integer element : nonbasicIndexes){
            matrix.setElement(rows-1, k, coefficients.getElement(element, 0));
            ++k;
        }
        return matrix;
    }

    /**
     *
     * @param nonBasicMatrix      A matrix containing the remaining columns from
     *                            constraintsMatrix (non-basic matrix P).
     * @param maxDeltaColumnIndex The index of the column with the maximum value in
     *                            deltas.
     * @return The column vector from the non-basic matrix corresponding to
     *         maxDeltaColumnIndex.
     */
    private static Matrix enteringVector(Matrix nonBasicMatrix, int maxDeltaColumnIndex) {
        Matrix enteringVector = new Matrix(nonBasicMatrix.getRows(), 1);
        for (int i = 0; i < enteringVector.getRows(); i++) {
            enteringVector.setElement(i, 0, nonBasicMatrix.getElement(i, maxDeltaColumnIndex));
        }

        return enteringVector;
    }

    /**
     *
     * @param basisVectorsValues    Values of the basic variables (XB).
     * @param denominatorsForRatios Values to which basisVectorsValues should be
     *                              divided
     * @return Vector of ratios (basisVectorsValues)/(denominatorsForRatios).
     */
    private static Matrix calculateRatios(Matrix basisVectorsValues, Matrix denominatorsForRatios) {
        Matrix ratios = new Matrix(basisVectorsValues.getRows(), 1);
        for (int i = 0; i < ratios.getRows(); i++) {
            ratios.setElement(i, 0, basisVectorsValues.getElement(i, 0) /
                    denominatorsForRatios.getElement(i, 0));
        }
        return ratios;
    }

    /**
     * Example:<br></br>
     * constraintsMatrix =<br></br>
     * 1 2 3 0 4 0<br></br>
     * 0 6 4 0 1 1<br></br>
     * 0 3 2 1 6 0<br></br>
     * basis =<br></br>
     * 1 0 0<br></br>
     * 0 1 0<br></br>
     * 0 0 1<br></br>
     * indices: {0, 5, 3}
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @return Array List of indices of columns that construct the basis in the
     *         constraints matrix as in the example above.
     */
    private static List<Integer> basisVectorsIndices(Matrix constraintsMatrix) {
        Map<Integer, Integer> mapIndices = new HashMap<Integer, Integer>();
        List<Integer> basisVIndices = new ArrayList<>();
        List<Integer> listRow = new ArrayList<>();
        int k = 0;
        for (int j = 0; j < constraintsMatrix.getColumns(); j++) {
            int numb0 = 0, numb1 = 0;
            for (int i = 0; i < constraintsMatrix.getRows(); i++) {
                if (constraintsMatrix.getElement(i, j) == 0) {
                    numb0++;
                } else if (constraintsMatrix.getElement(i, j) == 1) {
                    k = i;
                    numb1++;
                }
            }
            if ((numb1 == 1 && (numb0 + numb1) == constraintsMatrix.getRows())) {
                listRow.add(k);
                mapIndices.put(k, j);
                k++;
            }
        }
        listRow.sort(Comparator.naturalOrder());
        int size = listRow.size();
        for (int i = 0; i < size; i++) {
            basisVIndices.add(i, mapIndices.get(listRow.get(i)));
        }
        return basisVIndices;
    }

    /**
     * Example:<br></br>
     * constraintsMatrix =<br></br>
     * 1 2 3 0 4 0<br></br>
     * 0 6 4 0 1 1<br></br>
     * 0 3 2 1 6 0<br></br>
     * non-basic columns =<br></br>
     * 2 3 4<br></br>
     * 6 4 1<br></br>
     * 3 2 6<br></br>
     * indices: {1, 2, 4}
     *
     * @param constraintsMatrix A matrix representing the coefficients of the
     *                          constraints.
     * @return Array List of indices of columns that not construct the basis in the
     *         constraints matrix as in the example above.
     */
    private static List<Integer> nonBasisVectorsIndices(Matrix constraintsMatrix) {
        List<Integer> nonBasisVIndices = new ArrayList<Integer>();
        int k = 0;
        for (int j = 0; j < constraintsMatrix.getColumns(); j++) {
            int numb0 = 0, numb1 = 0;
            for (int i = 0; i < constraintsMatrix.getRows(); i++) {
                if (constraintsMatrix.getElement(i, j) == 0) {
                    numb0++;
                } else if (constraintsMatrix.getElement(i, j) == 1) {
                    numb1++;
                }
            }
            if (!(numb1 == 1 && (numb0 + numb1) == constraintsMatrix.getRows())) {
                nonBasisVIndices.add(k, j);
                k++;
            }
        }
        return nonBasisVIndices;
    }

    /**
     *
     * @param deltas Array List of doubles, each element is delta calculated for
     *               each column from constraints matrix
     * @return the index of maximum positive value among deltas, if no positive
     *         deltas, returns -1.
     */
    private static int maxDeltaColumnIndex(Matrix deltas) {
        double maxDelta = -1;
        int res = -1;
        for (int i = 0; i < deltas.getRows(); i++)
        {
            double deltaI = deltas.getElement(i, 0);
            if (maxDelta < deltaI && deltaI > 0) {
                maxDelta = deltaI;
                res = i;
            }
        }
        return res;
    }

    /**
     *
     * @param ratios Vector of ratios
     * @return The index of the row with the minimum ratio in ratios.
     */
    private static int minRatioRowIndex(Matrix ratios) {
        double minRatio = Double.MAX_VALUE;
        int res = -1;
        for (int i = 0; i < ratios.getRows(); i++)
        {
            double ratioI = ratios.getElement(i, 0);
            if (minRatio > ratioI && ratioI >= 0) {
                minRatio = ratioI;
                res = i;
            }
        }
        return res;
    }

    /**
     *
     * @param basisVectorsIndices    Indices of the basic vectors (XBi).
     * @param nonBasicVectorsIndices Indices of the non-basic vectors (XPi).
     * @param minRatioRowIndex       The index of the row with the minimum ratio in
     *                               ratios.
     * @param maxDeltaColumnIndex    The index of the column with the maximum value
     *                               in deltas.
     */
    private static void swapIndices(List<Integer> basisVectorsIndices, List<Integer> nonBasicVectorsIndices,
                                    int minRatioRowIndex, int maxDeltaColumnIndex) {
        int a = basisVectorsIndices.get(minRatioRowIndex);
        int b = nonBasicVectorsIndices.get(maxDeltaColumnIndex);
        basisVectorsIndices.set(minRatioRowIndex, b);
        nonBasicVectorsIndices.set(maxDeltaColumnIndex, a);
    }

    /**
     *
     * @param basisCoefficients    Coefficients corresponding to the basis columns
     *                             (Cb).
     * @param nonBasicCoefficients Coefficients corresponding to the non-basic
     *                             columns (Cp).
     * @param minRatioRowIndex     The index of the row with the minimum ratio in
     *                             ratios.
     * @param maxDeltaColumnIndex  The index of the column with the maximum value in
     *                             deltas.
     */
    private static void swapCoefficients(Matrix basisCoefficients, Matrix nonBasicCoefficients,
                                         int minRatioRowIndex, int maxDeltaColumnIndex) {
        double a = basisCoefficients.getElement(minRatioRowIndex, 0);
        double b = nonBasicCoefficients.getElement(maxDeltaColumnIndex, 0);
        basisCoefficients.setElement(minRatioRowIndex, 0, b);
        nonBasicCoefficients.setElement(maxDeltaColumnIndex, 0, a);
    }

    /**
     *
     * @param basisMatrix         A matrix containing a subset of columns from
     *                            constraintsMatrix (basis matrix B).
     * @param nonBasicMatrix      A matrix containing the remaining columns from
     *                            constraintsMatrix (non-basic matrix P).
     * @param minRatioRowIndex    The index of the row with the minimum ratio in
     *                            ratios.
     * @param maxDeltaColumnIndex The index of the column with the maximum value in
     *                            deltas.
     */
    private static void swapColumns(Matrix basisMatrix, Matrix nonBasicMatrix,
                                    int minRatioRowIndex, int maxDeltaColumnIndex) {
        for (int i = 0; i < basisMatrix.getRows(); i++) {
            double a = basisMatrix.getElement(i, minRatioRowIndex);
            double b = nonBasicMatrix.getElement(i, maxDeltaColumnIndex);
            basisMatrix.setElement(i, minRatioRowIndex, b);
            nonBasicMatrix.setElement(i, maxDeltaColumnIndex, a);
        }
    }

    /**
     *
     * @param basisCoefficients   Coefficients corresponding to the basis columns
     *                            (Cb).
     * @param basisVectorsIndices Indices of the basic vectors (XBi).
     * @param basisVectorsValues  Values of the basic variables (XB).
     * @param answer              The result of the objective function evaluation
     *                            (z).
     */
    private static void printAnswer(Matrix basisCoefficients, List<Integer> basisVectorsIndices,
                                    Matrix basisVectorsValues, double answer, int accuracy) {
        HashMap<Integer, Double> map = new HashMap<>();
        for (int i = 0; i < basisCoefficients.getRows(); i++)
        {
            if (basisCoefficients.getElement(i, 0) != 0){

                map.put(basisVectorsIndices.get(i)+1, basisVectorsValues.getElement(i, 0));

            }
        }
        StringBuilder patternBuilder = new StringBuilder("#."); // Start building the pattern

        for (int i = 0; i < accuracy; i++) {
            patternBuilder.append("#"); // Add a placeholder for each decimal place
        }

        DecimalFormat decimalFormat = new DecimalFormat(patternBuilder.toString());
        for (int i = 0; i < map.size(); i++){
            System.out.println("x" + decimalFormat.format(map.keySet().toArray()[i]) + " = " + decimalFormat.format(map.values().toArray()[i]));
        }
        System.out.println("z = " + decimalFormat.format(answer));
    }

    public static void main(String[] args) {

        Scanner in = new Scanner(System.in);

        System.out.println("Is it max or min problem? (0 for min, 1 for max)");
        int typeOfProblem = in.nextInt();
        System.out.println("Number of rows in the matrix of coefficients of constraint function - A:");
        int rows = in.nextInt();
        System.out.println("Number of columns in the matrix of coefficients of constraint function - A:");
        int columns = in.nextInt();
        System.out.println("A vector of coefficients of objective function - C:");
        Matrix coefficients = new Matrix(columns, 1); // c

        for (int i = 0; i < coefficients.getRows(); i++) {
            coefficients.setElement(i, 0, in.nextFloat());
        }

        System.out.println("A matrix of coefficients of constraint function - A:");
        Matrix constraintsMatrix = new Matrix(rows, columns); // A

        for (int i = 0; i < constraintsMatrix.getRows(); i++) {
            for (int j = 0; j < constraintsMatrix.getColumns(); j++) {
                constraintsMatrix.setElement(i, j, in.nextFloat());
            }
        }

        System.out.println("A vector of right-hand side numbers - b:");
        Matrix rightNumbers = new Matrix(rows, 1); // b

        for (int i = 0; i < rightNumbers.getRows(); i++) {
            rightNumbers.setElement(i, 0, in.nextFloat());
        }

        System.out.println("The approximation accuracy (integer number representing the number of decimal places):");
        int accuracy = in.nextInt();

        if(typeOfProblem == 1){
            for (int i = 0; i < coefficients.getRows(); i++){
                if(coefficients.getElement(i, 0) != 0.0){
                coefficients.setElement(i, 0, -1*coefficients.getElement(i, 0));
                }
            }
        }

        Matrix basisMatrix = basisMatrix(constraintsMatrix); // B
        Matrix identity = new Matrix(constraintsMatrix.getRows(), constraintsMatrix.getRows());
        for (int i = 0; i < constraintsMatrix.getRows(); i++) {
            identity.setElement(i, i, 1);
        }
        for (int i = 0; i < basisMatrix.getRows(); i++){
            for (int j = 0; j < basisMatrix.getColumns(); j++){
                if(basisMatrix.getElement(i, j) != identity.getElement(i, j)){
                    System.out.println("The method is not applicable!");
                    return;
                }
            }
        }
        Matrix nonBasicMatrix = nonBasicMatrix(constraintsMatrix); // P
        Matrix inverseBasisMatrix = basisMatrix.inverse(); // B^-1
        Matrix basisCoefficients = basisCoefficients(constraintsMatrix, coefficients); // Cb
        Matrix nonBasicCoefficients = nonBasicCoefficients(constraintsMatrix, coefficients); // Cp
        Matrix basisVectorsValues = Matrix.multiply(inverseBasisMatrix, rightNumbers); // XB
        List<Integer> basisVectorsIndices = basisVectorsIndices(constraintsMatrix); // XBi
        List<Integer> nonBasicVectorsIndices = nonBasisVectorsIndices(constraintsMatrix); // XPi
        double answer = Matrix.multiply(basisCoefficients, basisVectorsValues).getElement(0, 0); // z
        Matrix inverseBasisCoefficients = Matrix.multiply(basisCoefficients, inverseBasisMatrix); // Cb * B^-1

        Matrix deltas = Matrix.multiply(inverseBasisCoefficients, nonBasicMatrix).subtract(nonBasicCoefficients); // zj
        // -
        // Cj
        deltas = deltas.transpose();
        int maxDeltaColumnIndex = maxDeltaColumnIndex(deltas);

        if (maxDeltaColumnIndex == -1) {
            printAnswer(basisCoefficients, basisVectorsIndices, basisVectorsValues, answer, accuracy);
            return;
        }

        Matrix enteringVector = enteringVector(nonBasicMatrix, maxDeltaColumnIndex); // Pe

        Matrix denominatorsForRatios = Matrix.multiply(inverseBasisMatrix, enteringVector); // B^-1 * Pe
        Matrix ratios = calculateRatios(basisVectorsValues, denominatorsForRatios);
        int minRatioRowIndex = minRatioRowIndex(ratios);
        if (minRatioRowIndex == -1){
            System.out.println("The method is not applicable!");
            return;
        }

        while (true) {
            swapIndices(basisVectorsIndices, nonBasicVectorsIndices, minRatioRowIndex, maxDeltaColumnIndex);
            basisCoefficients = basisCoefficients.transpose();
            nonBasicCoefficients = nonBasicCoefficients.transpose();
            swapCoefficients(basisCoefficients, nonBasicCoefficients, minRatioRowIndex, maxDeltaColumnIndex);
            swapColumns(basisMatrix, nonBasicMatrix, minRatioRowIndex, maxDeltaColumnIndex);
            inverseBasisMatrix = basisMatrix.inverse(); // B^-1
            basisVectorsValues = Matrix.multiply(inverseBasisMatrix, rightNumbers); // XB
            basisCoefficients = basisCoefficients.transpose();
            answer = Matrix.multiply(basisCoefficients, basisVectorsValues).getElement(0, 0); // z

            inverseBasisCoefficients = Matrix.multiply(basisCoefficients, inverseBasisMatrix); // Cb * B^-1
            nonBasicCoefficients = nonBasicCoefficients.transpose();
            deltas = Matrix.multiply(inverseBasisCoefficients, nonBasicMatrix).subtract(nonBasicCoefficients); // zj -
            // Cj
            deltas = deltas.transpose();
            maxDeltaColumnIndex = maxDeltaColumnIndex(deltas);

            if (maxDeltaColumnIndex == -1) {
                break;
            }

            enteringVector = enteringVector(nonBasicMatrix, maxDeltaColumnIndex); // Pe

            denominatorsForRatios = Matrix.multiply(inverseBasisMatrix, enteringVector); // B^-1 * Pe
            ratios = calculateRatios(basisVectorsValues, denominatorsForRatios);
            minRatioRowIndex = minRatioRowIndex(ratios);
            if (minRatioRowIndex == -1){
                System.out.println("The method is not applicable!");
                return;
            }

        }
        basisCoefficients = basisCoefficients.transpose();
        if(typeOfProblem == 1){
            answer = -answer;
        }
        printAnswer(basisCoefficients, basisVectorsIndices, basisVectorsValues, answer, accuracy);

        in.close();
    }
}