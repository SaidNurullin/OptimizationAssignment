import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class Main {
    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the constraints.
     * @return A matrix containing a subset of columns from constraintsMatrix that construct identity matrix.
     */
    private static Matrix basisMatrix(Matrix constraintsMatrix) {
        return null;
    }

    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the constraints.
     * @return A matrix containing the columns from constraintsMatrix that not construct identity matrix
     */
    private static Matrix nonBasicMatrix(Matrix constraintsMatrix) {
        return null;
    }

    /**
     *
     * @param constraintsMatrix A matrix representing the coefficients of the constraints.
     * @param coefficients A matrix representing the coefficients of the objective function.
     * @return Coefficients corresponding to the basis columns
     */
    private static Matrix basisCoefficients(Matrix constraintsMatrix, Matrix coefficients) {
        return null;
    }

    /**
     *
     * @param coefficients A matrix representing the coefficients of the objective function.
     * @param basisCoefficients Coefficients corresponding to the basis columns (Cb).
     * @return Coefficients corresponding to the non-basic columns.
     */
    private static Matrix nonBasicCoefficients(Matrix coefficients, Matrix basisCoefficients) {
        return null;
    }

    /**
     *
     * @param nonBasicMatrix A matrix containing the remaining columns from constraintsMatrix (non-basic matrix P).
     * @param maxDeltaColumnIndex The index of the column with the maximum value in deltas.
     * @return The column vector from the non-basic matrix corresponding to maxDeltaColumnIndex.
     */
    private static Matrix enteringVector(Matrix nonBasicMatrix, int maxDeltaColumnIndex) {
        Matrix enteringVector = new Matrix(nonBasicMatrix.getRows(), 1);
        for (int i = 0; i < enteringVector.getRows(); i++){
            enteringVector.setElement(i, 0, nonBasicMatrix.getElement(i, maxDeltaColumnIndex));
        }

        return enteringVector;
    }

    /**
     *
     * @param basisVectorsValues Values of the basic variables (XB).
     * @param denominatorsForRatios Values to which basisVectorsValues should be divided
     * @return Vector of ratios (basisVectorsValues)/(denominatorsForRatios).
     */
    private static Matrix calculateRatios(Matrix basisVectorsValues, Matrix denominatorsForRatios) {
        Matrix ratios = new Matrix(basisVectorsValues.getRows(), 1);
        for (int i = 0; i < ratios.getRows(); i++){
            ratios.setElement(i, 0, basisVectorsValues.getElement(i, 0) /
                    denominatorsForRatios.getElement(i, 0));
        }
        return ratios;
    }

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
     * @param constraintsMatrix  A matrix representing the coefficients of the constraints.
     * @return Array List of indices of columns that construct the basis in the constraints matrix as in the example above.
     */
    private static List<Integer> basisVectorsIndices(Matrix constraintsMatrix) {
        return null;
    }

    /**
     * Example:<br>
     * constraintsMatrix =<br>
     * 1 2 3 0 4 0<br>
     * 0 6 4 0 1 1<br>
     * 0 3 2 1 6 0<br>
     * non-basic columns =<br>
     * 2 3 4<br>
     * 6 4 1<br>
     * 3 2 6<br>
     * indices: {1, 2, 4}
     * @param constraintsMatrix  A matrix representing the coefficients of the constraints.
     * @return Array List of indices of columns that not construct the basis in the constraints matrix as in the example above.
     */
    private static List<Integer> nonBasisVectorsIndices(Matrix constraintsMatrix) {
        return null;
    }

    /**
     *
     * @param deltas Array List of doubles, each element is delta calculated for each column from constraints matrix
     * @return the index of maximum positive value among deltas, if no positive deltas, returns -1.
     */
    private static int maxDeltaColumnIndex(Matrix deltas) {
        return 0;
    }

    /**
     *
     * @param ratios Vector of ratios
     * @return The index of the row with the minimum ratio in ratios.
     */
    private static int minRatioRowIndex(Matrix ratios) {
        return 0;
    }

    /**
     *
     * @param basisVectorsIndices Indices of the basic vectors (XBi).
     * @param nonBasicVectorsIndices Indices of the non-basic vectors (XPi).
     * @param minRatioRowIndex The index of the row with the minimum ratio in ratios.
     * @param maxDeltaColumnIndex The index of the column with the maximum value in deltas.
     */
    private static void swapIndices(List<Integer> basisVectorsIndices, List<Integer> nonBasicVectorsIndices,
                                    int minRatioRowIndex, int maxDeltaColumnIndex) {
    }

    /**
     *
     * @param basisCoefficients Coefficients corresponding to the basis columns (Cb).
     * @param nonBasicCoefficients Coefficients corresponding to the non-basic columns (Cp).
     * @param minRatioRowIndex The index of the row with the minimum ratio in ratios.
     * @param maxDeltaColumnIndex The index of the column with the maximum value in deltas.
     */
    private static void swapCoefficients(Matrix basisCoefficients, Matrix nonBasicCoefficients,
                                         int minRatioRowIndex, int maxDeltaColumnIndex) {
    }

    /**
     *
     * @param basisMatrix A matrix containing a subset of columns from constraintsMatrix (basis matrix B).
     * @param nonBasicMatrix A matrix containing the remaining columns from constraintsMatrix (non-basic matrix P).
     * @param minRatioRowIndex The index of the row with the minimum ratio in ratios.
     * @param maxDeltaColumnIndex The index of the column with the maximum value in deltas.
     */
    private static void swapColumns(Matrix basisMatrix, Matrix nonBasicMatrix,
                                    int minRatioRowIndex, int maxDeltaColumnIndex) {
    }

    /**
     *
     * @param basisCoefficients Coefficients corresponding to the basis columns (Cb).
     * @param basisVectorsIndices Indices of the basic vectors (XBi).
     * @param basisVectorsValues Values of the basic variables (XB).
     * @param answer The result of the objective function evaluation (z).
     */
    private static void printAnswer(Matrix basisCoefficients, List<Integer> basisVectorsIndices,
                                    Matrix basisVectorsValues, double answer) {
    }

    public static void main(String[] args) {

        Scanner in = new Scanner(System.in);



        Matrix coefficients = new Matrix(6, 1); //c

        for (int i = 0; i < coefficients.getRows(); i++){
            coefficients.setElement(i, 0, in.nextFloat());
        }


        Matrix constraintsMatrix = new Matrix(4,6); //A

        for (int i = 0; i < constraintsMatrix.getRows(); i++){
            for (int j = 0; j < constraintsMatrix.getColumns(); j++){
                constraintsMatrix.setElement(i, j, in.nextFloat());
            }
        }

        Matrix rightNumbers = new Matrix(4, 1); //b

        for (int i = 0; i < rightNumbers.getRows(); i++){
            rightNumbers.setElement(i, 0, in.nextFloat());
        }





        Matrix basisMatrix = basisMatrix(constraintsMatrix); //B
        Matrix nonBasicMatrix = nonBasicMatrix(constraintsMatrix); // P
        Matrix inverseBasisMatrix = basisMatrix.inverse(); // B^-1
        Matrix basisCoefficients = basisCoefficients(constraintsMatrix, coefficients); //Cb
        Matrix nonBasicCoefficients = nonBasicCoefficients(coefficients, basisCoefficients); //Cp
        Matrix basisVectorsValues = Matrix.multiply(inverseBasisMatrix, rightNumbers); //XB
        List<Integer> basisVectorsIndices = basisVectorsIndices(constraintsMatrix); //XBi
        List<Integer> nonBasicVectorsIndices = nonBasisVectorsIndices(constraintsMatrix); // XPi
        double answer = Matrix.multiply(basisCoefficients, basisVectorsValues).getElement(0, 0); //z
        Matrix inverseBasisCoefficients = Matrix.multiply(basisCoefficients, inverseBasisMatrix); //Cb * B^-1

        Matrix deltas = Matrix.multiply(inverseBasisCoefficients, nonBasicMatrix).subtract(nonBasicCoefficients); // zj - Cj
        int maxDeltaColumnIndex = maxDeltaColumnIndex(deltas);

        if(maxDeltaColumnIndex == -1){
            printAnswer(basisCoefficients, basisVectorsIndices, basisVectorsValues, answer);
        }

        Matrix enteringVector = enteringVector(nonBasicMatrix, maxDeltaColumnIndex); //Pe

        Matrix denominatorsForRatios = Matrix.multiply(inverseBasisMatrix, enteringVector); //B^-1 * Pe
        Matrix ratios = calculateRatios(basisVectorsValues, denominatorsForRatios);
        int minRatioRowIndex = minRatioRowIndex(ratios);




        while (true){
            swapIndices(basisVectorsIndices, nonBasicVectorsIndices, minRatioRowIndex, maxDeltaColumnIndex);
            swapCoefficients(basisCoefficients, nonBasicCoefficients, minRatioRowIndex, maxDeltaColumnIndex);
            swapColumns(basisMatrix, nonBasicMatrix, minRatioRowIndex, maxDeltaColumnIndex);
            inverseBasisMatrix = basisMatrix.inverse(); // B^-1
            basisVectorsValues = Matrix.multiply(inverseBasisMatrix, rightNumbers); //XB
            answer = Matrix.multiply(basisCoefficients, basisVectorsValues).getElement(0,0); //z


            inverseBasisCoefficients = Matrix.multiply(basisCoefficients, inverseBasisMatrix); //Cb * B^-1
            deltas = Matrix.multiply(inverseBasisCoefficients, nonBasicMatrix).subtract(nonBasicCoefficients); // zj - Cj
            maxDeltaColumnIndex = maxDeltaColumnIndex(deltas);

            if(maxDeltaColumnIndex == -1){
                break;
            }

            enteringVector = enteringVector(nonBasicMatrix, maxDeltaColumnIndex); //Pe


            denominatorsForRatios = Matrix.multiply(inverseBasisMatrix, enteringVector); //B^-1 * Pe
            ratios = calculateRatios(basisVectorsValues, denominatorsForRatios);
            minRatioRowIndex = minRatioRowIndex(ratios);

        }

        printAnswer(basisCoefficients, basisVectorsIndices, basisVectorsValues, answer);


    }




}