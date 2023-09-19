import java.util.ArrayList;
import java.util.List;

public class Main {

    public static List<Double> basisCoefficients(List<List<Double>> constraintsMatrix, List<Double> coefficients){
        return null;
    }

    public static List<Integer> basisVectorsIndices(List<List<Double>> constraintsMatrix){
        return null;
    }

    public static List<Double> calculateDeltas(List<List<Double>> constraintsMatrix,
                                        List<Double> coefficients, List<Double> basisCoefficients){
        return null;
    }

    public static int maxDeltaColumnIndex(List<Double> deltas){
        return -1;
    }

    public static List<Double> calculateRatios(int maxDeltaColumnIndex, List<Double> rightNumbers,
                                        List<List<Double>> constraintsMatrix){
        return null;
    }

    public static int minRatioRowIndex(List<Double> ratios){
        return 0;
    }

    public static List<Integer> newBasisVectorsIndices(List<Integer> basisVectorIndices,
                                                int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    public static List<Double> newBasisCoefficients(List<Double> basisCoefficients, List<Double> coefficients,
                                             int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    public static List<Double> newRightNumbers(List<Double> rightNumbers, List<List<Double>> constraintsMatrix,
                                        int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    public static List<List<Double>> newConstraintsMatrix(List<List<Double>> constraintsMatrix,
                                                   int maxDeltaColumnIndex, int minRatioRowIndex){
        return null;
    }

    public static List<Double> decisionVariables(List<Integer> basisVectorsIndices, List<Double> rightNumbers){
        return null;
    }

    public static double calculateAnswer(List<Double> rightNumbers, List<Double> basisCoefficients){
        return 0;
    }

    public static void main(String[] args) {
        List<Double> coefficients = new ArrayList<>();
        List<List<Double>> constraintsMatrix = new ArrayList<>();
        List<Double> rightNumbers = new ArrayList<>();

        List<Double> basisCoefficients = basisCoefficients(constraintsMatrix, coefficients);
        List<Integer> basisVectorsIndices = basisVectorsIndices(constraintsMatrix);

        List<Double> deltas;
        List<Double> ratios;

        int maxDeltaColumnIndex;
        int minRatioRowIndex;

        while (true){
            deltas = calculateDeltas(constraintsMatrix, coefficients, basisCoefficients);
            maxDeltaColumnIndex = maxDeltaColumnIndex(deltas);
            if(maxDeltaColumnIndex == -1){
                break;
            }
            ratios = calculateRatios(maxDeltaColumnIndex, rightNumbers, constraintsMatrix);
            minRatioRowIndex = minRatioRowIndex(ratios);
            basisVectorsIndices = newBasisVectorsIndices(basisVectorsIndices, maxDeltaColumnIndex, minRatioRowIndex);
            basisCoefficients = newBasisCoefficients(basisCoefficients, coefficients, maxDeltaColumnIndex, minRatioRowIndex);
            rightNumbers = newRightNumbers(rightNumbers, constraintsMatrix, maxDeltaColumnIndex, minRatioRowIndex);
            constraintsMatrix = newConstraintsMatrix(constraintsMatrix, maxDeltaColumnIndex, minRatioRowIndex);
        }

        List<Double> decisionVariables = decisionVariables(basisVectorsIndices, rightNumbers);
        double answer = calculateAnswer(rightNumbers, basisCoefficients);

        for (double var : decisionVariables){
            System.out.println(var);
        }
        System.out.println(answer);
    }
}