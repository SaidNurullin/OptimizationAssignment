import java.util.ArrayList;
import java.util.List;

public class Main {

    public static void main(String[] args) {
        List<Double> coefficients = new ArrayList<>();
        List<List<Double>> constraintsMatrix = new ArrayList<>();
        List<Double> rightNumbers = new ArrayList<>();

        List<Double> basisCoefficients = Said.basisCoefficients(constraintsMatrix, coefficients);
        List<Integer> basisVectorsIndices = Ezekiel.basisVectorsIndices(constraintsMatrix);

        List<Double> deltas;
        List<Double> ratios;

        int maxDeltaColumnIndex;
        int minRatioRowIndex;

        while (true){
            deltas = Ezekiel.calculateDeltas(constraintsMatrix, coefficients, basisCoefficients);
            maxDeltaColumnIndex = Timur.maxDeltaColumnIndex(deltas);
            if(maxDeltaColumnIndex == -1){
                break;
            }
            ratios = Said.calculateRatios(maxDeltaColumnIndex, rightNumbers, constraintsMatrix);
            minRatioRowIndex = Daniil.minRatioRowIndex(ratios);
            basisVectorsIndices = Said.newBasisVectorsIndices(basisVectorsIndices, maxDeltaColumnIndex, minRatioRowIndex);
            basisCoefficients = Mikhail.newBasisCoefficients(basisCoefficients, coefficients, maxDeltaColumnIndex, minRatioRowIndex);
            rightNumbers = Daniil.newRightNumbers(rightNumbers, constraintsMatrix, maxDeltaColumnIndex, minRatioRowIndex);
            constraintsMatrix = Said.newConstraintsMatrix(constraintsMatrix, maxDeltaColumnIndex, minRatioRowIndex);
        }

        List<Double> decisionVariables = Timur.decisionVariables(basisVectorsIndices, rightNumbers, coefficients.size());
        double answer = Mikhail.calculateAnswer(rightNumbers, basisCoefficients);

        for (double var : decisionVariables){
            System.out.println(var);
        }
        System.out.println(answer);
    }
}