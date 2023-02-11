import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class CostMatrixParser {

  double[][] parseCostMatrix(String inputFile) {
    Scanner scanner;
    double[][] costMatrix;
    try {
      scanner = new Scanner(new File(inputFile));
      String[] arr = scanner.nextLine().split("\\s+");
      int rows = Integer.parseInt(arr[0]);
      int cols = Integer.parseInt(arr[1]);
      costMatrix = new double[rows][cols];
      for (int i = 0; i < rows; i++) {
        arr = scanner.nextLine().split("\\s+");
        for (int j = 0; j < cols; j++) {
          costMatrix[i][j] = Double.parseDouble(arr[j]);
        }
      }
      return costMatrix;
    } catch (FileNotFoundException e) {
      System.out.println("Error while parsing the cost matrix");
      e.printStackTrace();
      costMatrix = new double[0][0];
      return costMatrix;
    }
  }

  ArrayList<ArrayList<Double>> parseCostMatrix1(String inputFile) {
    Scanner scanner;
    ArrayList<ArrayList<Double>> costMatrix = new ArrayList<>();
    String[] arr;
    try {
      scanner = new Scanner(new File(inputFile));
      while (scanner.hasNextLine()) {
        arr = scanner.nextLine().split("\\s+");
        ArrayList<Double> temp = new ArrayList<>();
        for (int i = 0; i < arr.length; i++) {
          temp.add(Double.parseDouble(arr[i]));
        }
        costMatrix.add(temp);
      }
      return costMatrix;
    } catch (FileNotFoundException e) {
      System.out.println("Error while parsing the cost matrix");
      e.printStackTrace();
      return costMatrix;
    }
  }
}
