import java.io.*;
import java.util.*;

public class GenerateDatasets {

  static PrintWriter writer;
  static String inputFile;
  static String outputFile;
  static int N;

  public static void generateDatasets() {
    try {
      Scanner scanner = new Scanner(new File(inputFile));
      List<Point> pickupPoints = new ArrayList<>();
      List<Point> dropOffPoints = new ArrayList<>();
      while (scanner.hasNext()) {
        String[] arr = scanner.nextLine().split("\\s+");
        pickupPoints.add(
          new Point(Double.parseDouble(arr[0]), Double.parseDouble(arr[1]), 1)
        );
        dropOffPoints.add(
          new Point(Double.parseDouble(arr[2]), Double.parseDouble(arr[3]), 1)
        );
      }
      scanner.close();
      Collections.shuffle(pickupPoints);
      writer = new PrintWriter(outputFile + "_pickup.txt", "UTF-8");
      for (int i = 0; i < 2 * N; i++) {
        writer.println(
          pickupPoints.get(i).getX() + " " + pickupPoints.get(i).getY()
        );
      }
      writer.close();

      Collections.shuffle(pickupPoints);
      writer = new PrintWriter(outputFile + "_dropoff.txt", "UTF-8");
      for (int i = 0; i < 2 * N; i++) {
        writer.println(
          dropOffPoints.get(i).getX() + " " + dropOffPoints.get(i).getY()
        );
      }
      writer.close();

      Collections.shuffle(pickupPoints);
      writer = new PrintWriter(outputFile + "_pick_and_drop.txt", "UTF-8");
      for (int i = 0; i < N; i++) {
        writer.println(
          pickupPoints.get(i).getX() + " " + pickupPoints.get(i).getY()
        );
      }
      for (int i = 0; i < N; i++) {
        writer.println(
          dropOffPoints.get(i).getX() + " " + dropOffPoints.get(i).getY()
        );
      }
      writer.close();
    } catch (FileNotFoundException | UnsupportedEncodingException e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    String[] months = {
      "jan_01_2014",
      "jan_24_2014",
      "mar_24_2014",
      "apr_13_2014",
      "jun_09_2014",
      "jul_17_2014",
      "aug_11_2014",
    };
    N = 10000;

    for (int i = 0; i < months.length; i++) {
      inputFile = "Datasets/" + months[i] + ".txt";
      outputFile = months[i];
      generateDatasets();
    }
  }
}
