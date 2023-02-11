import java.io.*;
import java.util.*;

public class Runner {

  static int N; // 
  static List<Point> A;
  static List<Point> B;
  static double[][] costMatrix;
  static final double xMin = 0;
  static final double xMax = 10;
  static final double yMin = 0;
  static final double yMax = 10;
  static final double zMin = 0;
  static final double zMax = 10;
  static int power = 2;
  static double INFINITY = Double.MAX_VALUE;
  static Boundary b = new Boundary(yMax, yMin, xMin, xMax, zMax, zMin); // We work with a 3D space even when we we have 2D poinst.
  static final double delta = 0.01;
  static String distanceType = "Manhattan"; // The distance between two points measured along axes at right angles. 
                                            //Distance =|x_2-x_1| + |y_2-y_1| + |z_2-z_1| between points (x_1,y_1,z_1), (x_2,y_2,z_2)
  static boolean generateIntegerPoints;
  static boolean generate3DPoints;  // Have the choice to generate 3D points instead of 2D points, 
  static boolean runOnlyDivideAndConquer;
  static int totalCalls;
  static int maxLevels;
  static long[] time;
  static int[] hungarianSearches;
  static double [] totalExcess;
  static int [] matchesToBoundary;
  static int [] directMatches;
  static int [][] totalPathLengths;
  static boolean generateRandomPoints = true;
  static String inputFile;

  /**
   * @return - A random point within the given region
   */
  static Point getRandomPoint() {
    double xValue = xMin + Math.random() * (xMax - xMin); // Generating random x coordinates of points between xMin and xMax.
    double yValue = yMin + Math.random() * (yMax - yMin); // Generating random y coordinates of points between yMin and yMax.
    if (generate3DPoints) { // Have a choice as to whether we want to generate 3D or 2D points. 
      double zValue = zMin + Math.random() * (zMax - zMin);
      return new Point(xValue, yValue, zValue);
    }
    //Random rand = new Random();
    //return new Point(xValue, yValue, 1.0);
    //System.out.println(Math.abs(rand.nextInt()));
    return new Point(xValue, yValue, Math.abs(0.0 + Math.random() * (10.0 - 0.0)));
  }

  /**
   * @return - A random point with integer coordinates within the given region
   */
  static Point getRandomIntegerPoint() {
    Random rand = new Random();
    double xValue = rand.nextInt((int) (xMax - xMin)) + xMin; // Gives a pseudorandom, uniformly distributed int value [xMin,xMax).
    double yValue = rand.nextInt((int) (yMax - yMin)) + yMin; // Gives a pseudorandom, uniformly distributed int value [yMin,yMax).
    if (generate3DPoints) {
      double zValue = rand.nextInt((int) (zMax - zMin)) + zMin;
      return new Point(xValue, yValue, zValue);
    }
    return new Point(xValue, yValue, 1);
  }

  /**
   * @param a - Point 1
   * @param b - Point 2
   * @return - Distance between the two points
   */
  static double getDistance(Point a, Point b) {
    if (distanceType.equals("Manhattan")) {
      return (
        Math.pow(
          (
            Math.abs(a.getX() - b.getX()) +
            Math.abs(a.getY() - b.getY()) +
            Math.abs(a.getZ() - b.getZ()) // Will return 0 if we are working with 2D points. 
          ),
          power // Math expression represented above: (|x_2-x_1| + |y_2-y_1| + |z_2-z_1|) ^ 3 for 2D

        )
      );
    } else { // Math expression represented below: (|x_2-x_1| + |y_2-y_1| + |z_2-z_1|)^(power/2) if power is even
             // Math expression represented below: (|x_2-x_1| + |y_2-y_1| + |z_2-z_1|)^( (power + 1) /2) if power is odd. 
      return (
        Math.pow(
          (
            ((a.getX() - b.getX()) * (a.getX() - b.getX())) +
            ((a.getY() - b.getY()) * (a.getY() - b.getY())) +
            ((a.getZ() - b.getZ()) * (a.getZ() - b.getZ()))
          ),
          power / 2 
        ) *
        (
          power % 2 == 1
            ? Math.sqrt(
              ((a.getX() - b.getX()) * (a.getX() - b.getX())) +
              ((a.getY() - b.getY()) * (a.getY() - b.getY())) +
              ((a.getZ() - b.getZ()) * (a.getZ() - b.getZ()))
            )
            : 1
        )
      );
    }
  }

  static double[] generateRandomCapacities(int sum, int count)
  {
    java.util.Random g = new java.util.Random();

    double vals[] = new double[count];
    sum -= count;

    for (int i = 0; i < count-1; ++i) {
        vals[i] = g.nextInt(sum);
    }
    vals[count-1] = sum;

    java.util.Arrays.sort(vals);
    for (int i = count-1; i > 0; --i) {
        vals[i] -= vals[i-1];
    }
    for (int i = 0; i < count; ++i) { ++vals[i]; }

    
    return vals;

  }

  static boolean checkCapacities(List<Point> A, List<Point> B)
  {
      for(int i = 0; i < A.size(); i++)
      {
        System.out.println(A.get(i));
        System.out.println(B.get(i));
        if( A.get(i).getCapacity() < 0)
        {
          System.out.println(A.get(i));
          return false;
        }
        if( B.get(i).getCapacity() < 0)
        {
          System.out.println(B.get(i));
          return false;
        }
      }
      return true;
  }

  /**
   * @return - Matching cost returned by OT calculation
   */
  static double getOptimalDistanceCost() {
    double[] supplies = new double[N];
    double[] demands = new double[N];
    for (int i = 0; i < N; i++) {
      supplies[i] = B.get(i).getCapacity();
      demands[i] = A.get(i).getCapacity();
    }
    Mapping mp = new Mapping(N, demands, supplies, costMatrix, delta);
    mp.verifyFlow(demands, supplies, mp.getFlow()); // Ensure that proposed matching is a feasible one.
    return mp.getTotalCost(); // Total cost of the matching.
  }

  /**
   * Prints the cost matrix
   */
  static void printCostMatrix() {
    if (costMatrix == null) {
      System.out.println("Cost Matrix is not yet initialized");
      return;
    }
    System.out.println("Cost Matrix:");
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        System.out.print(costMatrix[i][j] + " ");
      }
      //System.out.println();
    }
  }

  static void run() throws Exception{
    A = new ArrayList<>();
    B = new ArrayList<>();
    costMatrix = new double[N][N];
    Point p;
    double supplies[] = generateRandomCapacities(3000, 10);
    double demand[] = generateRandomCapacities(3000, 10);

    if (generateRandomPoints) {
      // Generate n points for the set A
      Set<Point> points = new HashSet<>();
      double ASum = 0;
      double BSum = 0;
      for (int i = 0; i < N; i++) {
        if (generateIntegerPoints) {
          p = getRandomIntegerPoint(); //Whether we need to test with random integer points or we need to test with double points
        } else {
          p = getRandomPoint();
        }
        if (points.contains(p)) { // Make sure that the points stored in set A are distinct, and we account for any accidental
                                  // points that could be the same.
          i--;
        } else {
          p.setId(i);
          ASum = ASum + p.getCapacity();
          A.add(p);
          points.add(p);
        } // Note that when getRandomPoint or getRandomIntegerPoint those generated points have a dual weight of 0.
      }

      // Generate n points for the set B
      for (int i = 0; i < N; i++) {
        if (generateIntegerPoints) {
          p = getRandomIntegerPoint();
        } else {
          p = getRandomPoint();
        }
        if (points.contains(p)) {
          i--;
        } else {
          p.setId(N + i);
          BSum = BSum + p.getCapacity();
          B.add(p);
          points.add(p);
        }
      }
      for(int i = 0; i < A.size(); i++)
      {
        A.get(i).setCapacity(1);
        B.get(i).setCapacity(1);
      }
      // Double sum1 = 0d;
      // Double sum2 = 0d;
      // for(int i = 0; i < A.size(); i++)
      // {
      //   A.get(i).setCapacity( (A.get(i).getCapacity())/ASum);
      //   //sum1 = sum1 + A.get(i).getCapacity();
      //   if(i == A.size() - 1)
      //   {
      //     System.out.println("sum1 = " + sum1);
      //     A.get(i).setCapacity( 1.0 - sum1);
      //     //A.get(i).setCapacity( (A.get(i).getCapacity())/ASum);
      //     sum1 = sum1 + A.get(i).getCapacity();
      //     System.out.println("A point: " + A.get(i));
      //     System.out.println("sum1 = " + sum1);
      //   }
      //   else
      //   {
      //     sum1 = sum1 + A.get(i).getCapacity();
      //   }
      //   B.get(i).setCapacity( (B.get(i).getCapacity())/BSum);
      //   if(i == B.size() - 1)
      //   {
      //     B.get(i).setCapacity( 1.0 - sum2);
      //    // B.get(i).setCapacity( (B.get(i).getCapacity())/BSum);
      //     sum2 = sum2 + B.get(i).getCapacity();
      //     System.out.println("B point: " + B.get(i));
      //     System.out.println("sum2 = " + sum2);
      //   }
      //   else
      //     sum2 = sum2 + B.get(i).getCapacity();
      // }
      // if(!checkCapacities(A, B))
      //   throw new Exception(" Negative Capacity Error " );
      // if(sum1 != 1.0)
      //   throw new Exception("The sum of all A capacities = " + sum1 );
      // if(sum2 != 1.0)
      //   throw new Exception("The sum of all B capacities = " + sum2 );
      // System.out.println();
    } else {
      try {
        Scanner scanner = new Scanner(new File(inputFile));
        for (int i = 0; i < N; i++) { //Adding points in set A 
          String[] arr = scanner.nextLine().split("\\s+"); // Splitting units by one or more whitespace.
          p = new Point(Double.parseDouble(arr[0]), Double.parseDouble(arr[1]), 1);
          p.setId(i);
          A.add(p);
        }
        for (int i = 0; i < N; i++) { //Adding points in set B.
          String[] arr = scanner.nextLine().split("\\s+");
          p = new Point(Double.parseDouble(arr[0]), Double.parseDouble(arr[1]), 1);
          p.setId(N + i);
          B.add(p);
        }
        scanner.close();
      } catch (FileNotFoundException e) {
        e.printStackTrace(); 
      }
    }

    //System.out.println("THE COST MATRIX GIVEN TO DIVIDE AND CONQUER METHOD");
    // Initialize the cost matrix
    for (int i = 0; i < N; i++) { // Filling the adjacency matrix with the cost of every edge. This is how we define our graph.
      for (int j = 0; j < N; j++) {
        costMatrix[i][j] = getDistance(A.get(i), B.get(j)); //The distance between given point A and given point B 
        //System.out.println("The cost of [" + i + "][" + (N+j) +"] = " + costMatrix[i][j] );

      }
    }
   // System.out.println();

    long startTime, endTime;
    double matchingCost1, matchingCost2;
    DivideAndConquer dac;

    startTime = System.currentTimeMillis();
    if (generate3DPoints) {
      dac = new DivideAndConquer3D(N, A, B, INFINITY, costMatrix, power); //S
    } else {
      dac = new DivideAndConquer(N, A, B, INFINITY, costMatrix, power);
    }
    matchingCost2 = dac.solver(b);
    System.out.println("Divide and Conquer cost: " + matchingCost2);
    endTime = System.currentTimeMillis();
    System.out.println(
      "Time taken by divide and conquer approach: " + (endTime - startTime)
    );
    maxLevels = dac.getRecursionTreeHeight();
    time = dac.getTimeTakenAtVariousLevels();
    hungarianSearches = dac.getHungarianSearchesAtVariousLevels();
    totalCalls = dac.getTotalHungarianSearches();
    totalExcess = dac.getTotalExcess();
    matchesToBoundary = dac.getMatchesToBoundary();
    directMatches = dac.getDirectMatches();
    totalPathLengths = dac.getTotalPathLengths();
    for(int i = 0; i < maxLevels; i++)
    {
       Arrays.sort(totalPathLengths[i], 0, hungarianSearches[i]);
    }


    if (!runOnlyDivideAndConquer) {
      for (int i = 0; i < N; i++) {
        A.get(i).setDual(0.0);
        A.get(i).setMatchId(-1);
        B.get(i).setDual(0.0);
        B.get(i).setMatchId(-1);
      }
      startTime = System.currentTimeMillis();
      Traditional trad = new Traditional(N, A, B, INFINITY, costMatrix, power);
      matchingCost1 = trad.solver();
      System.out.println("Traditional cost: " + matchingCost1);
      endTime = System.currentTimeMillis();
      System.out.println(
        "Time taken by traditional approach: " + (endTime - startTime)
      );

      startTime = System.currentTimeMillis();
      double otCost = getOptimalDistanceCost();
      System.out.println("OT Cost: " + otCost);
      endTime = System.currentTimeMillis();
      System.out.println("Time taken by OT: " + (endTime - startTime));

      if (matchingCost1 - otCost > 0) {
        System.out.println("Traditional Matching is not optimal");
      }
      if (matchingCost2 - otCost > 0) {
        System.out.println("Divide and Conquer Matching is not optimal");
      }
      if (
        Math.abs(matchingCost1 - matchingCost2) > Math.pow(10, (power + 8) * -1)
      ) {
        System.out.println("--------------------------------------------");
        System.out.println("--------------------------------------------");
        System.out.println("Matching costs do not match");
        System.out.println("--------------------------------------------");
        System.out.println("--------------------------------------------");
      }
    }
    System.out.println("Recursion tree height: " + maxLevels);
    System.out.println("Average number of excess points (excess = B - A) at various levels:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.print(totalExcess[i] + " ");
    }
    System.out.println("\nAugmenting path lengths searches done at various levels in the order max, min and median:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.println("Level " + i);
      int max = (hungarianSearches[i] - 1) > 0 ? hungarianSearches[i] - 1 : 0;
      System.out.println(totalPathLengths[i][max] + " " + totalPathLengths[i][0] + " " + totalPathLengths[i][ (int)(0.5 * (hungarianSearches[i])) ] );
    }
    System.out.println("Hungarian searches done at various levels:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.print(hungarianSearches[i] + " ");
    }
    System.out.println("\nNumber of Augmenting paths which end at the boundary at various levels:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.print(matchesToBoundary[i] + " ");
    }
    System.out.println("\nNumber of Augmenting paths of length 1 that directly match to the boundary:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.print(directMatches[i] + " ");
    }
    System.out.println("\nTime taken at various levels:");
    for (int i = 0; i <= maxLevels; i++) {
      System.out.print(time[i] + " ");
    }
    System.out.println("\nTotal hungarian searches: " + totalCalls);
  }

  public static void main(String[] args) throws Exception {
    boolean isDebugMode = false;
    if (isDebugMode) {
      runOnlyDivideAndConquer = false;
      distanceType = "Euclidean";
      generate3DPoints = false;
      generateIntegerPoints = false;
      power = 2;
      N = 1000;
      for (int i = 0; i < 15; i++) {
        run();
        System.out.println();
        System.out.println();
      }
      return;
    }
    if (generateRandomPoints) {
      distanceType = args[0];
      generateIntegerPoints = Boolean.parseBoolean(args[1]);
      generate3DPoints = Boolean.parseBoolean(args[2]);
      int minValue = Integer.parseInt(args[3]);
      int maxValue = Integer.parseInt(args[4]);
      int increment = Integer.parseInt(args[5]);
      int numOfRuns = Integer.parseInt(args[6]);
      runOnlyDivideAndConquer = Boolean.parseBoolean(args[7]);
      power = Integer.parseInt(args[8]);
      System.out.println("Distance Type: " + distanceType);
      System.out.println(
        "Points type: " + (generateIntegerPoints ? "Integer" : "Double")
      );
      System.out.println("Dimensions: " + (generate3DPoints ? 3 : 2));
      System.out.println();
      System.out.println("---------------------------------------------");
      for (N = minValue; N <= maxValue; N += increment) {
        System.out.println("Points size: " + N);
        System.out.println("---------------------------------------------");
        System.out.println();
        for (int j = 0; j < numOfRuns; j++) {
          run();
          System.out.println();
        }
        System.out.println("---------------------------------------------");
      }
    } else {
      N = 1000;
      INFINITY = Double.MAX_VALUE;
      b.setLeft(-75.0);
      b.setRight(-72.0);
      b.setTop(42.0);
      b.setBottom(38.0);
      runOnlyDivideAndConquer = true;
      String[] months = {
        "jan_01_2014",
        "jan_24_2014",
        "mar_24_2014",
        "apr_13_2014",
        "jun_09_2014",
        "jul_17_2014",
        "aug_11_2014",
      };
      distanceType = args[1];
      int start = Integer.parseInt(args[2]);
      int end = Integer.parseInt(args[3]);
      power = Integer.parseInt(args[4]);
      for (int j = start; j < end; j++) {
        String month = months[j];
        System.out.println(
          "--------------------Month: " + month + "--------------------"
        );
        System.out.println();
        if (args[0].charAt(0) == 'P') {
          System.out.println("--------------------Pickup--------------------");
          inputFile = "Datasets/10000/" + month + "_pickup.txt";
          for (int i = 0; i < 5; i++) {
            run();
            System.out.println("");
          }
        } else if (args[0].charAt(0) == 'D') {
          System.out.println("--------------------Dropoff--------------------");
          inputFile = "Datasets/10000/" + month + "_dropoff.txt";
          for (int i = 0; i < 5; i++) {
            run();
            System.out.println("");
          }
        } else {
          System.out.println(
            "--------------------Pickup and Dropoff--------------------"
          );
          inputFile = "Datasets/10000/" + month + "_pick_and_drop.txt";
          for (int i = 0; i < 5; i++) {
            run();
            System.out.println("");
          }
        }
      }
    }
  }
}
