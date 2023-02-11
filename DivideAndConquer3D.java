import java.util.*;

public class DivideAndConquer3D extends DivideAndConquer {

  int N;
  List<Point> setA;
  List<Point> setB;
  double[][] costMatrix;
  double[][] distance;
  Boolean[] visited;
  double INFINITY;
  int totalCalls;
  int maxLevels;
  long[] time;
  int[] hungarianSearches;
  int[] excess;
  int[] nodes;
  double left, right, top, bottom, front, back;
  int power;
  double threshold;

  /**
   * Constructor
   * @param N - Number of points in setA and setB
   * @param setA - List of points in setA
   * @param setB - List of points in setB
   * @param INFINITY - INFINITY
   * @param costMatrix - N * N costMatrix
   * @param power - Power of the distance
   */
  DivideAndConquer3D(
    int N,
    List<Point> setA,
    List<Point> setB,
    double INFINITY,
    double[][] costMatrix,
    int power
  ) {
    super();
    this.N = N;
    this.setA = setA;
    this.setB = setB;
    this.costMatrix = costMatrix;
    this.INFINITY = INFINITY;
    this.power = power;
    this.threshold = Math.pow(10, (power + 8) * -1);

    // Row 1 will store the shortest distance from source
    // Row 2 will store the parent pointers
    distance = new double[2][2 * N + 7];
    visited = new Boolean[2 * N + 7];

    totalCalls = 0;
    maxLevels = 0;
    time = new long[50];
    hungarianSearches = new int[50];
    excess = new int[50];
    nodes = new int[50];
    for (int i = 0; i < 50; i++) {
      time[i] = 0;
      hungarianSearches[i] = 0;
      excess[i] = 0;
      nodes[i] = 0;
    }
  }

  /**
   * @return - Array of time taken for the conquer step at various levels
   */
  long[] getTimeTakenAtVariousLevels() {
    return time;
  }

  /**
   * @return - Array of nodes at various levels
   */
  int[] getNodesAtVariousLevels() {
    return nodes;
  }

  /**
   * @return - Array of number of hungarian searches done at various levels
   */
  int[] getHungarianSearchesAtVariousLevels() {
    return hungarianSearches;
  }

  /**
   * @return - Total number of hungarian searches done
   */
  int getTotalHungarianSearches() {
    return totalCalls;
  }

  /**
   * @return - Array of excess at various levels
   */
  int[] getExcessAtVariousLevels() {
    return excess;
  }

  /**
   * @return - Recursion tree height
   */
  int getRecursionTreeHeight() {
    return maxLevels;
  }

  /**
   * @param n - Number of points under consideration
   * @return - Index of the point that is not yet visited and has the minimum distance
   */
  int getMinDistanceNode(int n) {
    double minDist = INFINITY;
    int minIdx = -1;
    for (int i = 0; i < n; i++) {
      if (!visited[i] && distance[0][i] < minDist) {
        minDist = distance[0][i];
        minIdx = i;
        if (minDist == 0) {
          break;
        }
      }
    }
    return minIdx;
  }

  /**
   * Augments the matching and updates the dual weights
   * @param A - List of points in set A
   * @param B - List of points in set B
   * @param b - Boundary
   * @throws Exception
   */
  boolean hungarianSearch(List<Point> A, List<Point> B, Boundary b)
    throws Exception {
    int m = A.size();
    int k = m + B.size() + 7;

    // Initialize the distance and visited arrays
    for (int i = 0; i < k; i++) {
      distance[0][i] = INFINITY;
      distance[1][i] = -1;
      visited[i] = false;
    }

    // Distance to the virtual source is 0
    distance[0][0] = 0.0;
    double lMin = INFINITY;
    int freeMinDistanceIdxInA = -1;

    visited[0] = true;
    // Slack of the edge from the virtual source to free points in set B is 0
    for (int u = m + 7; u < k; u++) {
      if (B.get(u - m - 7).getMatchId() == -1) {
        distance[0][u] = 0;
        distance[1][u] = 0;
      }
    }

    // Conduct Dijsktra's until a free point in set A or a boundary point is found
    while (freeMinDistanceIdxInA == -1) {
      int u = getMinDistanceNode(k);
      if (u == -1) {
        throw new Exception(
          "No minimum index found while running the Dijsktra's algorithm"
        );
      }

      // Stop as soon as a free point in set A or a boundary point is found
      if (u < m + 7 && (u < 7 || A.get(u - 7).getMatchId() == -1)) {
        freeMinDistanceIdxInA = u;
        lMin = distance[0][u];
        break;
      }

      // Mark u as visited
      visited[u] = true;

      // Update the distances of the neighbours of u
      if (u < m + 7) {
        // u is a matched point in A and will have only one neighbor
        for (int v = m + 7; v < k; v++) {
          if (A.get(u - 7).getMatchId() == B.get(v - m - 7).getId() - N) {
            double slack =
              costMatrix[A.get(u - 7).getId()][B.get(v - m - 7).getId() - N] -
              A.get(u - 7).getDual() -
              B.get(v - m - 7).getDual();
            if (Math.abs(slack) <= threshold) {
              slack = 0;
            }
            if (distance[0][v] > distance[0][u] + slack) {
              distance[0][v] = distance[0][u] + slack;
              distance[1][v] = u;
            }
            break;
          }
        }
      } else {
        // u is a point of B
        // All the boundaries and unmatched points in A can be reached from u
        double slack =
          Math.pow((b.getTop() - B.get(u - m - 7).getY()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][1] > distance[0][u] + slack) {
          distance[0][1] = distance[0][u] + slack;
          distance[1][1] = u;
        }
        slack =
          Math.pow((b.getRight() - B.get(u - m - 7).getX()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][2] > distance[0][u] + slack) {
          distance[0][2] = distance[0][u] + slack;
          distance[1][2] = u;
        }
        slack =
          Math.pow((B.get(u - m - 7).getY() - b.getBottom()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][3] > distance[0][u] + slack) {
          distance[0][3] = distance[0][u] + slack;
          distance[1][3] = u;
        }
        slack =
          Math.pow((B.get(u - m - 7).getX() - b.getLeft()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][4] > distance[0][u] + slack) {
          distance[0][4] = distance[0][u] + slack;
          distance[1][4] = u;
        }
        slack =
          Math.pow((b.getFront() - B.get(u - m - 7).getZ()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][5] > distance[0][u] + slack) {
          distance[0][5] = distance[0][u] + slack;
          distance[1][5] = u;
        }
        slack =
          Math.pow((B.get(u - m - 7).getZ() - b.getBack()), power) -
          B.get(u - m - 7).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][6] > distance[0][u] + slack) {
          distance[0][6] = distance[0][u] + slack;
          distance[1][6] = u;
        }
        for (int v = 7; v < m + 7; v++) {
          if (B.get(u - m - 7).getMatchId() != A.get(v - 7).getId()) {
            slack =
              costMatrix[A.get(v - 7).getId()][B.get(u - m - 7).getId() - N] -
              A.get(v - 7).getDual() -
              B.get(u - m - 7).getDual();
            if (Math.abs(slack) <= threshold) {
              slack = 0;
            }
            if (distance[0][v] > distance[0][u] + slack) {
              distance[0][v] = distance[0][u] + slack;
              distance[1][v] = u;
            }
          }
        }
      }
    }

    // Get the augmenting path
    List<Integer> path = new ArrayList<>();
    path.add(freeMinDistanceIdxInA);
    while (freeMinDistanceIdxInA != 0) {
      path.add((int) distance[1][freeMinDistanceIdxInA]);
      freeMinDistanceIdxInA = (int) distance[1][freeMinDistanceIdxInA];
    }

    // Update the matching (Augmenting the path)
    for (int i = 0; i < path.size() - 1; i += 2) {
      if (path.get(i) < 7) {
        B.get(path.get(i + 1) - m - 7).setMatchId(path.get(i) - 8);
      } else {
        A
          .get(path.get(i) - 7)
          .setMatchId(B.get(path.get(i + 1) - m - 7).getId() - N);
        B
          .get(path.get(i + 1) - m - 7)
          .setMatchId(A.get(path.get(i) - 7).getId());
      }
    }

    // Update the dual weights
    for (int i = 7; i < k; i++) {
      if (visited[i] && distance[0][i] < lMin) {
        if (i >= m + 7) {
          B
            .get(i - m - 7)
            .setDual(B.get(i - m - 7).getDual() + (lMin - distance[0][i])); // dual weight uodate of points in set B
        } else {
          A
            .get(i - 7)
            .setDual(A.get(i - 7).getDual() - (lMin - distance[0][i])); // dual weight uodate of points in set A
        }
      }
    }
    return true;
  }

  /**
   * Print the matching
   */
  void printMatching() {
    System.out.println("Matching: ");
    Set<Integer> set = new HashSet<>();
    for (int i = 0; i < N; i++) {
      if (setA.get(i).getMatchId() == -1) {
        System.out.println(setA.get(i).toString() + " is umatched");
      } else {
        if (set.contains(setA.get(i).getMatchId())) {
          System.out.println(
            setB.get(setA.get(i).getMatchId()).toString() +
            " matched more than once"
          );
        }
        set.add(setA.get(i).getMatchId());
        System.out.println(
          setA.get(i).toString() +
          " matched to: " +
          setB.get(setA.get(i).getMatchId()).toString()
        );
      }
    }
  }

  /**
   * Prints the dual weights
   */
  void printDualWeights() {
    System.out.println("Dual weights of A: ");
    for (int i = 0; i < N; i++) {
      System.out.println(setA.get(i).getDual());
    }

    System.out.println("Dual weights of B: ");
    for (int i = 0; i < N; i++) {
      System.out.println(setB.get(i).getDual());
    }
  }

  /**
   * Validates the matching
   * @throws Exception
   */
  void validateMatching() throws Exception {
    Set<Integer> set = new HashSet<>();
    for (int i = 0; i < N; i++) {
      if (setA.get(i).getMatchId() == -1) {
        throw new Exception(setA.get(i).toString() + " is umatched");
      } else {
        if (set.contains(setA.get(i).getMatchId())) {
          throw new Exception(
            setB.get(setA.get(i).getMatchId()).toString() +
            " matched more than once"
          );
        }
        set.add(setA.get(i).getMatchId());
      }
    }
    validateDualWeights();
  }

  /**
   * Validates the dual weights
   * @throws Exception
   */
  void validateDualWeights() throws Exception {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (setA.get(i).getMatchId() == N + j) {
          if (
            costMatrix[i][j] -
            setA.get(i).getDual() -
            setB.get(j).getDual() >=
            threshold
          ) {
            System.out.println("Edge cost: " + costMatrix[i][j]);
            System.out.println(
              "Dual Weight of point from set A: " + setA.get(i).getDual()
            );
            System.out.println(
              "Dual Weight of point from set B: " + setB.get(j).getDual()
            );
            System.out.println(
              "Sum of dual weights: " +
              (setA.get(i).getDual() + setB.get(j).getDual())
            );
            System.out.println(
              "Slack: " +
              (costMatrix[i][j] - setA.get(i).getDual() - setB.get(j).getDual())
            );
            throw new Exception(
              "Infeasible dual weights for matched edge: " + i + " -> " + j
            );
          }
        } else if (
          costMatrix[i][j] <
          setA.get(i).getDual() +
          setB.get(j).getDual() -
          threshold
        ) {
          System.out.println("Edge cost: " + costMatrix[i][j]);
          System.out.println(
            "Dual Weight of point from set A: " + setA.get(i).getDual()
          );
          System.out.println(
            "Dual Weight of point from set B: " + setB.get(j).getDual()
          );
          System.out.println(
            "Sum of dual weights: " +
            (setA.get(i).getDual() + setB.get(j).getDual())
          );
          System.out.println(
            "Slack: " +
            (costMatrix[i][j] - setA.get(i).getDual() - setB.get(j).getDual())
          );
          throw new Exception(
            "Infeasible dual weights for unmatched edge: " + i + " -> " + j
          );
        }
      }
    }
  }

  /**
   * @return - Matching cost
   */
  double getMatchingCost() throws Exception {
    int unMatchedPoints = 0;
    double cost = 0;
    for (int i = 0; i < N; i++) {
      if (setA.get(i).getMatchId() == -1) {
        unMatchedPoints++;
      } else {
        cost += costMatrix[i][setA.get(i).getMatchId()];
      }
    }
    if (unMatchedPoints > 0) {
      throw new Exception(
        "Matching has " + unMatchedPoints + "unmatched points"
      );
    }
    return cost;
  }

  /**
   * Solver function
   * @param A - List of points in set A
   * @param B - List of points in set B
   * @param b - Boundary
   * @throws Exception
   */
  void solverHelper(List<Point> A, List<Point> B, Boundary b, int level)
    throws Exception {
    // Update the count of nodes
    nodes[level]++;

    // Update the excess
    excess[level] += Math.abs(A.size() - B.size());

    maxLevels = Math.max(maxLevels, level);
    Boundary bNew = new Boundary(
      b.getTop(),
      b.getBottom(),
      b.getLeft(),
      b.getRight(),
      b.getFront(),
      b.getBack()
    );
    if (b.getLeft() == left) {
      bNew.setLeft(-INFINITY);
    }
    if (b.getRight() == right) {
      bNew.setRight(INFINITY);
    }
    if (b.getTop() == top) {
      bNew.setTop(INFINITY);
    }
    if (b.getBottom() == bottom) {
      bNew.setBottom(-INFINITY);
    }
    if (b.getFront() == front) {
      bNew.setFront(INFINITY);
    }
    if (b.getBack() == back) {
      bNew.setBack(-INFINITY);
    }

    // If no points of set B are present or both left and right or both top and bottom boundaries are the same, do nothing
    if (B.size() == 0) {
      return;
    }

    // If two opposite boundaries coincide, do nothing
    if (
      b.getTop() - b.getBottom() < threshold ||
      b.getRight() - b.getLeft() < threshold ||
      b.getFront() - b.getBack() < threshold
    ) {
      return;
    }

    // If no points of set A are present, update the dual weights of all points in set B to be the shortest distance to the boundary
    if (A.size() == 0) {
      for (int i = 0; i < B.size(); i++) {
        double topDist = Math.pow(bNew.getTop() - B.get(i).getY(), power);
        double rightDist = Math.pow(bNew.getRight() - B.get(i).getX(), power);
        double bottomDist = Math.pow(B.get(i).getY() - bNew.getBottom(), power);
        double leftDist = Math.pow(B.get(i).getX() - bNew.getLeft(), power);
        double frontDist = Math.pow(bNew.getFront() - B.get(i).getZ(), power);
        double backDist = Math.pow(B.get(i).getZ() - bNew.getBack(), power);
        B
          .get(i)
          .setDual(
            Math.min(
              Math.min(
                Math.min(topDist, bottomDist),
                Math.min(rightDist, leftDist)
              ),
              Math.min(frontDist, backDist)
            )
          );
      }
      return;
    }

    // If only one point of set B is present, match it to the nearest point or the boundary
    if (B.size() == 1) {
      double topDist = Math.pow(bNew.getTop() - B.get(0).getY(), power);
      double rightDist = Math.pow(bNew.getRight() - B.get(0).getX(), power);
      double bottomDist = Math.pow(B.get(0).getY() - bNew.getBottom(), power);
      double leftDist = Math.pow(B.get(0).getX() - bNew.getLeft(), power);
      double frontDist = Math.pow(bNew.getFront() - B.get(0).getZ(), power);
      double backDist = Math.pow(B.get(0).getZ() - bNew.getBack(), power);

      double minDist = Math.min(
        Math.min(Math.min(topDist, bottomDist), Math.min(rightDist, leftDist)),
        Math.min(frontDist, backDist)
      );

      int matchedPointIdx = -1;
      for (int i = 0; i < A.size(); i++) {
        if (costMatrix[A.get(i).getId()][B.get(0).getId() - N] <= minDist) {
          minDist = costMatrix[A.get(i).getId()][B.get(0).getId() - N];
          matchedPointIdx = i;
        }
      }
      if (matchedPointIdx != -1) {
        B.get(0).setMatchId(A.get(matchedPointIdx).getId());
        A.get(matchedPointIdx).setMatchId(B.get(0).getId() - N);
      }
      B.get(0).setDual(minDist);
      return;
    }

    // Else split the box into 8 smaller blocks and recursively solve for the smaller subproblems
    // Divide Step
    double xSplit = (b.getLeft() + b.getRight()) / 2;
    double ySplit = (b.getTop() + b.getBottom()) / 2;
    double zSplit = (b.getFront() + b.getBack()) / 2;
    List<Point> A1 = new ArrayList<>();
    List<Point> A2 = new ArrayList<>();
    List<Point> A3 = new ArrayList<>();
    List<Point> A4 = new ArrayList<>();
    List<Point> A5 = new ArrayList<>();
    List<Point> A6 = new ArrayList<>();
    List<Point> A7 = new ArrayList<>();
    List<Point> A8 = new ArrayList<>();
    List<Point> B1 = new ArrayList<>();
    List<Point> B2 = new ArrayList<>();
    List<Point> B3 = new ArrayList<>();
    List<Point> B4 = new ArrayList<>();
    List<Point> B5 = new ArrayList<>();
    List<Point> B6 = new ArrayList<>();
    List<Point> B7 = new ArrayList<>();
    List<Point> B8 = new ArrayList<>();

    // Assign the boxes to each point in set A
    for (Point p : A) {
      if (p.getZ() < zSplit) {
        if (p.getX() < xSplit) {
          if (p.getY() > ySplit) {
            A1.add(p);
          } else {
            A3.add(p);
          }
        } else {
          if (p.getY() > ySplit) {
            A2.add(p);
          } else {
            A4.add(p);
          }
        }
      } else {
        if (p.getX() < xSplit) {
          if (p.getY() > ySplit) {
            A5.add(p);
          } else {
            A7.add(p);
          }
        } else {
          if (p.getY() > ySplit) {
            A6.add(p);
          } else {
            A8.add(p);
          }
        }
      }
    }

    // Assign the boxes to each point in set B
    for (Point p : B) {
      if (p.getZ() < zSplit) {
        if (p.getX() < xSplit) {
          if (p.getY() > ySplit) {
            B1.add(p);
          } else {
            B3.add(p);
          }
        } else {
          if (p.getY() > ySplit) {
            B2.add(p);
          } else {
            B4.add(p);
          }
        }
      } else {
        if (p.getX() < xSplit) {
          if (p.getY() > ySplit) {
            B5.add(p);
          } else {
            B7.add(p);
          }
        } else {
          if (p.getY() > ySplit) {
            B6.add(p);
          } else {
            B8.add(p);
          }
        }
      }
    }

    // Get the new boundaries for the 4 subproblems
    Boundary b1 = new Boundary(
      b.getTop(),
      ySplit,
      b.getLeft(),
      xSplit,
      zSplit,
      b.getBack()
    );
    Boundary b2 = new Boundary(
      b.getTop(),
      ySplit,
      xSplit,
      b.getRight(),
      zSplit,
      b.getBack()
    );
    Boundary b3 = new Boundary(
      ySplit,
      b.getBottom(),
      b.getLeft(),
      xSplit,
      zSplit,
      b.getBack()
    );
    Boundary b4 = new Boundary(
      ySplit,
      b.getBottom(),
      xSplit,
      b.getRight(),
      zSplit,
      b.getBack()
    );
    Boundary b5 = new Boundary(
      b.getTop(),
      ySplit,
      b.getLeft(),
      xSplit,
      b.getFront(),
      zSplit
    );
    Boundary b6 = new Boundary(
      b.getTop(),
      ySplit,
      xSplit,
      b.getRight(),
      b.getFront(),
      zSplit
    );
    Boundary b7 = new Boundary(
      ySplit,
      b.getBottom(),
      b.getLeft(),
      xSplit,
      b.getFront(),
      zSplit
    );
    Boundary b8 = new Boundary(
      ySplit,
      b.getBottom(),
      xSplit,
      b.getRight(),
      b.getFront(),
      zSplit
    );

    // Solve the 8 subproblems independently
    solverHelper(A1, B1, b1, level + 1);
    solverHelper(A2, B2, b2, level + 1);
    solverHelper(A3, B3, b3, level + 1);
    solverHelper(A4, B4, b4, level + 1);
    solverHelper(A5, B5, b5, level + 1);
    solverHelper(A6, B6, b6, level + 1);
    solverHelper(A7, B7, b7, level + 1);
    solverHelper(A8, B8, b8, level + 1);

    // Perform hungarian searches for the points matched to the boundaries
    // Conquer Step
    // Capture the number of hungarian searches and the time taken
    long startTime = System.currentTimeMillis();

    // Get the number of points matched to the boundaries
    int unMatchedPoints = 0;
    for (int i = 0; i < B.size(); i++) {
      if (B.get(i).getMatchId() < 0) {
        B.get(i).setMatchId(-1);
        unMatchedPoints++;
      }
    }

    // Do not perform hungarian searches for points already matched to the current boundaries
    for (int i = 0; i < B.size(); i++) {
      // Do nothing if the point is already matched
      if (B.get(i).getMatchId() >= 0) {
        continue;
      }

      double topDist = Math.pow(bNew.getTop() - B.get(i).getY(), power);
      double rightDist = Math.pow(bNew.getRight() - B.get(i).getX(), power);
      double bottomDist = Math.pow(B.get(i).getY() - bNew.getBottom(), power);
      double leftDist = Math.pow(B.get(i).getX() - bNew.getLeft(), power);
      double frontDist = Math.pow(bNew.getFront() - B.get(0).getZ(), power);
      double backDist = Math.pow(B.get(0).getZ() - bNew.getBack(), power);

      double minDist = Math.min(
        Math.min(Math.min(topDist, bottomDist), Math.min(rightDist, leftDist)),
        Math.min(frontDist, backDist)
      );

      // Match to the boundary if possible
      if (Math.abs(minDist - B.get(i).getDual()) <= threshold) {
        unMatchedPoints--;
        if (minDist == topDist) {
          B.get(i).setMatchId(-7);
        } else if (minDist == rightDist) {
          B.get(i).setMatchId(-6);
        } else if (minDist == bottomDist) {
          B.get(i).setMatchId(-5);
        } else if (minDist == leftDist) {
          B.get(i).setMatchId(-4);
        } else if (minDist == frontDist) {
          B.get(i).setMatchId(-3);
        } else {
          B.get(i).setMatchId(-2);
        }
      }
    }

    // Update the number of hungarian searches performed
    hungarianSearches[level] += unMatchedPoints;
    totalCalls += unMatchedPoints;

    // Perform hungarian search for the points matched to the boundaries
    while (unMatchedPoints > 0) {
      unMatchedPoints--;
      hungarianSearch(A, B, bNew);
    }
    long endTime = System.currentTimeMillis();
    time[level] += endTime - startTime;
  }

  /**
   * @param b - Boundary
   * @return - Minimum cost for perfect matching of points in setA and setB
   */
  double solver(Boundary b) {
    try {
      left = b.getLeft();
      right = b.getRight();
      top = b.getTop();
      bottom = b.getBottom();
      front = b.getFront();
      back = b.getBack();
      solverHelper(setA, setB, b, 0);
      validateMatching();
      return getMatchingCost();
    } catch (Exception e) {
      e.printStackTrace();
      return 0;
    }
  }
}
