import java.util.*;

public class DivideAndConquer {

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
  int[] hungarianSearches; // The number of hungarian searches per level
  int[] matchesToBoundary; // The number of augmenting paths that end at the boundary for each level
  int [] directMatches; // The number of augmenting paths that end at the boundary and are of length 1 for each level
  double[] totalExcess; 
  double [] ExcessSum;
  double [] ExcessCount;
  int pathLength; // Temporary variable that needs global scope in order to get path lengths.
  int [] counter; // Counter for all the levels will help in storing the lengths of all the augmenting paths.
  int [][] totalPathLengths; // A list of all the augmenting paths at each level.
  double left, right, top, bottom;
  int power;
  double threshold;
  int funcCount;
  int boundaryCount;
  int directMatchCount;

  /**
   * Constructor
   * 
   * @param N          - Number of points in setA and setB
   * @param setA       - List of points in setA
   * @param setB       - List of points in setB
   * @param INFINITY   - INFINITY
   * @param costMatrix - N * N costMatrix
   * @param power      - Power of the distance
   */
  DivideAndConquer(
      int N,
      List<Point> setA,
      List<Point> setB,
      double INFINITY,
      double[][] costMatrix,
      int power) {
    this.N = N;
    this.setA = setA;
    this.setB = setB;
    this.costMatrix = costMatrix;
    this.INFINITY = INFINITY;
    this.power = power;
    this.threshold = Math.pow(10, (power + 8) * -1); // threshold = 1/10^(power + 8)

    // Row 1 will store the shortest distance from source
    // Row 2 will store the parent pointers
    distance = new double[2][2 * N + 5]; // 2*N+ 5: Distance from a ficticious source to 2 * N points + 4 boundary
                                         // points + s itself.
    visited = new Boolean[2 * N + 5];
    funcCount = 0;
    boundaryCount = 0;
    totalCalls = 0;
    maxLevels = 0;
    time = new long[50];
    hungarianSearches = new int[50]; // 50 is the max number of levels is what we are expecting for the recursion
                                     // tree to go down to.
    totalExcess = new double[50];
    ExcessSum = new double[50];
    ExcessCount = new double[50];
    directMatches = new int[50];
    matchesToBoundary = new int[50];
    counter = new int[50];
    totalPathLengths = new int[50][10 * N];
    for (int i = 0; i < 50; i++) {
      time[i] = 0;
      hungarianSearches[i] = 0;
      totalExcess[i] = 0;
      ExcessSum[i] = 0.0;
      ExcessCount[i] = 0.0;
      directMatches[i] = 0;
      matchesToBoundary[i] = 0;
      counter[i] = 0;
    }

    for(int i = 0; i < 50; i++)
    {
      for(int j = 0; j < 2 * N; j++)
      {
        totalPathLengths[i][j] = 0;
      }
    }

    

  }

  // Default constructor
  DivideAndConquer() {
  }

  /**
   * @return - Array of time taken for the conquer step at various levels
   */
  long[] getTimeTakenAtVariousLevels() {
    return time;
  }

  /**
   * @return - Array of number of hungarian searches done at various levels
   */
  int[] getHungarianSearchesAtVariousLevels() {
    return hungarianSearches;
  }

 /**
   * @return - Array of number of excess points at various levels
   */
  double[] getTotalExcess() {
    return totalExcess;
  } 

  /**
   * @return - Array of number of augmenting paths who sink is a boundary point at various levels
   */
  int[] getMatchesToBoundary() {
    return matchesToBoundary;
  } 

  /**
   * @return - Array of number of augmenting paths who sink is a boundary point at various levels
   */
  int[] getDirectMatches() {
    return directMatches;
  } 

  /**
   * @return - Array of augmenting path lengths at various levels
   */
  int[][] getTotalPathLengths() {
    return totalPathLengths;
  } 

  /**
   * @return - Total number of hungarian searches done
   */
  int getTotalHungarianSearches() {
    return totalCalls;
  }

  /**
   * @return - Recursion tree height
   */
  int getRecursionTreeHeight() {
    return maxLevels;
  }

  /**
   * @param n - Number of points under consideration
   * @return - Index of the point that is not yet visited and has the minimum
   *         distance
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
   * 
   * @param A - List of points in set A
   * @param B - List of points in set B
   * @param b - Boundary
   * @throws Exception
   */
  boolean hungarianSearch(List<Point> A, List<Point> B, Boundary b) // One iteration of the hungarian algorithm.
      throws Exception {
    funcCount++;
    int m = A.size(); // Number of all the A points
    int k = m + B.size() + 5; // Number of all the A + B points + all 4 boundary points + virtual source S point.
   
    // Initialize the distance and visited arrays
    for (int i = 0; i < k; i++) {
      distance[0][i] = INFINITY; // First row of the distance array stores the distance from the source to one of
                                 // the k points.
      distance[1][i] = -1; // The preceeding point id to the ith id point in the shortest path from s to to
                           // ith point.
      visited[i] = false;
    }
    // Distance to the virtual source is 0
    distance[0][0] = 0.0;
    double lMin = INFINITY; // l_v from the Hungarian Algorithm
    int freeMinDistanceIdxInA = -1; // Point id of the free vertex. It is integer because it is id.
    visited[0] = true;

    // Freeing all the points that are matched to the boundary (for now not differentiating between inner and outer boundarries)
    for (int u = m + 5; u < k; u++) { // Checking to see which B points were initially conncted to the boundary.
        if(B.get(u - m - 5).getMap().containsKey(-1)) // Free all the B points that had positive flow into the inner boundary points.
        {
            distance[0][u] = 0; 
            distance[1][u] = 0;
            double flow = B.get(u - m - 5).getMap().get(-1);
            B.get(u - m - 5).setTotalFlow(B.get(u - m - 5).getTotalFlow() - flow);
            B.get(u - m - 5).getMap().remove(-1);
        }
          // if(b.getTop() == INFINITY && b.getBottom() == -INFINITY && b.getLeft() == -INFINITY && b.getRight() == INFINITY)
          // { 
            // Initially clearing out any B points that match to the boundary.
            if(B.get(u - m - 5).getMap().containsKey(-5)) 
            {
                distance[0][u] = 0; 
                distance[1][u] = 0;
                double flow = B.get(u - m - 5).getMap().get(-5);
                B.get(u - m - 5).setTotalFlow(B.get(u - m - 5).getTotalFlow() - flow);
                B.get(u - m - 5).getMap().remove(-5);
            }

            if(B.get(u - m - 5).getMap().containsKey(-4)) 
            {
                distance[0][u] = 0; 
                distance[1][u] = 0;
                double flow = B.get(u - m - 5).getMap().get(-4);
                B.get(u - m - 5).setTotalFlow(B.get(u - m - 5).getTotalFlow() - flow);
                B.get(u - m - 5).getMap().remove(-4);
            }

            if(B.get(u - m - 5).getMap().containsKey(-3)) 
            {
                distance[0][u] = 0; 
                distance[1][u] = 0;
                double flow = B.get(u - m - 5).getMap().get(-3);
                B.get(u - m - 5).setTotalFlow(B.get(u - m - 5).getTotalFlow() - flow);
                B.get(u - m - 5).getMap().remove(-3);
            }

            if(B.get(u - m - 5).getMap().containsKey(-2))
            {
                distance[0][u] = 0; 
                distance[1][u] = 0;
                double flow = B.get(u - m - 5).getMap().get(-2);
                B.get(u - m - 5).setTotalFlow(B.get(u - m - 5).getTotalFlow() - flow);
                B.get(u - m - 5).getMap().remove(-2);
            } 
          //}
        if ((B.get(u - m - 5).getCapacity() - B.get(u - m - 5).getTotalFlow() - threshold ) > 0) // The free B points. 
        { 
          distance[0][u] = 0; 
          distance[1][u] = 0;
        }
    }
    // Used for index referencing  within cost matrix and slack calculation from Point Id to i.
    HashMap<Integer,Integer> idMap = new HashMap<Integer, Integer>(); //UPDATE HERE
    for(int i = 0; i < B.size(); i++)
    {
      idMap.put(B.get(i).getId(), i); // UPDATE HERE
    }
    // Conduct Dijsktra's until a free point in set A or a boundary point is found
    while (freeMinDistanceIdxInA == -1) 
    {
      int u = getMinDistanceNode(k); // The node id of the l_v node in the hungarian algorithm
      if (u == -1) { // If -1 then all the nodes have been visited.
        throw new Exception(
            "No minimum index found while running the Dijsktra's algorithm");
      }

      // Stop as soon as a free point in setA or a boundary point is found
      if (  (u < 5) || ( (u >= 5 && u < m + 5) && ((A.get(u - 5).getCapacity() - A.get(u - 5).getTotalFlow() - threshold) > 0) ) ) // UPDATE HERE
      {
        freeMinDistanceIdxInA = u;
        lMin = distance[0][u];
        break;
      }

      // Mark u as visited
      visited[u] = true;

      // Update the distances of the neighbours of u
      if (u < m + 5) { // If u is an A point
        // Getting all B points that have a positive flow to "this" A point
        for(int v: A.get(u - 5).getMap().keySet()) { // Finding the matched b point for A point, u.
           if( (A.get(u - 5).getDual() + B.get(idMap.get(v)).getDual()) + threshold  < costMatrix[A.get(u - 5).getId()][v - N] )
           {
                System.out.println("INFEASIBILITY CONDITIONS IN UPDATE NEIGHBOR DISTANCE WHEN U IS A POINT");
                System.out.println("The boundary is " + b.toString());
                System.out.println("a = " + A.get(u - 5));
                System.out.println("b = " + B.get(idMap.get(v)));
                System.out.println("y("+A.get(u - 5).getId()+") = " + A.get(u - 5).getDual());
                System.out.println("y(" + B.get(idMap.get(v)).getId() + ") = " + B.get(idMap.get(v)).getDual());
                System.out.println("c(a,b) = " + costMatrix[A.get(u - 5).getId()][v - N] );
                throw new Exception("INFEASIBLE CONDITION WHEN U IS A POINT");
           }
           double slack = Math.abs(costMatrix[A.get(u - 5).getId()][v - N] - A.get(u - 5).getDual() - B.get(idMap.get(v)).getDual()); // UPDATE HERE
            if (Math.abs(slack) <= threshold) {
              slack = 0;
            }
            if (distance[0][idMap.get(v) + m + 5] > distance[0][u] + slack) { // UPDATE HERE 
              distance[0][idMap.get(v) + m + 5] = distance[0][u] + slack;
              distance[1][idMap.get(v) + m + 5] = u;
            }
        }
      } 
      else 
      {
        // u is a point of B
        // All the boundary points and free points in A (free = capacity - totalFlow > 0) can be reached from u
        double slack = Math.pow((b.getTop() - B.get(u - m - 5).getY()), power) -
            B.get(u - m - 5).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][1] > distance[0][u] + slack) {
          distance[0][1] = distance[0][u] + slack;
          distance[1][1] = u;
        }
        slack = Math.pow((b.getRight() - B.get(u - m - 5).getX()), power) -
            B.get(u - m - 5).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][2] > distance[0][u] + slack) {
          distance[0][2] = distance[0][u] + slack;
          distance[1][2] = u;
        }
        slack = Math.pow((B.get(u - m - 5).getY() - b.getBottom()), power) -
            B.get(u - m - 5).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][3] > distance[0][u] + slack) {
          distance[0][3] = distance[0][u] + slack;
          distance[1][3] = u;
        }
        slack = Math.pow((B.get(u - m - 5).getX() - b.getLeft()), power) -
            B.get(u - m - 5).getDual();
        if (Math.abs(slack) <= threshold) {
          slack = 0;
        }
        if (distance[0][4] > distance[0][u] + slack) {
          distance[0][4] = distance[0][u] + slack;
          distance[1][4] = u;
        }
        for (int i = 5; i < m + 5; i++) {

            double residualCapacity = Math.min(B.get(u - m - 5).getCapacity(), A.get(i - 5).getCapacity());
            if( !B.get(u - m - 5).getMap().containsKey(A.get(i - 5).getId()) || B.get(u - m - 5).getMap().get(A.get(i - 5).getId()) < residualCapacity)
            {
              if(!B.get(u - m - 5).getMap().containsKey(A.get(i - 5).getId()))
              {
                  if( (A.get(i - 5).getDual() + B.get(u - m - 5 ).getDual()) > costMatrix[A.get(i - 5).getId()][B.get(u - m - 5).getId() - N] + threshold) 
                  {
                    System.out.println("INFEASIBLE CONDITION 1 ");
                    System.out.println("The boundary is " + b.toString());
                    System.out.println("a = " + A.get(i - 5));
                    System.out.println("b = " + B.get(u - m - 5 ));
                    System.out.println(" y(" + A.get(i - 5).getId() + ") = " + A.get(i - 5).getDual());
                    System.out.println(" y(" + B.get(u - m - 5 ).getId() + ") = " + B.get(u - m - 5 ).getDual());
                    System.out.println("c(a,b) = " + costMatrix[A.get(i - 5).getId()][B.get(u - m - 5).getId() - N] );
                    throw new Exception("INFEASIBLE CONDITION 1");
                  }
              } 

              else
              {
                  if( (A.get(i - 5).getDual() + B.get(u - m - 5 ).getDual()) + threshold  < costMatrix[A.get(i - 5).getId()][B.get(u - m - 5).getId() - N] ) 
                  {
                    System.out.println("INFEASIBLE CONDITION 2 ");
                    System.out.println("The boundary is " + b.toString());
                    System.out.println(" y(" + A.get(i - 5).getId() + ") = " + A.get(i - 5).getDual());
                    System.out.println(" y(" + B.get(u - m - 5 ).getId() + ") = " + B.get(u - m - 5 ).getDual());
                    System.out.println("c(a,b) = " + costMatrix[A.get(i - 5).getId()][B.get(u - m - 5).getId() - N] );
                    System.out.println("The distance array for A point " + i +  " = " + distance[1][i]);
                    System.out.println("The distance for A vertex " + i + " which is " + A.get(i - 5) + " = " + distance[0][i]);
                    System.out.println("The distance for B vertex " + u + " which is " + B.get(u - m - 5) + " = " + distance[0][u]);
                    throw new Exception("INFEASIBLE CONDITION 2");
                  }
              }
              
              slack = Math.abs(costMatrix[A.get(i - 5).getId()][B.get(u - m - 5).getId() - N] -
                  A.get(i - 5).getDual() -
                  B.get(u - m - 5 ).getDual());
              if (Math.abs(slack) <= threshold) {
                slack = 0;
              }
              
              if (distance[0][i] > distance[0][u] + slack) 
              {
                distance[0][i] = distance[0][u] + slack;
                distance[1][i] = u;
              }
          }
          
        }
      }
    }

    // Get the augmenting path
    List<Integer> path = new ArrayList<>();
    double bottleneckCapacity = INFINITY;
    path.add(freeMinDistanceIdxInA);
    while (freeMinDistanceIdxInA != 0) {
      int parent = (int) distance[1][freeMinDistanceIdxInA];
      if( freeMinDistanceIdxInA >= 5 && freeMinDistanceIdxInA < m + 5) // A point
      { 
        double residualCapacity = Math.min(A.get(freeMinDistanceIdxInA - 5).getCapacity(), B.get(parent - m - 5).getCapacity()); 
        if(A.get(freeMinDistanceIdxInA - 5).getMap().get(B.get(parent  - m - 5).getId()) == null) // If Parent Point is not yet in the neighbor hashmap of the A point
        {
          if(bottleneckCapacity > residualCapacity) 
            bottleneckCapacity = residualCapacity;
        }
        else 
        {
          if(bottleneckCapacity > (residualCapacity - A.get(freeMinDistanceIdxInA - 5).getMap().get(B.get(parent  - m - 5).getId())) ) 
            bottleneckCapacity = (residualCapacity - A.get(freeMinDistanceIdxInA - 5).getMap().get(B.get(parent  - m - 5).getId()));
        }
          
      }
      else if(freeMinDistanceIdxInA >= m + 5) //  B point
      {
        if(parent != 0 && bottleneckCapacity > B.get(freeMinDistanceIdxInA - m - 5).getMap().get(A.get(parent - 5).getId())) 
          bottleneckCapacity = B.get(freeMinDistanceIdxInA - m - 5).getMap().get(A.get(parent - 5).getId()); 
      }
      else // Boundary point
      {
        if(bottleneckCapacity > B.get(parent - m - 5).getCapacity())
          bottleneckCapacity = B.get(parent - m - 5).getCapacity();
      }
      path.add(parent);
      freeMinDistanceIdxInA = parent;
    }
    path.remove(path.size() - 1); // Remove the source vertex from the augmenting path
    pathLength = path.size();

    if(path.get(0) > 0 && path.get(0) < 4)
    {
      boundaryCount++;
      if(path.size() == 2)
        directMatchCount++;
    }
    double endRCapacity  = 0;
    if(path.get(0) > 0 && path.get(0) < 5)  endRCapacity = INFINITY; 
    else endRCapacity = A.get(path.get(0) - 5).getCapacity() - A.get(path.get(0) - 5).getTotalFlow();
    double startRCapacity = B.get(path.get(path.size() - 1) - m - 5).getCapacity() - B.get(path.get(path.size() - 1) - m - 5).getTotalFlow();
    double newFlow = Math.min(bottleneckCapacity, Math.min(startRCapacity, endRCapacity)); // Flow to inject into augmenting path.


    /**
     * 1) Add the flow through the augmenting path (Increase forward Decrease back). 
     * 2) Update the neighbors hashmap of all of these vertices after pushing flow.
     * 3) Update the total flow that is flowing across all of these vertices on augmenting path.
     **/

    // Updating(adding) the total flow of the start and end vertex of the augmenting path.
    if(path.get(0) >= 5)
      A.get(path.get(0) - 5).setTotalFlow(newFlow + A.get(path.get(0) - 5).getTotalFlow()); // changing total flow on the end point.
    B.get(path.get(path.size() - 1) - m - 5).setTotalFlow(newFlow + B.get(path.get(path.size() - 1) - m - 5).getTotalFlow()); // changing the flow on the beginning point.

    for (int i = 0; i < path.size(); i = i + 2) //Iterating from the end of the augmenting path to the start.
    {
      //  The next arrayelement is a B point: the flow from next vertex to this vertex is a forward edge hence increase flow.
      if(i != path.size() - 1 && path.get(i) >= 5)
      {
        // if HashMap does not contain a flow from the previous vertex then add the neighbor and the flow
        if( !A.get(path.get(i)- 5).getMap().containsKey(B.get(path.get(i + 1) - m - 5).getId()))
        {
          A.get(path.get(i) - 5).getMap().put(B.get(path.get(i + 1) - m - 5).getId() , newFlow); 
          B.get(path.get(i + 1) - m - 5).getMap().putIfAbsent(A.get(path.get(i) - 5).getId(), newFlow);
        }
        else
        {
          double oldFlow = A.get(path.get(i) - 5 ).getMap().get(B.get(path.get(i + 1) - m - 5).getId());
          A.get(path.get(i) - 5).getMap().replace(B.get(path.get(i + 1) - m - 5).getId(), oldFlow + newFlow);
          B.get(path.get(i + 1) - m - 5).getMap().replace(A.get(path.get(i) - 5).getId(), oldFlow + newFlow);
        }
      }

      if(path.get(i) < 5)
      {
        B.get(path.get(i + 1) - m - 5).getMap().putIfAbsent(path.get(i) - 6, newFlow);
      }
      
      //The previous arrayelement is a B point: the flow from previous vertex to this vertex is a backward edge hence decrease flow.
      if(i != 0 && path.get(i) >= 5) // Free A vertex in the augmenting path does not have a next node to update if i = 0.
      {
        //if(A.get(path.get(i) - 1).getMap().get(path.get(i - 1)) == null) continue;
        double oldFlow = A.get(path.get(i) - 5).getMap().get(B.get(path.get(i - 1) - m - 5).getId());
        if (oldFlow > newFlow)
        {
          A.get(path.get(i) - 5).getMap().replace(B.get(path.get(i - 1) - m - 5).getId(), oldFlow - newFlow);
          B.get(path.get(i - 1) - m - 5).getMap().replace(A.get(path.get(i) - 5).getId(), oldFlow - newFlow);
        }
        else 
        {
          A.get(path.get(i) - 5).getMap().remove(B.get(path.get(i - 1) - m - 5).getId());
          B.get(path.get(i - 1) - m - 5).getMap().remove(A.get(path.get(i) - 5).getId());
        }
      }
    }
    // Update the dual weights
    for (int i = 5; i < k; i++) {
      if (distance[0][i] < lMin) {
        if (i >= m + 5) {
          B
              .get(i - m - 5)
              .setDual(B.get(i - m - 5).getDual() + (lMin - distance[0][i])); // dual weight update of points in set B
        } else {
          A
              .get(i - 5)
              .setDual(A.get(i - 5).getDual() - (lMin - distance[0][i])); // dual weight update of points in set A
        }
      }
    }

    if ( (B.get(path.get(path.size() - 1) - m - 5).getCapacity() - B.get(path.get(path.size() - 1) - m - 5).getTotalFlow() - threshold> 0))
      return true;
    else return false;
  }

  /**
   * Print the matching
   * @throws Exception
   */
  void printFlow() throws Exception {
    System.out.println("Transport Plan for Divide and Conquer: ");
    for (int i = 0; i < N; i++) {
      if (setA.get(i).getCapacity() - setA.get(i).getTotalFlow() - threshold > 0) {
        throw new Exception (setA.get(i).toString() + " is free");
      }
      if(setB.get(i).getCapacity() - setB.get(i).getTotalFlow() -threshold > 0)
      {
        throw new Exception(setB.get(i).toString() + " is free");
      }
        System.out.println(
          setA.get(i).toString());
        System.out.println(
          setB.get(i).toString());
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
    * 
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
                   " matched more than once");
         }
        set.add(setA.get(i).getMatchId());
     }
     }
     validateDualWeights();
  }

  /**
   * Validates the dual weights
   * 
   * @throws Exception
   */
  void validateDualWeights() throws Exception {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (setA.get(i).getMap().get(N + j) != null) {
          if (costMatrix[i][j] -
              setA.get(i).getDual() -
              setB.get(j).getDual() > threshold) {
            System.out.println("Edge cost: " + costMatrix[i][j]);
            System.out.println(
                "Dual Weight of point from set A: " + setA.get(i).getDual());
            System.out.println(
                "Dual Weight of point from set B: " + setB.get(j).getDual());
            System.out.println(
                "Sum of dual weights: " +
                    (setA.get(i).getDual() + setB.get(j).getDual()));
            System.out.println(
                "Slack: " +
                    (costMatrix[i][j] - setA.get(i).getDual() - setB.get(j).getDual()));
            System.out.println("I am in the if");
            throw new Exception(
                "Infeasible dual weights for matched edge: " + i + " -> " + (N + j) );
          }
        } else if (costMatrix[i][j] < setA.get(i).getDual() +
            setB.get(j).getDual() -
            threshold) {
          System.out.println("I am in the else");
          System.out.println("Edge cost: " + costMatrix[i][j]);
          System.out.println(
              "Dual Weight of point from set A: " + setA.get(i).getDual());
          System.out.println(
              "Dual Weight of point from set B: " + setB.get(j).getDual());
          System.out.println(
              "Sum of dual weights: " +
                  (setA.get(i).getDual() + setB.get(j).getDual()));
          System.out.println(
              "Slack: " +
                  (costMatrix[i][j] - setA.get(i).getDual() - setB.get(j).getDual()));
          throw new Exception(
              "Infeasible dual weights for unmatched edge: " + i + " -> " + (N + j) );
        }
      }
    }
  }

  /**
   * @return - Transport Cost
  **/
  double getTransportCost() throws Exception {
    double cost = 0;
    //double cost1 = 0;
    for (int i = 0; i < setA.size(); i++) {
      for(Map.Entry<Integer,Double> entry: setB.get(i).getMap().entrySet())
      {
        cost += entry.getValue() * costMatrix[entry.getKey()][i];
      }
    }
    return cost;

  }

  /**
   * Solver function
   * 
   * @param A - List of points in set A
   * @param B - List of points in set B
   * @param b - Boundary
   * @throws Exception
   */
  void solverHelper(List<Point> A, List<Point> B, Boundary b, int level)
      throws Exception {
    maxLevels = Math.max(maxLevels, level); // In the first run both levels and maxLevels are 0.
    Boundary bNew = new Boundary(
        b.getTop(),
        b.getBottom(),
        b.getLeft(),
        b.getRight());
    if (b.getLeft() == left) { // if the given boundary is the outermost boundary where all the points live
      bNew.setLeft(-INFINITY); // For calculations regarding matching points b to the boundary.
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

    // If no points of set B are present or both left and right or both top and
    // bottom boundaries are the same, do nothing
    if (B.size() == 0) { // If at this level there are no B points then return.
      return;
    }

    // If two opposite boundaries coincide, do nothing
    if ( // Check to see if the left and right boundaries OR top and bottom boundaries
         // coincide
    b.getTop() - b.getBottom() < threshold ||
        b.getRight() - b.getLeft() < threshold) {
      return; // This is to make sure none of the boundaries are the same size as a point.
    }

    // If no points of set A are present, update the dual weights of all points in
    // set B to be the shortest distance to the boundary
    if (A.size() == 0) {
      for (int i = 0; i < B.size(); i++) { // Get the perpendicular distance for every from the point to the boundary
                                           // line.
        double topDist = Math.pow(bNew.getTop() - B.get(i).getY(), power); // Note Math.pow has a time complexity of
                                                                           // O(1).
        double rightDist = Math.pow(bNew.getRight() - B.get(i).getX(), power);
        double bottomDist = Math.pow(B.get(i).getY() - bNew.getBottom(), power);
        double leftDist = Math.pow(B.get(i).getX() - bNew.getLeft(), power);
        B
            .get(i)
            .setDual(
                Math.min(
                    Math.min(topDist, bottomDist), 
                    Math.min(rightDist, leftDist)));
        B.get(i).getMap().put(-1, B.get(i).getCapacity());
        B.get(i).setTotalFlow(B.get(i).getCapacity());
      }
      return;
    }

    //If only one point of set B is present, match it to the nearest point or the
    //boundary
    // if (B.size() == 1) {
      
    //   double topDist = Math.pow(bNew.getTop() - B.get(0).getY(), power);
    //   double rightDist = Math.pow(bNew.getRight() - B.get(0).getX(), power);
    //   double bottomDist = Math.pow(B.get(0).getY() - bNew.getBottom(), power);
    //   double leftDist = Math.pow(B.get(0).getX() - bNew.getLeft(), power);

    //   double minDist = Math.min(
    //       Math.min(topDist, bottomDist),
    //       Math.min(rightDist, leftDist)); // Getting the minimum distance from the point B to the boundary.
      
    //   for(int i = 0; i < A.size(); i++)
    //   {
    //       if(costMatrix[A.get(i).getId()][B.get(i).getId() - N] < minDist)
    //       {
    //         double startResidualCapacity = B.get(0).getCapacity() - B.get(0).getTotalFlow();
    //         double endResidualCapacity = A.get(i).getCapacity() - A.get(i).getTotalFlow();
    //         double flow = Math.min(startResidualCapacity, endResidualCapacity);
    //         if(flow != 0)
    //         {
    //             A.get(i).setTotalFlow(flow + A.get(i).getTotalFlow());
    //             B.get(i).setTotalFlow(flow + B.get(0).getTotalFlow());
    //             A.get(i).put()
    //         }
            
    //       }
    //   }
      
    // }

    // Else split the box into 4 smaller boxes and recursively solve for the smaller
    // subproblems
    // Divide Step
    double xSplit = (b.getLeft() + b.getRight()) / 2;
    double ySplit = (b.getTop() + b.getBottom()) / 2;
    List<Point> A1 = new ArrayList<>();
    List<Point> A2 = new ArrayList<>();
    List<Point> A3 = new ArrayList<>();
    List<Point> A4 = new ArrayList<>();
    List<Point> B1 = new ArrayList<>();
    List<Point> B2 = new ArrayList<>();
    List<Point> B3 = new ArrayList<>();
    List<Point> B4 = new ArrayList<>();

    // Assign the boxes to each point in set A
    for (Point p : A) {
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
    }

    // Assign the boxes to each point in set B
    for (Point p : B) {
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
    }

    // Get the new boundaries for the 4 subproblems
    Boundary b1 = new Boundary(b.getTop(), ySplit, b.getLeft(), xSplit);
    Boundary b2 = new Boundary(b.getTop(), ySplit, xSplit, b.getRight());
    Boundary b3 = new Boundary(ySplit, b.getBottom(), b.getLeft(), xSplit);
    Boundary b4 = new Boundary(ySplit, b.getBottom(), xSplit, b.getRight());
    // Solve the 4 subproblems independently
    solverHelper(A1, B1, b1, level + 1);
    solverHelper(A2, B2, b2, level + 1);
    solverHelper(A3, B3, b3, level + 1);
    solverHelper(A4, B4, b4, level + 1);
     long startTime = System.currentTimeMillis();
    // Get the number of points matched to the boundaries
    int freePoints = 0;
    for (int i = 0; i < B.size(); i++) {
      if (B.get(i).getMap().containsKey(-1) || ( (B.get(i).getCapacity() - B.get(i).getTotalFlow()) > 0) ) { // If point B is matched to point then match id is the id of that point A
                                       // otherwise it is -1
       // B.get(i).setMatchId(-1); // Seems redundant? - Don't be so sure there is something you need to understand
                                 // first.
        freePoints++;
      }
      if(bNew.getTop() == INFINITY && bNew.getBottom() == -INFINITY && bNew.getLeft() == -INFINITY && bNew.getRight() == INFINITY)
      {
        if (B.get(i).getMap().containsKey(-2) || B.get(i).getMap().containsKey(-3) || B.get(i).getMap().containsKey(-4) || B.get(i).getMap().containsKey(-5))
        {
          freePoints++;
        }
      }
    }
    // Do not perform hungarian searches for points already matched to the current
    // boundaries
    for (int i = 0; i < B.size(); i++) {
      // Do nothing if the point is not free
      if (!B.get(i).getMap().containsKey(-1) || (B.get(i).getCapacity() - B.get(i).getTotalFlow() <= 0) ) { // UPDATE HERE
        continue;
      }

      double topDist = Math.pow(bNew.getTop() - B.get(i).getY(), power); // Get distance of the point to the top
                                                                         // boundary.
      double rightDist = Math.pow(bNew.getRight() - B.get(i).getX(), power);
      double bottomDist = Math.pow(B.get(i).getY() - bNew.getBottom(), power);
      double leftDist = Math.pow(B.get(i).getX() - bNew.getLeft(), power);

      double minDist = Math.min(
          Math.min(topDist, bottomDist),
          Math.min(rightDist, leftDist)); // The actual minimum boundary distance we had figured out before.

      // Match to the boundary if possible
      if (Math.abs(minDist - B.get(i).getDual()) <= threshold) { // Check to see if boundary distance is equal to dual
                                                                 // weight.
        freePoints--;
        if (minDist == topDist) {
          double flow = B.get(i).getMap().get(-1);
          B.get(i).getMap().remove(-1);
          B.get(i).getMap().put(-5, flow);
        } else if (minDist == rightDist) {
          double flow = B.get(i).getMap().get(-1);
          B.get(i).getMap().remove(-1);
          B.get(i).getMap().put(-4, flow);
        } else if (minDist == bottomDist) {
          double flow = B.get(i).getMap().get(-1);
          B.get(i).getMap().remove(-1);
          B.get(i).getMap().put(-3, flow);
        } else {
          double flow = B.get(i).getMap().get(-1);
          B.get(i).getMap().remove(-1);
          B.get(i).getMap().put(-2, flow);
        }
      }
    }
    // Update the number of hungarian searches performed
    hungarianSearches[level] += freePoints;
    totalCalls += freePoints;
   ExcessSum[level] += (B.size() - A.size());
   ExcessCount[level]++;
    while ( freePoints > 0) { //Performing Hungarian untill there are no free points.
      if (!hungarianSearch(A, B, bNew))
        freePoints--;
      totalPathLengths[level][counter[level]] = pathLength;
      counter[level]++;
    }
    matchesToBoundary[level] += boundaryCount;
    directMatches[level] += directMatchCount;
    boundaryCount = 0;
    directMatchCount = 0;
    long endTime = System.currentTimeMillis();
    time[level] += endTime - startTime;
  }

  /**
   * @param b - Boundary
   * @return - Minimum cost for perfect matching of points in setA and setB
   */
  double solver(Boundary b) {
    try {
      left = b.getLeft(); // Gets the left side of the boundary b which is the outermost boundary where
                          // all points live.
      right = b.getRight(); // Gets the right side of the boundary b which is the outermost boundary where
                            // all points live.
      top = b.getTop(); // Gets the top side of the boundary b which is the outermost boundary where all
                        // points live.
      bottom = b.getBottom(); // Gets the bottom side of the boundary b which is the outermost boundary where
                              // all points live.
      solverHelper(setA, setB, b, 0);
      for(int i = 0; i < totalExcess.length; i++)
      {   
          if(ExcessCount[i] != 0)
            totalExcess[i] = (double)(ExcessSum[i] / ExcessCount[i]);
          else 
            totalExcess[i] = 0;  
      }
      //validateMatching();
      validateDualWeights();
      printFlow();
      return getTransportCost();
    } catch (Exception e) {
      e.printStackTrace();
      return 0;
    }
  }
}
