import java.util.*;

public class Traditional {

  int N;
  List<Point> A;
  List<Point> B;
  double[][] costMatrix;
  double[][] distance;
  Boolean[] visited;
  HashMap<Integer, Point> mapIdToPoint;
  double INFINITY;
  double delta = 0.01;
  int power;
  double threshold;
  boolean flag;

  /**
   * Constructor
   * @param N - Number of points in set A and set B
   * @param A - List of points in set A
   * @param B - List of points in set B
   * @param INFINITY - INFINITY
   * @param costMatrix - N * N costMatrix
   * @param power - Power of the distance
   */
  public Traditional(
    int N,
    List<Point> A,
    List<Point> B,
    double INFINITY,
    double[][] costMatrix,
    int power
  ) {
    this.N = N;
    this.A = A;
    this.B = B;
    this.costMatrix = costMatrix;
    this.INFINITY = INFINITY;
    this.power = power;
    this.threshold = Math.pow(10, (power + 8) * -1);
    flag = true;
    // Row 1 will store the shortest distance from source
    // Row 2 will store the parent pointers
    distance = new double[2][2 * N + 1];
    visited = new Boolean[2 * N + 1];
  }

  /**
   * @return - Index of the point that is not yet visited and has the minimum distance
   */
  int getMinDistanceNode() {
    double minDist = INFINITY;
    int minIdx = -1;
    int n = 2 * N + 1;
    for (int i = 0; i < n; i++) {
      if (!visited[i] && distance[0][i] < minDist) {
        minDist = distance[0][i];
        minIdx = i;
      }
    }
    return minIdx;
  }

  /**
   * Augments the matching and updates the dual weights
   * @throws Exception
   */
  boolean hungarianSearch() throws Exception {
    int n = 2 * N + 1;
    if(flag) // Reset Flow values so traditional runs properly.
    {
      flag = false;

      for(int i = 0; i < A.size(); i++)
      {
        A.get(i).setTotalFlow(0);
        A.get(i).getMap().clear();
      }
  
      for(int i = 0; i < B.size(); i++)
      {
        B.get(i).setTotalFlow(0);
        B.get(i).getMap().clear();
      }
    }
   
    // Initialize the distance and visited arrays
    for (int i = 0; i < n; i++) {
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
    for (int u = N + 1; u < n; u++) 
    {

      if ( (B.get(u - N - 1).getCapacity() - B.get(u - N - 1).getTotalFlow()) - threshold > 0) {
        distance[0][u] = 0;
        distance[1][u] = 0;
      }
    }
    // Conduct Dijsktra's until a free point in set A is found
    while (true) 
    {
      int u = getMinDistanceNode();
      if (u == -1) {
        throw new Exception(
          "No minimum index found while running the Dijsktra's algorithm"
        );
      }

      // Mark u as visited
      visited[u] = true;

      // Stop as soon as a free point in set A is found
      if ( (u > 0 && u <= N) && (A.get(u - 1).getCapacity() - A.get(u - 1).getTotalFlow() - threshold) > 0) {
        lMin = distance[0][u];
        freeMinDistanceIdxInA = u;
        break;
      }

      

      // Update the distances of the neighbours of u
      if (u <= N) {
        // u is a non free point in A and all of it's B neighbors that have backedge will have distances that are updated. 
        for(int v: A.get(u - 1).getMap().keySet())
        {
            double slack =
              Math.abs(costMatrix[u - 1][v - N - 1] -
              A.get(u - 1).getDual() -
              B.get(v - N - 1).getDual());
            if (Math.abs(slack) <= threshold) {
              slack = 0;
            }
            if( (A.get(u - 1).getDual() + B.get(v - N - 1).getDual()) + threshold  < costMatrix[A.get(u - 1).getId()][v - N - 1] )
            {
                System.out.println("INFEASIBILITY CONDITIONS IN UPDATE NEIGHBOR DISTANCE WHEN U IS A POINT");
                System.out.println(" y(a) = " + A.get(u - 1).getDual());
                System.out.println(" y(b) = " + B.get(v - N - 1).getDual());
                System.out.println("c(a,b) = " + costMatrix[A.get(u - 1).getId()][v - N - 1] );
                throw new Exception("INFEASIBLE CONDITION WHEN U IS A POINT");
            }
            if (distance[0][v] > distance[0][u] + slack) {
              distance[0][v] = distance[0][u] + slack;
              distance[1][v] = u;
            }
        }
        //System.out.println();
      }
       else {
        // u is a point of B
        // Then the distance to all A points for which we have a forward edge all edges where flow(a,b) < min{sb - da} 
        // All A points not in u's map will have distance updated and all A point in u's map that satisfy the above inequality. 
        for (int i = 1; i <= N; i++) 
        { 
            double residualCapacity = Math.min(B.get(u - N - 1).getCapacity(), A.get(i - 1).getCapacity());
            //System.out.println("Point contains the key " + B.get(u - N - 1).getMap().containsKey(i) );
            if( !B.get(u - N - 1).getMap().containsKey(i) || B.get(u - N - 1).getMap().get(i) < residualCapacity)
            {

              if(!B.get(u - N - 1).getMap().containsKey(i))
              {
                  if( (A.get(i - 1).getDual() + B.get(u - N - 1 ).getDual()) > costMatrix[i - 1][u - N - 1] + threshold) 
                  {
                    System.out.println("INFEASIBLE CONDITION 1 ");
                    System.out.println(" y(a) = " + A.get(i - 1).getDual());
                    System.out.println(" y(b) = " + B.get(u - N - 1 ).getDual());
                    System.out.println("c(a,b) = " + costMatrix[i - 1][u - N - 1] );
                    throw new Exception("INFEASIBLE CONDITION 1");
                  }
              } 

              else
              {
                  if( (A.get(i - 1).getDual() + B.get(u - N- 1 ).getDual()) + threshold  < costMatrix[i - 1][u - N - 1] ) 
                  {
                    System.out.println("INFEASIBLE CONDITION 2 ");
                    System.out.println(" y(a) = " + A.get(i - 1).getDual());
                    System.out.println(" y(b) = " + B.get(u - N - 1 ).getDual());
                    System.out.println("c(a,b) = " + costMatrix[i - 1][u - N - 1] );
                    throw new Exception("INFEASIBLE CONDITION 2");
                  }
              }
              double slack =
              Math.abs(costMatrix[i - 1][u - N - 1] -
              A.get(i - 1).getDual() -
              B.get(u - N - 1).getDual());
              if (Math.abs(slack) <= threshold) 
              {
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

    if (freeMinDistanceIdxInA == -1) {
      throw new Exception("No free vertex found");
    }

    // Get the augmenting path and the bottleneck capacity of the augmenting path. For forward edges need to find min{sb,da} - flow on that edge.
    //For backward edges we want just flow on that edge. Of all those values get the minimum value.
    List<Integer> path = new ArrayList<>();
    double bottleneckCapacity = INFINITY;
    path.add(freeMinDistanceIdxInA);
    while (freeMinDistanceIdxInA != 0) {
      // If the current point is an A point then the edge that connects this vertex to the previous B vertex is a forward edge.
      // Hence we check if the current minimum bottleneck capacity is less than 
      int parent = (int) distance[1][freeMinDistanceIdxInA];
      if(freeMinDistanceIdxInA <= N) // Flow goes from parent to freeMinDistanceIdxInA: parent = B, fMDIIA = A => forward edge.
      {
        double residualCapacity = Math.min(A.get(freeMinDistanceIdxInA - 1).getCapacity(), B.get(parent - N - 1).getCapacity());
        if(A.get(freeMinDistanceIdxInA - 1).getMap().get(parent) == null) 
        {
          if(bottleneckCapacity > residualCapacity) 
            bottleneckCapacity = residualCapacity;
        }
        else 
        {
          if(bottleneckCapacity > (residualCapacity - A.get(freeMinDistanceIdxInA - 1).getMap().get(parent)) )
            bottleneckCapacity = (residualCapacity - A.get(freeMinDistanceIdxInA - 1).getMap().get(parent));
        }
          
      }
      else 
      {
        if(parent != 0 && bottleneckCapacity > B.get(freeMinDistanceIdxInA - N - 1).getMap().get(parent))
          bottleneckCapacity = B.get(freeMinDistanceIdxInA - N - 1).getMap().get(parent);
      }
      path.add(parent);
      freeMinDistanceIdxInA = parent;
    }
    path.remove(path.size() - 1); // Remove the source vertex from the augmenting path

    double endRCapacity = A.get(path.get(0) - 1).getCapacity() - A.get(path.get(0) - 1).getTotalFlow();
    //System.out.println("This is the end residual capacity " + endRCapacity);
    double startRCapacity = B.get(path.get(path.size() - 1) - N - 1).getCapacity() - B.get(path.get(path.size() - 1) - N - 1).getTotalFlow();
    //System.out.println("This is the the start residual capacity " + startRCapacity);
    double newFlow = Math.min(bottleneckCapacity, Math.min(startRCapacity, endRCapacity)); // Flow to inject into augmenting path.
    //System.out.println("This is the new flow that is being augmented into this path " + newFlow);
    /**
     * 1) Add the flow through the augmenting path (Increase forward Decrease back). 
     * 2) Update the neighbors hashmap of all of these vertices after pushing flow.
     * 3) Update the total flow that is flowing across all of these vertices on augmenting path.
     **/

    // Updating(adding) the total flow of the start and end vertex of the augmenting path.
    A.get(path.get(0) - 1).setTotalFlow(newFlow + A.get(path.get(0) - 1).getTotalFlow()); // changing total flow on the end point.
    B.get(path.get(path.size() - 1) - N - 1).setTotalFlow(newFlow + B.get(path.get(path.size() - 1) - N - 1).getTotalFlow()); // changing the flow on the beginning point.

    
    for (int i = 0; i < path.size(); i = i + 2) //Iterating from the end of the augmenting path to the start.
    {
      //  The next arrayelement is a B point: the flow from next vertex to this vertex is a forward edge hence increase flow.
      if(i != path.size() - 1)
      {
        // if HashMap does not contain a flow from the previous vertex then add the neighbor and the flow
        if( !A.get(path.get(i)- 1).getMap().containsKey(path.get(i + 1)))
        {
          A.get(path.get(i) - 1).getMap().put(path.get(i + 1) , newFlow); 
          B.get(path.get(i + 1) - N - 1).getMap().put(path.get(i), newFlow);
        }
        else
        {
          double oldFlow = A.get(path.get(i) - 1 ).getMap().get(path.get(i + 1));
          A.get(path.get(i) - 1).getMap().replace(path.get(i + 1), oldFlow + newFlow);
          B.get(path.get(i + 1) - N - 1).getMap().replace(path.get(i), oldFlow + newFlow);
        }
      }
      
      //The previous arrayelement is a B point: the flow from previous vertex to this vertex is a backward edge hence decrease flow.
      if(i != 0) // Free A vertex in the augmenting path does not have a next node to update.
      {
        //if(A.get(path.get(i) - 1).getMap().get(path.get(i - 1)) == null) continue;
        double oldFlow = A.get(path.get(i) - 1).getMap().get(path.get(i - 1));
        if (oldFlow > newFlow)
        {
          A.get(path.get(i) - 1).getMap().replace(path.get(i - 1), oldFlow - newFlow);
          B.get(path.get(i - 1) - N - 1).getMap().replace(path.get(i), oldFlow - newFlow);
        }
        else 
        {
          A.get(path.get(i) - 1).getMap().remove(path.get(i - 1));
          //System.out.println("This is the point " + B.get(path.get(i - 1) - N - 1).getMap().get(path.get(i)));
          B.get(path.get(i - 1) - N - 1).getMap().remove(path.get(i));
        }
      }
    }

    // Update the dual weights
    for (int i = 1; i < n; i++) {
      if (visited[i] && distance[0][i] < lMin) {
        if (i > N) {
          B
            .get(i - N - 1)
            .setDual(B.get(i - N - 1).getDual() + (lMin - distance[0][i])); // dual weight update of points in set B
        } else {
          A
            .get(i - 1)
            .setDual(A.get(i - 1).getDual() - (lMin - distance[0][i])); // dual weight update of points in set A
        }
      }
      //System.out.println("The point1 is " + A.get(path.get(0)).getId() + " with total flow1 of " + A.get(path.get(0)).getTotalFlow());

    }

    //System.out.println();

        if( (B.get(path.get(path.size() - 1) - N - 1).getCapacity() - B.get(path.get(path.size() - 1) - N - 1).getTotalFlow()) - threshold > 0)
          return true;
        else return false;
  }






  /**
   * Prints the matching
   * @throws Exception
   **/
   
  void printFlow() throws Exception {
    System.out.println("Transport Plan for Traditional: ");
    for (int i = 0; i < N; i++) 
    {
      if (A.get(i).getCapacity() - A.get(i).getTotalFlow() - threshold > 0) {
        System.out.println(A.get(i).toString() + " is free");
      }
      if(A.get(i).getCapacity() < A.get(i).getTotalFlow())
        throw new Exception("Capacity of A point " + A.get(i) + " is greater than the TotalFlow");
      if(B.get(i).getCapacity() - B.get(i).getTotalFlow() -threshold > 0)
      {
        System.out.println(B.get(i).toString() + " is free");
      }
      if(B.get(i).getCapacity() < B.get(i).getTotalFlow())
        throw new Exception("Capacity of B point" + B.get(i) + " is greater than the TotalFlow");
        System.out.println(
          A.get(i).toString());
    }
  }
  

  /**
   * Validates the flow. Check if it is a max flow. 
   * @throws Exception
  
   
  void validateFlow () throws Exception {
    Set<Integer> set = new HashSet<>();
    for (int i = 0; i < N; i++) {
      if (A.get(i).getMatchId() == -1) {
        throw new Exception(A.get(i).toString() + " is umatched");
      } else {
        if (set.contains(A.get(i).getMatchId())) {
          throw new Exception(
            B.get(A.get(i).getMatchId()).toString() + " matched more than once"
          );
        }
        set.add(A.get(i).getMatchId());
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
        if (A.get(i).getMap().get(N + j + 1) != null) {
          if (
            costMatrix[i][j] -
            A.get(i).getDual() -
            B.get(j).getDual() >
            threshold
          ) {
            System.out.println("Edge cost: " + costMatrix[i][j]);
            System.out.println(
              "Dual Weight of point from set A: " + A.get(i).getDual()
            );
            System.out.println(
              "Dual Weight of point from set B: " + B.get(j).getDual()
            );
            System.out.println(
              "Sum of dual weights: " +
              (A.get(i).getDual() + B.get(j).getDual())
            );
            System.out.println(
              "Slack: " +
              (costMatrix[i][j] - A.get(i).getDual() - B.get(j).getDual())
            );
            throw new Exception(
              "Infeasible dual weights for matched edge: " + i + " -> " + j
            );
          }
        } else if (
          costMatrix[i][j] < A.get(i).getDual() + B.get(j).getDual() - threshold
        ) {
          System.out.println("Edge cost: " + costMatrix[i][j]);
          System.out.println(
            "Dual Weight of point from set A: " + A.get(i).getDual()
          );
          System.out.println(
            "Dual Weight of point from set B: " + B.get(j).getDual()
          );
          System.out.println(
            "Sum of dual weights: " + (A.get(i).getDual() + B.get(j).getDual())
          );
          System.out.println(
            "Slack: " +
            (costMatrix[i][j] - A.get(i).getDual() - B.get(j).getDual())
          );
          throw new Exception(
            "Infeasible dual weights for unmatched edge: " + i + " -> " + j
          );
        }
      }
    }
  }

  /**
   * @return - Transport Cost
  **/
  double getTransportCost() throws Exception {
    double cost = 0;
    for (int i = 0; i < N; i++) {
      for(Map.Entry<Integer,Double> entry: B.get(i).getMap().entrySet())
      {
        cost += entry.getValue() * costMatrix[entry.getKey() - 1][i];
      }
    }
    return cost;
  }

  /**
   * @return - Minimum cost for perfect matching of points in setA and setB
   **/
  
  double solver() {
    try {
      int freePoints = N;
      while (freePoints > 0) {
        if(!hungarianSearch())
          freePoints--;
      }
      validateDualWeights();
      //printFlow();
      return getTransportCost();
    } catch (Exception e) {
      e.printStackTrace();
      return 0;
    }
  }
} 
