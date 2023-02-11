import java.util.*;
public class Point {

  private double x, y, z, dual, capacity, totalFlow;
  private int id, matchId;
  private HashMap<Integer,Double> neighbours;

  Point(double x, double y, double capacity) {
    this.x = x;
    this.y = y;
    this.z = 0;
    this.dual = 0.0;
    this.id = -1;
    this.capacity = capacity;
    neighbours = new HashMap<Integer, Double>(); //Map of neighbors(key) of this point with their respective positive flow(value).
    this.totalFlow = 0.0;
    this.matchId = -1;
  }

  Point(double x, double y, double z, double capacity) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.dual = 0.0;
    this.id = -1;
    this.capacity = capacity;
    neighbours = new HashMap<Integer, Double>(); //Map of neighbors(key) of this point with their respective positive flow(value).
    this.totalFlow = 0.0;
  }

  public double getX() {
    return this.x;
  }

  public void setX(double x) {
    this.x = x;
  }

  public double getY() {
    return this.y;
  }

  public void setY(double y) {
    this.y = y;
  }

  public double getZ() {
    return this.z;
  }

  public void setZ(double z) {
    this.z = z;
  }

  public double getDual() {
    return this.dual;
  }

  public void setDual(double dual) {
    this.dual = dual;
  }

  public int getId() {
    return this.id;
  }

  public void setId(int id) {
    this.id = id;
  }

  public double getCapacity() {
    return this.capacity;
  }

  public void setCapacity(double capacity) {
    this.capacity = capacity;
  }

  public double getTotalFlow() {
    return this.totalFlow;
  }

  public void setTotalFlow(double totalFlow) {
    this.totalFlow = totalFlow;
  }

  public HashMap<Integer, Double> getMap()
  {
    return neighbours;
  }

  public int getMatchId() {
    return this.matchId;
  }

  public void setMatchId(int matchId) {
    this.matchId = matchId;
  }
  @Override
  public String toString() {
    return String.format(
      "X: " +
      this.x +
      " Y: " +
      this.y +
      " Z: " +
      this.z +
      " id: " +
      this.id +
      " capacity: " +
      this.capacity +
      " totalFlow: " +
      this.totalFlow +
      " neighbors map "+
      this.neighbours
      +" dual Weights " +
      this.getDual()
    );
  }
}
