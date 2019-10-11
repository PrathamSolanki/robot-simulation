package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class DefaultTeam {
  public ArrayList<Point> calculAngularTSP(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
    //REMOVE >>>>>
    ArrayList<Point> result = (ArrayList<Point>)hitPoints.clone();
    //<<<<< REMOVE

    int N = points.size();

    double dist[][] = new double[N][N];
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++) {
        double distance = points.get(i).distance(points.get(j));
        if (distance <= edgeThreshold) dist[i][j] = distance;
        else dist[i][j] = Double.MAX_VALUE;
      }

    for (int p = 0; p < N; p++) 
      for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) {
          if (dist[i][p] + dist[p][j] < dist[i][j]) 
            dist[i][j] = dist[i][p] + dist[p][j]; 
        }

    System.out.println(dist[51][13]);

    return result;
  }
}
