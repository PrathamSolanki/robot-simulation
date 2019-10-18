package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.TreeMap;

public class DefaultTeam {
	double matrix[][];
	int path[][];
	int startingPointIndex;
	
  public ArrayList<Point> calculAngularTSP(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {

	  TreeMap<Double, ArrayList<Point>> solution = new TreeMap<>();
	  
	  for (int k=0; k < hitPoints.size(); k++) {
		  startingPointIndex = k;
		  ArrayList<Point> result = new ArrayList<>();
		    
		    ArrayList<Point> redPoints = calculTSP(hitPoints);
		    
		    computeMatrixDistances(points, edgeThreshold);
		    
		    
//			Point hitPoint = redPoints.get(0);
			Point startPoint = redPoints.get(0);
			
			int redPointsSize = redPoints.size();
			
			for (int i=0; i < redPointsSize; i++) {
				Point hitPoint = redPoints.get(i);
				int index = points.indexOf(hitPoint);
				
//				double minDistance = Double.MAX_VALUE; 
//				Point closestPoint = redPoints.get(i+1);
//				for (Point p : redPoints) {
//					if (p.equals(hitPoint)) continue;
//					int index2 = points.indexOf(p);
////					System.out.println(String.valueOf(index) + ", " + String.valueOf(index2));
//					double currentDistance = matrix[index][index2];
//					if (currentDistance < minDistance) {
//						minDistance = currentDistance;
//						closestPoint = p;
//					}
//				}
				
				ArrayList<Point> shortestPath;
				if (i == redPointsSize-1) {
					shortestPath = getShortestPath(points, index, points.indexOf(startPoint));
				}
				else {
					Point closestPoint = redPoints.get(i+1);
					shortestPath = getShortestPath(points, index, points.indexOf(closestPoint));
				}
				
				result.addAll(shortestPath);
				
//				redPoints.remove(hitPoint);
//				hitPoint = closestPoint;	
			}
			
			solution.put(angleScore(result), result);
	  }
    
	  return solution.firstEntry().getValue();
  }
  
  public static double angleScore(ArrayList<Point> inpts){
	    ArrayList<Point> pts = (ArrayList<Point>)inpts.clone();
	    for (int i=1;i<pts.size()+1;i++){
	      if (pts.get(i-1).equals(pts.get(i%pts.size()))){
	        pts.remove(i%pts.size());
	        i--;
	      }
	    }
	    if (pts.size()<=1) return 0;
	    double result = angle(pts.get(pts.size()-1),pts.get(0),pts.get(1))+pts.get(pts.size()-1).distance(pts.get(0));
	    for (int i=0;i<pts.size()-1;i++) {
	      result+=pts.get(i).distance(pts.get(i+1))+angle(pts.get(i),pts.get(i+1),pts.get((i+2)%pts.size()));
	    }
	    return result;
	  }
	  public static double angle(Point p, Point q, Point r){
	    if (p.distance(q)*q.distance(r)==0) return 0;
	    double dot = ((q.getX()-p.getX()) * (r.getX()-q.getX()) + (q.getY()-p.getY()) * (r.getY()-q.getY()));
	    double unitDot = dot/ (p.distance(q)*q.distance(r));
	    if (unitDot<-1) unitDot=-1;
	    if (unitDot>1) unitDot=1;
	    double result = (100/Math.PI)*Math.acos(unitDot);
	    return result;
	  }
  
  
  private ArrayList<Point> calculTSP(ArrayList<Point> points) {
		if (points.size()<4) {
          return points;
      }

      ArrayList<Point> return_points = new ArrayList<Point>();
      ArrayList<Point> rest = new ArrayList<Point>();
      rest.addAll(points);
      
      Point point = rest.remove(startingPointIndex);
      return_points.add(point);
      
      
      while (!rest.isEmpty()) {
      	double min_dist = Double.MAX_VALUE;
      	int nearest_point = 0;
      	
      	for(int i=0; i < rest.size(); i++) {
      		if(rest.get(i).distance(point) < min_dist) {
      			nearest_point = i;
      			min_dist = rest.get(i).distance(point);
      		}
      	}
      	
      	point = rest.remove(nearest_point);
      	return_points.add(point);    	
      }
      
      // Now let's improve the solution by removing the overlapping paths
      // Local optimization
      double current_score = score(return_points);
      double old_score = current_score +1;
 
      
      while (old_score > current_score) {
      	return_points = improve(return_points);
      	
      	old_score = current_score;
      	current_score = score(return_points);
      }

      /******************
       * TO BE MODIFIED *
       ******************/
      return return_points;
	  }

  
  private double score(ArrayList<Point> points) {
  	double score = 0;
  	for (int i=1; i<points.size(); i++) {
  		score += points.get(i).distance(points.get(i-1));
  	}
  	
  	score += points.get(points.size()-1).distance(points.get(0));
	
  	return score;
  }
  
  private ArrayList<Point> improve (ArrayList<Point> points) {
	  for (int i=0;i<points.size();i++){
          for (int j=i+2;j<points.size() ;j++){
              double a=points.get(i).distance(points.get((i+1)%points.size()));
              double b=points.get(j%points.size()).distance(points.get((j+1)%points.size()));
              double c=points.get(i).distance(points.get(j%points.size()));
              double d=points.get((i+1)%points.size()).distance(points.get((j+1)%points.size()));
              if (a+b>c+d) {
                  ArrayList<Point> p=new ArrayList<Point>();
                  for (int k=0;k<=i;k++) p.add(points.get(k));
                  for (int k=j;k>i;k--) p.add(points.get(k));
                  for (int k=j+1;k<points.size();k++) p.add(points.get(k));
                  return p;
              }
          }
      }
      return points;
  }
  
  
  private ArrayList<Point> getShortestPath(ArrayList<Point> points, int index1, int index2) {
	  
	  ArrayList<Point> result = new ArrayList<>();
	  int nextIndex = index1;
	  while (nextIndex != index2) {
		  result.add(points.get(nextIndex));
		  nextIndex = path[nextIndex][index2];
	  }
	  
	  return result;
  }
  
  
  private void computeMatrixDistances(ArrayList<Point> points, int edgeThreshold) {
	  int N = points.size();

	    matrix = new double[N][N];
	    path = new int[N][N];

	    for (int i = 0; i < N; i++) 
	      for (int j = 0; j < N; j++) {
	        double distance = points.get(i).distance(points.get(j));
	        if (distance <= edgeThreshold) {
	          matrix[i][j] = distance;
	          path[i][j] = j;
	        }
	        else {
	          matrix[i][j] = Double.MAX_VALUE;
	          path[i][j] = -1;
	        }
	      }

	    for (int p = 0; p < N; p++) 
	      for (int i = 0; i < N; i++) 
	        for (int j = 0; j < N; j++) {
	          if (matrix[i][p] < Double.MAX_VALUE && matrix[p][j] < Double.MAX_VALUE && matrix[i][p] + matrix[p][j] < matrix[i][j]) { 
	            matrix[i][j] = matrix[i][p] + matrix[p][j];
	            path[i][j] = path[i][p];
	          }
	        }
  }
}
