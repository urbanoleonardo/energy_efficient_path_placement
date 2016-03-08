import java.util.ArrayList;

public class Trajectory {
	
/*
 * Trajectory class.
 * It will contain a list (or a vector) of targets and will be generated from the method
 * INTERPOLATE inside the class path.
 * trajectory should have inside methods to retrieve some points or to shift the whole trajectory 
 * of a certain vector.
 */
	public ArrayList<Target> points = new ArrayList<Target>();
	public ArrayList<Double> timeInstants = new ArrayList<Double>();
	public ArrayList<double[]> extTorques = new ArrayList<double[]>();
	public ArrayList<double[]> extForces = new ArrayList<double[]>();
	
	public void shiftTime(double timeOffset){
		for(int i = 0 ; i < timeInstants.size(); i++){
			timeInstants.set(i, (timeInstants.get(i) + timeOffset));
			
		}
	}
	
}
