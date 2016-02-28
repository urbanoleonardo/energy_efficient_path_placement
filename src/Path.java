import com.sun.javafx.geom.Quat4f;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Quat4d;

public class Path {
/*
 * Path class that contains initial and final TARGETS, and then has a Trajectory attribute
 * that will be filled if the INTERPOLATE method will be called (still have to see if it will 
 * automatically fill the attribute or whatnot).
 * The "Father" method INTERPOLATE will generate basically a straight line
 * while for other paths like the square one (that will be sons) the method will behave differently 
 */
	private Target initialPosition;
	private Target finalPosition;
	private Trajectory interpolatedPath; //the Trajectory will have also the initial and 
										  // final position within of course.
										  	
	
	private double tSample; //time resolution for the trajectory
	private double xSample; //space resolution for the trajectory
	
	private double maxVel; //these two parameters MAYBE have to belong to the robot
	private double maxAcc;
	
	//Different constructors depending on how many parameters are specified
	
	public Path(Target initialPosition, Target finalPosition, double tSample, double xSample){
		this.initialPosition = initialPosition;
		this.finalPosition = finalPosition;
		this.tSample = tSample;
		this.xSample = xSample;
	}
	
	public Path(Target initialPosition, Target finalPosition)
	{
		this.initialPosition = initialPosition;
		this.finalPosition = finalPosition;
	}
	//END part of constructor
	
	//Part of GET methods
	public Target[] getPathPositions()
	{
		
		/*
		 * It returns a vector of Target with initial and final position
		 */
		
		Target[] pathPositions = new Target[2];
		pathPositions[0] = this.initialPosition;
		pathPositions[1] = this.finalPosition;
		return pathPositions;
	}
	
	public Trajectory getTrajectory(){
		return this.interpolatedPath;
	}
	
	public double getTimeSample(){
		return this.tSample;
	}
	
	public double getSpaceSample(){
		return this.xSample;
	}
	//END part of the GET methods
	
	
	//Part of the SET methods, where you can maybe set some parameters (not sure if needed)
	public void setMaxAcc(double acceleration){
		this.maxAcc = acceleration;
	}
	
	public void setMaxVel(double velocity){
		this.maxVel = velocity;
	}
	
	public void setXSample(double x_sample){
		this.xSample = x_sample;
	}
	
	public void setTSample(double t_sample){
		this.tSample = t_sample;
	}
	
	//END of the SET methods
	
	//OTHER PUBLIC methods
	
	public Trajectory interpolate(){
		/*
		 * This method interpolates the Path according to X_sample and a 
		 * Trapezoidal Velocity Profile (TVP). 
		 * (see literature at http://home.deib.polimi.it/rocco/cir/Motion%20planning.pdf)
		 * The number of points of the Trajectory strictly depends on the value of X_sample.
		 * This is the FATHER Interpolate method, that interpolates along a straight line
		 * from an initial position to a final position. 
		 */
		Trajectory interpolatedPath = new Trajectory();
		
		
		double accTime = this.maxVel/this.maxAcc;
		double[] initPos = this.initialPosition.getPosition();
		double[] finPos = this.finalPosition.getPosition();
		
		double distance = Math.sqrt(Math.pow(finPos[0] - initPos[0],2) + Math.pow(finPos[1] - initPos[1],2) + Math.pow(finPos[2] - initPos[2],2));
		
		
		//Log for debug
		//System.out.println("Distance : " + distance + " and acceleration time : " + acc_time);
		int N = 0;
		
		double tot_time = 0;
		double time = 0;
		double max_vel;
		double acc_distance;
		
		boolean acc_phase;
		boolean dec_phase;
		
		double new_position;
		
//		interpolatedPath.points.add(this.initialPosition);
//		interpolatedPath.timeInstants.add(time);
		
		if(((2*accTime)*this.maxVel/2.0) < distance)
		{
			tot_time = (distance - (2*accTime)*this.maxVel/2.0)/this.maxVel + 2*accTime;
			max_vel = this.maxVel;
			System.out.println("The maximum velocity is: " + this.maxVel + " m/s");
		}
		else
		{
			max_vel = Math.sqrt(distance*2*this.maxVel/(2*accTime));
			tot_time = distance*2/max_vel;
			accTime = accTime*max_vel/this.maxVel;
			System.out.println("The velocity profile is triangular, so the new maximum velocity is: " + max_vel + " m/s");
		}
		
		//log for debugging
		System.out.println("Distance: " + distance + " total time : " + tot_time + " and maximum velocity : " + max_vel);
		//System.out.println("Acceleration time: " + acc_time);
		
		acc_distance = max_vel*accTime/2;
		
		if(accTime > 0)
		{
			acc_phase = true;
		}
		else
		{
			acc_phase = false;
		}
		
		dec_phase = false;
		
		while(time <= tot_time)
		{
			++N;
			
			//
			Quat4d q = new Quat4d();
//			 q.interpolate(arg0, arg1, arg2);
			//
			
			
			if(time > accTime)
			{
				acc_phase = false;
			}
			if(time > (tot_time - accTime))
			{
				dec_phase = true;
			}
			
			if(acc_phase)
			{
				new_position = time*max_vel/accTime*time/2;
			}
			else
			{
				if(dec_phase)
				{
					new_position = distance - (tot_time - time)*(tot_time - time)*max_vel/(2*accTime);
				}
				else
				{
					new_position = acc_distance + (max_vel*(time - accTime));
				}
			}
			
			
			
			double[] position_vector = Matrix.subtractVectors(this.finalPosition.getPosition(), this.initialPosition.getPosition());
			position_vector = Matrix.multiplyScalMatr(new_position/distance, position_vector);
			position_vector = Matrix.addMatrices(position_vector, this.initialPosition.getPosition());
			
			//for the moment I'm not changing the rotation matrix
			Target newPointInTrajectory = new Target(position_vector, this.initialPosition.getRotation());
			
			//add the new target to the trajectory
			interpolatedPath.points.add(newPointInTrajectory);
			
//			time += this.tSample;
			
			//add the time instant to the trajectory's list
			interpolatedPath.timeInstants.add(time);
			time += this.tSample;
			
			
		}//end of the WHILE loop
		
		time = tot_time;
		
		interpolatedPath.points.add(this.finalPosition);
		interpolatedPath.timeInstants.add(time);
		//the time offset still has to be considered (maybe a method?)
		
		System.out.println("Sampling time : " + this.tSample);
		//System.out.println("Number of points N: " + N);
		
		this.interpolatedPath = interpolatedPath;
		return interpolatedPath;
	}
	
	public void translateTrajectory(double[] vector){
		for(Target i : this.interpolatedPath.points)
		{
			i.translateTarget(vector);
		}
	}
	
	public void translateTargets(double[] vector){
		this.initialPosition.translateTarget(vector);
		this.finalPosition.translateTarget(vector);
	}
	
	
	//END of PUBLIC methods
	
	
	//PRIVATE methods
	
	public static Quat4d eulerToQuaternion( float eulerX, float eulerY, float eulerZ ){
		/*
		 * X -> PSI
		 * y -> THETA
		 * Z -> PHI
		 */
		
		double c1 = Math.cos(eulerX/2.0);
		double s1 = Math.sin(eulerX/2.0);
		double c2 = Math.cos(eulerY/2.0);
		double s2 = Math.sin(eulerY/2.0);
		double c3 = Math.cos(eulerZ/2.0);
		double s3 = Math.sin(eulerZ/2.0);
		
		Quat4d q = new Quat4d((c3*c2*c1 + s3*s2*s1), (s3*c2*c1 - c3*s2*s1) , (c3*s2*c1 + s3*c2*s1) , (c3*c2*s1 - s3*s2*c1));
		
		
		return q;
	}
	
	public static List<double[]> rotm2eul (double[][] matrix){
		ArrayList<double[]> angles = new ArrayList<double[]>();
		
		//Copy the matrix elements to make formulas clear
		double r11 = matrix[0][0];
		double r12 = matrix[0][1];
		double r13 = matrix[0][2];
		double r21 = matrix[1][0];
		double r22 = matrix[1][1];
		double r23 = matrix[1][2];
		double r31 = matrix[2][0];
		double r32 = matrix[2][1];
		double r33 = matrix[2][2];
		//-----------------------------------
		
		if( Math.abs(r31) < 0.99999 || Math.abs(r31) > 1.00001){ //R31 is NOT +1 or -1
			double[] angles1 = new double[3];
			double[] angles2 = new double[3];
			
			angles1[1] = - Math.asin(r31); //THETA 1 
			angles2[1] = Math.PI - angles1[0]; // THETA 2
			
			angles1[0] = Math.atan2(r32/Math.cos(angles1[0]), r33/Math.cos(angles1[0])); //PSI 1
			angles2[0] = Math.atan2(r32/Math.cos(angles2[0]), r33/Math.cos(angles2[0])); //PSI 2
			
			angles1[2] = Math.atan2(r21/Math.cos(angles1[0]), r11/Math.cos(angles1[0])); //PHI 1
			angles2[2] = Math.atan2(r21/Math.cos(angles2[0]), r11/Math.cos(angles2[0])); //PHI 2
			
			angles.add(angles1);
			angles.add(angles2);
		} else{
			double[] angles1 = new double[3]; //In this case there will be only one solution
			
			angles1[2] = 0.0;

			if(0.99999 < r31 && r31 < 1.00001){ //R31 = 1
				angles1[1] = -Math.PI/2.0;
				
				angles1[0] = - angles1[2] + Math.atan2(-r12, -r13);
				
			}else{ //R31 = -1
				angles1[1] = Math.PI/2.0;
				
				angles1[0] = angles1[2] + Math.atan2(r12, r13);
			}
			
			angles.add(angles1);
		}
		
		
		
		
		return angles;
	}
}
