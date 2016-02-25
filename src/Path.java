import com.sun.javafx.geom.Quat4f;

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
		
		Target[] PathPositions = new Target[2];
		PathPositions[0] = this.initialPosition;
		PathPositions[1] = this.finalPosition;
		return PathPositions;
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
	
	public void TranslateTrajectory(double[] vector)
	{
		for(Target i : this.interpolatedPath.points)
		{
			i.translateTarget(vector);
		}
	}
	
	
	//END of PUBLIC methods
	
	
	//PRIVATE methods
	
	public static Quat4f eulerToQuaternion( float eulerX, float eulerY, float eulerZ ){
		Quat4f q = new Quat4f();
		
		
		
		return q;
	}
}
