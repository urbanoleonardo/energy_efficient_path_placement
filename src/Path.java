
public class Path {
/*
 * Path class that contains initial and final TARGETS, and then has a Trajectory attribute
 * that will be filled if the INTERPOLATE method will be called (still have to see if it will 
 * automatically fill the attribute or whatnot).
 * The "Father" method INTERPOLATE will generate basically a straight line
 * while for other paths like the square one (that will be sons) the method will behave differently 
 */
	private Target intial_position;
	private Target final_position;
	private Trajectory interpolated_path; //the Trajectory will have also the initial and 
										  // final position within of course.
										  	
	
	private double t_sample; //time resolution for the trajectory
	private double x_sample; //space resolution for the trajectory
	
	private double max_vel; //these two parameters MAYBE have to belong to the robot
	private double max_acc;
	
	//Different constructors depending on how many parameters are specified
	
	public Path(Target initial_position, Target final_position, double t_sample, double x_sample)
	{
		this.intial_position = initial_position;
		this.final_position = final_position;
		this.t_sample = t_sample;
		this.x_sample = x_sample;
	}
	
	public Path(Target initial_position, Target final_position)
	{
		this.intial_position = initial_position;
		this.final_position = final_position;
	}
	//END part of constructor
	
	//Part of GET methods
	public Target[] GetPathPositions()
	{
		Target[] PathPositions = new Target[2];
		PathPositions[0] = this.intial_position;
		PathPositions[1] = this.final_position;
		return PathPositions;
	}
	
	public Trajectory GetTrajectory(){
		return this.interpolated_path;
	}
	
	public double GetTimeSample(){
		return this.t_sample;
	}
	
	public double GetSpaceSample(){
		return this.x_sample;
	}
	//END part of the GET methods
	
	
	//Part of the SET methods, where you can maybe set some parameters (not sure if needed)
	public void SetMaxAcc(double acceleration){
		this.max_acc = acceleration;
	}
	
	public void SetMaxVel(double velocity){
		this.max_vel = velocity;
	}
	
	//END of the SET methods
	
	//OTHER PUBLIC methods
	
	public Trajectory Interpolate(){
		/*
		 * This method interpolates the Path according to X_sample and a 
		 * Trapezoidal Velocity Profile (TVP). The number of points of the Trajectory
		 * strictly depends on the value of X_sample.
		 * This is the FATHER Interpolate method, that interpolates along a straight line
		 * from an initial position to a final position. 
		 */
		Trajectory interpolatedPath = new Trajectory();
		
		double acc_time = this.max_acc/this.max_vel;
		double[][] init_pos = this.intial_position.getPosition();
		double[][] fin_pos = this.final_position.getPosition();
		
		//double distance = Math.sqrt(Math.pow(fin_pos[0] - init_pos[0],2) + Math.pow(fin_pos[1] - init_pos[1],2) + Math.pow(fin_pos[2] - init_pos[2],2));
		interpolatedPath.Points.add(this.intial_position);
		/*
		if(((2*acc_time)*this.max_vel/2.0) < distance)
		{
			
		}
		*/
		
		return interpolatedPath;
	}
	
	
	
	//END of PUBLIC methods
	
	
	//PRIVATE methods
	
	
	
	

	
}
