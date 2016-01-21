
public class Path {
/*
 * Path class that contains initial and final TARGETS, and then has a Trajectory attribute
 * that will be filled if the INTERPOLATE method will be called (still have to see if it will 
 * automatically fill the attribute or whatnot).
 * The "Father" method INTERPOLATE will generate basically a straight line
 * while for other paths like the square one (that will be sons) the method will behave differently 
 */
	private Target initial_position;
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
		this.initial_position = initial_position;
		this.final_position = final_position;
		this.t_sample = t_sample;
		this.x_sample = x_sample;
	}
	
	public Path(Target initial_position, Target final_position)
	{
		this.initial_position = initial_position;
		this.final_position = final_position;
	}
	//END part of constructor
	
	//Part of GET methods
	public Target[] GetPathPositions()
	{
		
		/*
		 * It returns a vector of Target with initial and final position
		 */
		
		Target[] PathPositions = new Target[2];
		PathPositions[0] = this.initial_position;
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
	
	public void SetXSample(double x_sample){
		this.x_sample = x_sample;
	}
	
	public void SetTSample(double t_sample){
		this.t_sample = t_sample;
	}
	
	//END of the SET methods
	
	//OTHER PUBLIC methods
	
	public Trajectory Interpolate(){
		/*
		 * This method interpolates the Path according to X_sample and a 
		 * Trapezoidal Velocity Profile (TVP). 
		 * (see literature at http://home.deib.polimi.it/rocco/cir/Motion%20planning.pdf)
		 * The number of points of the Trajectory strictly depends on the value of X_sample.
		 * This is the FATHER Interpolate method, that interpolates along a straight line
		 * from an initial position to a final position. 
		 */
		Trajectory interpolatedPath = new Trajectory();
		
		
		double acc_time = this.max_acc/this.max_vel;
		double[] init_pos = this.initial_position.getPosition();
		double[] fin_pos = this.final_position.getPosition();
		
		double distance = Math.sqrt(Math.pow(fin_pos[0] - init_pos[0],2) + Math.pow(fin_pos[1] - init_pos[1],2) + Math.pow(fin_pos[2] - init_pos[2],2));
		
		//Log for debug
		//System.out.println("Distance : " + distance + " and acceleration time : " + acc_time);
		
		double tot_time = 0;
		double time = 0;
		double max_vel;
		double acc_distance;
		
		boolean acc_phase;
		boolean dec_phase;
		
		double new_position;
		
		interpolatedPath.Points.add(this.initial_position);
		interpolatedPath.TimeInstants.add(time);
		
		if(((2*acc_time)*this.max_vel/2.0) < distance)
		{
			tot_time = (distance - (2*acc_time)*this.max_vel/2.0)/this.max_vel + 2*acc_time;
			max_vel = this.max_vel;
			System.out.println("The maximum velocity is: " + this.max_vel + " m/s");
		}
		else
		{
			max_vel = Math.sqrt(distance*2*this.max_vel/(2*acc_time));
			tot_time = distance*2/max_vel;
			acc_time = acc_time*max_vel/this.max_vel;
			System.out.println("The velocity profile is triangular, so the new maximum velocity is: " + max_vel + " m/s");
		}
		
		acc_distance = max_vel*acc_time/2;
		
		if(acc_time > 0)
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
			if(time > acc_time)
			{
				acc_phase = false;
			}
			if(time > (tot_time - acc_time))
			{
				dec_phase = true;
			}
			
			if(acc_phase)
			{
				new_position = time*max_vel/acc_time*time/2;
			}
			else
			{
				if(dec_phase)
				{
					new_position = distance - (tot_time - time)*(tot_time - time)*max_vel/(2*acc_time);
				}
				else
				{
					new_position = acc_distance + (max_vel*(time - acc_time));
				}
			}
			
			double[] position_vector = vector_subtract(this.final_position.getPosition(), this.initial_position.getPosition());
			position_vector = vector_MultiplyConstant(position_vector, new_position/distance);
			position_vector = vector_sum(position_vector, this.initial_position.getPosition());
			
			//for the moment I'm not changing the rotation matrix
			Target newPointInTrajectory = new Target(position_vector, this.initial_position.getRotation());
			
			//add the new target to the trajectory
			interpolatedPath.Points.add(newPointInTrajectory);
			
			time += this.t_sample;
		}//end of the WHILE loop
		
		interpolatedPath.Points.add(this.final_position);
		interpolatedPath.TimeInstants.add(time);
		//the time offset still has to be considered (maybe a method?)
		
		this.interpolated_path = interpolatedPath;
		return interpolatedPath;
	}
	
	public void TranslateTrajectory(double[] vector)
	{
		for(Target i : this.interpolated_path.Points)
		{
			i.TranslateTarget(vector);
		}
	}
	
	
	//END of PUBLIC methods
	
	
	//PRIVATE methods
	
	private double[] vector_sum(double[] a, double[] b)
	{
		double[] result = new double[a.length];
		for(int i = 0; i < a.length; i++)
		{
			result[i] = a[i] + b[i];
		}
		return result;
	}
	
	private double[] vector_subtract(double[] a, double[] b)
	{
		double[] result = new double[a.length];
		for(int i = 0; i < a.length; i++)
		{
			result[i] = a[i] - b[i];
		}
		return result;
	}
	
	private double[] vector_MultiplyConstant(double[] a, double constant)
	{
		
		for(int i = 0; i < a.length; i++)
		{
			a[i] = a[i] * constant;
		}
		return a;
	}
	
	

	
}
