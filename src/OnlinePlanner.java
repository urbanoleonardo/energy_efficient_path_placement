/*
TO DO 
decide how to get the current position of the robot, so the starting point that is
needed for the inverse kinematics. 
*/
public class OnlinePlanner implements Runnable{

	private Target currPosition;
	private int targetsLength;
	private Target[] targets;
	private Robot robot;
	
	private double[][] link_length;
	private double[][] joint_limits;
	
	
	public OnlinePlanner(Target[] targets, Robot robot)
	{
		
		targetsLength = copytargets(targets);
		this.robot = robot;
		this.link_length = robot.getParameters("link_length");
		this.joint_limits = robot.getParameters("joint_limits");
	}
	
	public void run(){
		if (targetsLength == 1)
	{
	
	}
	else
	{
		
	}
		
		
	}
	
	
	private int copytargets(Target[] targetsTOstore)
	{
		int length = targetsTOstore.length;
		this.targets = new Target[length];
		for(int i = 0; i<length; i++)
		{
			this.targets[i] = targetsTOstore[i];
		}
		return length;
	}
	
	
	private void Online_Kinematics(Target Current_position, Target Next_position)
	{
		/*
		 * Need to get the coordinates of the WIRST center Pw
		 */
	}
	
	
	private void Online_Dynamics()
	{
		
	}
	
	private double[] solve_theta1(double[] Pw)
	{
		/*
		 * Input Pw: vector 3X1 of wirst coordinates
		 */
		
		//Instantiate the vector
		double[] theta_1 = new double[2];
		
		//calculate first theta1
		theta_1[0] = Math.atan2(Pw[1], Pw[0]);
		
		//calculate second theta1
		theta_1[1] = Math.atan2(-Pw[1], -Pw[0]);
		
		return theta_1;
	}
	
	/*
	 * Our notation is different form the paper.
	 * We adopt the iterative DH notation.
	 * If it's confronted with MATLAB module we have:
	 * MATLAB  ------   JAVA
	 * a1				a2
	 * a2				a4
	 * a3				a5
	 * a4				a6
	 * -------------------
	 * d1				d1
	 * d2				d3
	 * 
	 */
	
	private double[][] solve_theta2_3(double[] Pw, double theta_1 , int index)
	{
		
		double a2 = link_length[1][0];
		double a4 = link_length[3][0];
		double a5 = link_length[4][0];
		
		double d1 = link_length[0][1];
		double d3 = link_length[2][1];
		
		//standard values for initialization if target is out of reach
		double[] theta_3 = {1000 , 0}; 
		double[] theta_2 = {1000, 0};
		
		double[][] theta2_3 = new double[2][2];
	
		double rr;
		double a;
		
		if(index == 1)
		{
			rr = Math.sqrt(Math.pow((Pw[0] - a2*Math.cos(theta_1)), 2) + Math.pow((Pw[2] - a2*Math.sin(theta_1)), 2));
		}
		else
		{
			rr = - Math.sqrt(Math.pow((Pw[0] - a2*Math.cos(theta_1)), 2) + Math.pow((Pw[2] - a2*Math.sin(theta_1)), 2));
		}
		//a is the cosine of THETA 3
		a = ((-(d3*d3) - (a4 + a5)*(a4 + a5) + (Pw[2] - d1)*(Pw[2] - d1) + rr*rr ))/(2*d3*(a4 + a5));
		if(a < -1 || a > 1 )
		{
			theta2_3[0] = theta_2;
			theta2_3[1] = theta_3;
			System.out.println("Impossible to reach that target.");
			return theta2_3;
		}
		else
		{
			theta_3[0] = Math.atan2(Math.sqrt(1 - a4*a4), a);
			theta_3[1] = Math.atan2(-Math.sqrt(1 - a4*a4), a);
			
			theta2_3[1] = theta_3;
		}
		for(int i = 0; i < theta_3.length; i++)
		{
			theta_2[i] = Math.atan2((Pw[2] - d1),rr) + Math.atan2((a4 + a5)*Math.sin(theta_3[i]), d3 + (a4 + a5)*Math.cos(theta_3[i])); 
		}
		
		theta2_3[0] = theta_2;
		
		return theta2_3;
		
	}
	
	private double[][] solve_theta4_5_6(double[] Pw, double theta_1, double theta_2, double theta_3){
		
		/*
		 * TODO
		 * 1: get the target rotation
		 * 2: get from the robot the rotational matrix from 0 to TOOL
		 * 3: build the matrix from which the Euler's angles will be calculated
		 * 
		 * IMPORTANT: needed a way to read the joint angles of the previous configuration 
		 */
			
		double[][] R_zyz = {
				{0.45,  0.0,  0.20},
				{0.20, -0.90, 0.10},
				{0.30, -0.10, 1.00}
		};
		
		double previous_theta4_TEST = 0; //has to be copied from the previous configuration
		
		double theta_46;
		double[] theta_4 = new double[2];
		double[] theta_5 = new double[2];
		double[] theta_6 = new double[2];
		
		//double[][] theta_456;
		
 		double temp = Math.sqrt(R_zyz[2][0]*R_zyz[2][0] + R_zyz[2][1]*R_zyz[2][1]);
		
		theta_5[0] = Math.atan2(temp, R_zyz[2][2]);
		theta_5[1] = Math.atan2(-temp, R_zyz[2][2]);
		
		for(int i = 0; i < theta_5.length ; i++)
		{
			if(-0.0001 < theta_5[i] && theta_5[i] < 0.0001)
			{//theta_5 is 0 and we're in a singularity case
				theta_4[i] = previous_theta4_TEST;
				theta_46 = Math.atan2(R_zyz[1][0], R_zyz[0][0]);
				
				theta_6[i] = theta_46 - theta_4[i];
				
			}
			else
			{//theta_5 is NOT 0 and we can calculate the other angles
				temp = Math.sin(theta_5[i]);
				theta_4[i] = Math.atan2(R_zyz[0][2], -R_zyz[1][2])/temp;
				theta_6[i] = Math.atan2(R_zyz[2][0], -R_zyz[2][1])/temp;
			}
		}
		
		double[][] theta_456 = new double[3][2];
		theta_456[0] = theta_4;
		theta_456[1] = theta_5;
		theta_456[2] = theta_6;
		
		return theta_456;
		/*
		 * REMINDER
		 * in case performance change significantly it could be possible to assign
		 * the theta values directly to the 3x2 matrix and then return it.
		 */
	}
	
}
