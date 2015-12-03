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
		 */
		double[][] dummy = new double[3][3];
		
		double[][] R_zyz = {
				{0.45,  0.0,  0.20},
				{0.20, -0.90, 0.10},
				{0.30, -0.10, 1.00}
		};
		
		double[] theta_5 = new double[2];
		double temp = Math.sqrt(R_zyz[2][0]*R_zyz[2][0] + R_zyz[2][1]*R_zyz[2][1]);
		
		theta_5[0] = Math.atan2(temp, R_zyz[2][2]);
		theta_5[1] = Math.atan2(-temp, R_zyz[2][2]);
		
		for(int i = 0; i < theta_5.length ; i++)
		{
			
			
			
			
			
		}
		
		
		return dummy;
	}
	/*
	 * %euler angles Z-X-Z
theta5(1)=atan2(sqrt(Rzyz(3,1)^2+Rzyz(3,2)^2),Rzyz(3,3));
theta5(2)=atan2(-sqrt(Rzyz(3,1)^2+Rzyz(3,2)^2),Rzyz(3,3));

for cont10=1:2
    % in case of being a singularity
    if (theta5(cont10)<0.0001 && theta5(cont10)>-0.0001)
        if t==1
            theta_4f=0;
        else
            if theta5(cont10)<0
                for i=1:er(t-1)
                    if thetaM(5,i,t-1)<0
                        theta_sing=[theta_sing,abs(thetaM(1:3,i,t-1)-[theta1;theta2;theta3])];
                    end
                end
                a=size(theta_sing,2);
                [A,B]=min(sum(theta_sing(:,2:a),1));
            else %theta5 could be < or = 0
                if theta5(cont10)== 0
                    
                    for i=1:er(t-1)
                        
                        theta_sing=[theta_sing,abs(thetaM(1:3,i,t-1)-[theta1;theta2;theta3])];
                        
                    end
                    
                    a=size(theta_sing,2);
                    [A,B]=min(sum(theta_sing(:,2:a),1));
                    
                else %theta5 is <0
                    for i=1:er(t-1)
                        if thetaM(5,i,t-1)>0
                            theta_sing=[theta_sing,abs(thetaM(1:3,i,t-1)-[theta1;theta2;theta3])];
                        end
                    end
                    a=size(theta_sing,2);
                    [A,B]=min(sum(theta_sing(:,2:a),1)); 
                    %I am looking basically for the minimum theta4 needed to reach the desired orientation starting from 
                    %the previous one (t-1). In this way theta_4f is
                    %calculated. The overall rotation needed from the
                    %joints 4 and 6 is calculated as theta_46. Then if
                    %theta_4f doesn't exists it's assumed to be zero, and
                    %the sixth joint will do all the rotation, otherwise we
                    %rotate joint-4 of theta_4f and the joint-6 of the
                    %remaining angle (theta_46 - theta_4f)
                    
                    
                end
                
            end
            
            theta_4f=thetaM(4,B,t-1);
        end
        
        theta_46(1)=atan2(Rzyz(2,1),Rzyz(1,1));
        
        try
        theta6(cont10)= theta_46- theta_4f;
        theta4(cont10)=theta_4f;%assumption: the previous theta4
        catch %% if there is no previous theta4f, error rise, then 0 is assumed in catch
            theta6(cont10)= theta_46;
            theta4(cont10)=0;%assumption: the previous theta4
        end
        
    else %IS NOT a singularity
        %then we calculate the other angles with the Euler equations
        theta4(cont10)=atan2(Rzyz(1,3)/sin(theta5(cont10)),-Rzyz(2,3)/sin(theta5(cont10)));
        theta6(cont10)=atan2(Rzyz(3,1)/sin(theta5(cont10)),Rzyz(3,2)/sin(theta5(cont10)));
    end
end
	 */
}
