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
	protected Target initialPosition;
	protected Target finalPosition;
	protected Trajectory interpolatedPath; //the Trajectory will have also the initial and 
										  // final position within of course.
										  	
	
	protected double tSample; //time resolution for the trajectory
	protected double xSample; //space resolution for the trajectory
	
	protected double maxVel; //these two parameters MAYBE have to belong to the robot
	protected double maxAcc;
	
	protected boolean slerpOn;
	
	//Different constructors depending on how many parameters are specified
	
	public Path(){
		this.slerpOn = false;
	}
	
	public Path(Target initialPosition, Target finalPosition, double tSample, double xSample){
		this.initialPosition = initialPosition;
		this.finalPosition = finalPosition;
		this.tSample = tSample;
		this.xSample = xSample;
		
		if(Matrix.equals(this.initialPosition.getRotation(), this.finalPosition.getRotation())){
			this.slerpOn = false;
		} else{
			this.slerpOn = true;
		}
	}
	
	public Path(Target initialPosition, Target finalPosition)
	{
		this.initialPosition = initialPosition;
		this.finalPosition = finalPosition;
		
		this.tSample = 0.01;		//If a different resolution is required it has to be specified
		this.xSample = 1E-3;
		
		if(Matrix.equals(this.initialPosition.getRotation(), this.finalPosition.getRotation())){
			this.slerpOn = false;
		} else{
			this.slerpOn = true;
		}
	}
	//END part of constructor
	
	//Part of GET methods
	public Target getInitialPosition(){
		return this.initialPosition;
	}
	
	public Target getFinalPosition(){
		return this.finalPosition;
	}
	
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
	
	public double getMaxAcc(){
		return this.maxAcc;
	}
	
	public double getMaxVel(){
		return this.maxVel;
	}
	
	public boolean isSlerpOn() {
		return slerpOn;
	}

	//END part of the GET methods
	
	
	//Part of the SET methods, where you can maybe set some parameters (not sure if needed)
	public void setMaxAcc(double acceleration){
		this.maxAcc = acceleration;
	}
	
	public void setMaxVel(double velocity){
		this.maxVel = velocity;
	}
	
	public void setSlerpOn(boolean slerpOn) {
		this.slerpOn = slerpOn;
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
		
		double totTime = 0;
		double time = 0;
		double maxVel;
		double accDistance;
		
		boolean accPhase;
		boolean decPhase;
		
		double new_position;
		
		
 
		Quat4d q0 = this.rotmToQuat(this.initialPosition.getRotation());
		Quat4d q1 = this.rotmToQuat(this.finalPosition.getRotation());
		
//		interpolatedPath.points.add(this.initialPosition);
//		interpolatedPath.timeInstants.add(time);
		
		if(((2*accTime)*this.maxVel/2.0) < distance)
		{
			totTime = (distance - (2*accTime)*this.maxVel/2.0)/this.maxVel + 2*accTime;
			maxVel = this.maxVel;
//			System.out.println("The maximum velocity is: " + this.maxVel + " m/s");
		}
		else
		{
			maxVel = Math.sqrt(distance*2*this.maxVel/(2*accTime));
			totTime = distance*2/maxVel;
			accTime = accTime*maxVel/this.maxVel;
//			System.out.println("The velocity profile is triangular, so the new maximum velocity is: " + maxVel + " m/s");
		}
		
		//log for debugging
//		System.out.println("Distance: " + distance + " total time : " + totTime + " and maximum velocity : " + maxVel);
		//System.out.println("Acceleration time: " + acc_time);
		
		accDistance = maxVel*accTime/2;
		
		if(accTime > 0)
		{
			accPhase = true;
		}
		else
		{
			accPhase = false;
		}
		
		decPhase = false;
		
		while(time <= totTime)
		{
			++N;
			
			
			if(time > accTime)
			{
				accPhase = false;
			}
			if(time > (totTime - accTime))
			{
				decPhase = true;
			}
			
			if(accPhase)
			{
				new_position = time*maxVel/accTime*time/2;
			}
			else
			{
				if(decPhase)
				{
					new_position = distance - (totTime - time)*(totTime - time)*maxVel/(2*accTime);
				}
				else
				{
					new_position = accDistance + (maxVel*(time - accTime));
				}
			}
			
			
			
			double[] position_vector = Matrix.subtractVectors(this.finalPosition.getPosition(), this.initialPosition.getPosition());
			position_vector = Matrix.multiplyScalMatr(new_position/distance, position_vector);
			position_vector = Matrix.addMatrices(position_vector, this.initialPosition.getPosition());
			
			//
			double[][] rotm;
			
			if(slerpOn){
			
				Quat4d qm = this.slerp(q0, q1, time/totTime);
				rotm = this.quatToRotm(qm);			
			//this.quatToRotm(q0);
			}else{
				 rotm = this.initialPosition.getRotation();
			}
			//
			
			
			//The Targety is created.
			//Rotm is changed only if slerpOn was TRUE, that is if a rotation between initial position
			// and final position occurred.
			Target newPointInTrajectory = new Target(position_vector, rotm);
			
			//add the new target to the trajectory
			interpolatedPath.points.add(newPointInTrajectory);
			
//			time += this.tSample;
			
			//add the time instant to the trajectory's list
			interpolatedPath.timeInstants.add(time);
			time += this.tSample;
			
			
		}//end of the WHILE loop
		
		time = totTime;
		
		interpolatedPath.points.add(this.finalPosition);
		interpolatedPath.timeInstants.add(time);
		//the time offset still has to be considered (maybe a method?)
		
//		System.out.println("Sampling time : " + this.tSample);
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
	
	private Quat4d eulerToQuaternion( float eulerX, float eulerY, float eulerZ ){
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
	
	private List<double[]> rotm2eul (double[][] matrix){
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
	
	private Quat4d rotmToQuat (double[][] matrix){
		
		double m00 = matrix[0][0];
		double m01 = matrix[0][1];
		double m02 = matrix[0][2];
		double m10 = matrix[1][0];
		double m11 = matrix[1][1];
		double m12 = matrix[1][2];
		double m20 = matrix[2][0];
		double m21 = matrix[2][1];
		double m22 = matrix[2][2];
		
		double qx = 0;
		double qy = 0;
		double qz = 0;
		double qw = 0;
		
		double tr = m00 + m11 + m22;
		 
		if (tr > 0) { 
			double S = Math.sqrt(tr+1.0) * 2; // S=4*qw 
			qw = 0.25 * S;
			qx = (m21 - m12) / S;
			qy = (m02 - m20) / S; 
			qz = (m10 - m01) / S; 
		} else if ((m00 > m11)&(m00 > m22)) { 
			double S = Math.sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
			qw = (m21 - m12) / S;
			qx = 0.25 * S;
			qy = (m01 + m10) / S; 
			qz = (m02 + m20) / S; 
		} else if (m11 > m22) { 
			double S = Math.sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
			qw = (m02 - m20) / S;
			qx = (m01 + m10) / S; 
			qy = 0.25 * S;
			qz = (m12 + m21) / S; 
		} else { 
			double S = Math.sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
			qw = (m10 - m01) / S;
			qx = (m02 + m20) / S;
			qy = (m12 + m21) / S;
			qz = 0.25 * S;
		}
		 
		return new Quat4d(qx,qy,qz, qw);
		
	}
	
	private double[][] quatToRotm (Quat4d q){
		double[][] rotm = new double[3][3];
		
		    double sqw = q.w*q.w;
		    double sqx = q.x*q.x;
		    double sqy = q.y*q.y;
		    double sqz = q.z*q.z;

		    // invs (inverse square length) is only required if quaternion is not already normalised
		    double invs = 1 / (sqx + sqy + sqz + sqw);
		    rotm[0][0] = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
		    rotm[1][1] = (-sqx + sqy - sqz + sqw)*invs ;
		    rotm[2][2] = (-sqx - sqy + sqz + sqw)*invs ;
		    
		    double tmp1 = q.x*q.y;
		    double tmp2 = q.z*q.w;
		    rotm[1][0] = 2.0 * (tmp1 + tmp2)*invs ;
		    rotm[0][1] = 2.0 * (tmp1 - tmp2)*invs ;
		    
		    tmp1 = q.x*q.z;
		    tmp2 = q.y*q.w;
		    rotm[2][0] = 2.0 * (tmp1 - tmp2)*invs ;
		    rotm[0][2] = 2.0 * (tmp1 + tmp2)*invs ;
		    tmp1 = q.y*q.z;
		    tmp2 = q.x*q.w;
		    rotm[2][1] = 2.0 * (tmp1 + tmp2)*invs ;
		    rotm[1][2] = 2.0 * (tmp1 - tmp2)*invs ;      
		
		
		return rotm;
	}
	
	private Quat4d slerp(Quat4d qa, Quat4d qb, double t) {
	    // quaternion to return
		Quat4d qm = new Quat4d();
	    // Calculate angle between them.
	    double cosHalfTheta = qa.w * qb.w + qa.x * qb.x + qa.y * qb.y + qa.z * qb.z;
	    // if qa=qb or qa=-qb then theta = 0 and we can return qa
	    if (Math.abs(cosHalfTheta) >= 1.0){
	        qm.w = qa.w; 
	        qm.x = qa.x; 
	        qm.y = qa.y; 
	        qm.z = qa.z;
	        return qm;
	    }
	    // Calculate temporary values.
	    double halfTheta = Math.acos(cosHalfTheta);
	    double sinHalfTheta = Math.sqrt(1.0 - cosHalfTheta*cosHalfTheta);
	    // if theta = 180 degrees then result is not fully defined
	    // we could rotate around any axis normal to qa or qb
	    if (Math.abs(sinHalfTheta) < 0.001){ // fabs is floating point absolute
	        qm.w = (qa.w * 0.5 + qb.w * 0.5);
	        qm.x = (qa.x * 0.5 + qb.x * 0.5);
	        qm.y = (qa.y * 0.5 + qb.y * 0.5);
	        qm.z = (qa.z * 0.5 + qb.z * 0.5);
	        return qm;
	    }
	    double ratioA = Math.sin((1 - t) * halfTheta) / sinHalfTheta;
	    double ratioB = Math.sin(t * halfTheta) / sinHalfTheta; 
	    //calculate Quaternion.
	    qm.w = (qa.w * ratioA + qb.w * ratioB);
	    qm.x = (qa.x * ratioA + qb.x * ratioB);
	    qm.y = (qa.y * ratioA + qb.y * ratioB);
	    qm.z = (qa.z * ratioA + qb.z * ratioB);
	    return qm;
	}
	
}
