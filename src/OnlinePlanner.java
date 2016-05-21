import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

//import org.apache.poi.ss.usermodel.Row;
//import org.apache.poi.ss.usermodel.Sheet;
//import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;


/*
TODO 
decide how to get the current position of the robot, so the starting point that is
needed for the inverse kinematics. 
 */
public class OnlinePlanner implements Runnable{

	private static final String FILE_PATH = "src/Excel Output Files/testTrajectoryData.xlsx";

	private Target currPosition;
	private int targetsLength;
	private Target[] targets;
	private Path path = null;
	private Robot robot;

	private double[][] linkLength;
	private double[][] jointLimits;
	private double[] speedLimits;
	private double[] mass;
	private List<double[][]> inertiaTens;

	private double[][][] thetaM;

	private double acceleration;
	private double velocity;

	private double[][] rici; //it has to be structured as 6x3 matrix to work down in the dynamic analysis
	private double[][] rc; //6x3 matrix as well. The reason is because in loops 1x3 vectors are needed and will be taken by having r_c[i] instructions
	private double[][] r;

	private List<EnergyPoint> energyList;
	private Point3D point3D = null;
	private List<Point3D> points3D = null;


	public OnlinePlanner(Robot robot)
	{

		this.robot = robot;
		this.copyRobotParameters();	

		this.energyList = new LinkedList<EnergyPoint>();
	}

	public OnlinePlanner(Target[] targets, Robot robot)
	{

		targetsLength = copyTargets(targets);
		this.robot = robot;
		this.copyRobotParameters();		

		this.currPosition = targets[0];

		this.energyList = new LinkedList<EnergyPoint>();
	}

	public OnlinePlanner(Target target1, Target target2, Robot robot)
	{

		this.robot = robot;
		this.copyRobotParameters();		

		this.path = new Path(target1, target2);
		path.setMaxAcc(acceleration);
		path.setMaxVel(velocity);

		this.energyList = new LinkedList<EnergyPoint>();
	}

	public OnlinePlanner(Target target1, Target target2, Robot robot, Point3D point3D)
	{

		this.robot = robot;
		this.copyRobotParameters();		

		this.path = new Path(target1, target2);
		path.setMaxAcc(acceleration);
		path.setMaxVel(velocity);

		this.energyList = new LinkedList<EnergyPoint>();
		this.point3D = point3D;
	}

	public OnlinePlanner (Path p, Robot robot){
		
		this.robot = robot;
		this.copyRobotParameters();		
		
		this.path = p;
		this.energyList = new LinkedList<EnergyPoint>();
	}
	
	public OnlinePlanner(Path p, Robot robot, List<Point3D> points3D)
	{

		this.robot = robot;
		this.copyRobotParameters();		

		this.path = p;
		path.setMaxAcc(acceleration);
		path.setMaxVel(velocity);

		this.energyList = new LinkedList<EnergyPoint>();
		this.points3D = points3D;
	}

	public void run(){
		if(this.points3D == null){
			if(this.path == null){
				double xSample = 1E-3;
				double tSample = 0.01;
				boolean canSolve = true;

				for( int i = 0; i < targetsLength && canSolve; i++){
					//Line to be removed
					if(currPosition == targets[i]){
						//					System.out.println("Current position is equal to the next target.");
						continue;
					}

					this.path = new Path(currPosition, targets[i], tSample, xSample);
					path.setMaxAcc(acceleration);
					path.setMaxVel(velocity);

					canSolve = solver(path);
					currPosition = targets[i];

					if(point3D != null){
						//I need to add the position vector to the energyPoint just calculated
						EnergyPoint newPoint = this.energyList.get(this.energyList.size() - 1); 
						newPoint.setPosition(point3D.getPosition());
						this.energyList.set(this.energyList.size() - 1, newPoint);
					}
				}


			}else {

				solver(path);

				if(point3D != null){
					//I need to add the position vector to the energyPoint just calculated
					EnergyPoint newPoint = this.energyList.get(this.energyList.size() - 1); 
					newPoint.setPosition(point3D.getPosition());
					this.energyList.set(this.energyList.size() - 1, newPoint);
				}
			}
		} else {

			Target[] targets = path.getPathPositions();

			double[][] inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
			double[][] finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

			for(Point3D point : points3D){

				Target[] curr = new Target[2];
				curr[0] = new Target(inPosH);
				curr[1] = new Target(finPosH);
				double[] currPos = Matrix.copyVector(point.getPosition());
				curr[0].translateTarget(currPos);
				curr[1].translateTarget(currPos);


				Path currentPath = new Path(curr[0], curr[1]);
				currentPath.setMaxAcc(acceleration);
				currentPath.setMaxVel(velocity);

				solver(currentPath);

				EnergyPoint newPoint = this.energyList.get(this.energyList.size() - 1); 
				newPoint.setPosition(point.getPosition());
				this.energyList.set(this.energyList.size() - 1, newPoint);
			}

		}


		this.path = null;
		this.point3D = null;
		this.points3D = null;
	}

	public void run(Target newTarget){
		targets = new Target[1];
		targetsLength = 1;
		targets[0] = newTarget;
		this.run();

	}

	public void run(Target target1,Target target2 ){
		double xSample = 1E-3;
		double tSample = 0.01;

		targets = new Target[2];
		targetsLength = 2;
		currPosition = target1;
		targets[0] = target1;
		targets[1] = target2;

		this.path = new Path(target1, target2, tSample, xSample);
		path.setMaxAcc(acceleration);
		path.setMaxVel(velocity);

		solver(path);
		this.path = null;

		//		this.run();
	}

	public void run(Target target1,Target target2, Point3D point3D){
		double xSample = 1E-3;
		double tSample = 0.01;

		targets = new Target[2];
		targetsLength = 2;
		currPosition = target1;
		targets[0] = target1;
		targets[1] = target2;

		this.path = new Path(target1, target2, tSample, xSample);
		path.setMaxAcc(acceleration);
		path.setMaxVel(velocity);

		solver(path);

		//I need to add the position vector to the energyPoint just calculated
		EnergyPoint newPoint = this.energyList.get(this.energyList.size() - 1); 
		newPoint.setPosition(point3D.getPosition());
		this.energyList.set(this.energyList.size() - 1, newPoint);

		this.path = null;

		//		this.run();
	}

	public void run(Path path){
		this.path = path;
		this.run();
	}

	public void run(Path path, Point3D point3D){
		this.path = path;

		this.solver(path);

		//I need to add the position vector to the energyPoint just calculated
		EnergyPoint newPoint = this.energyList.get(this.energyList.size() - 1); 
		newPoint.setPosition(point3D.getPosition());
		this.energyList.set(this.energyList.size() - 1, newPoint);

		this.path = null;
	}

	private int copyTargets(Target[] targetsToStore)
	{
		int length = targetsToStore.length;
		this.targets = new Target[length];
		for(int i = 0; i<length; i++)
		{
			this.targets[i] = targetsToStore[i];
		}
		return length;
	}



	/**
	 * solver substitutes run() that previously was doing all the work.
	 * In solver the following actions are performed:
	 * 	- interpolation of the path given in input								//SECTION 1
	 * 	- inverse kinematics of the trajectory (obtained from path)				//SECTION 2
	 * 	- limitation are applied to the resulted configuration					//SECTION 3
	 * 	- configuration are also checked if they respect speed limitations		//SECTION 4
	 * 	- inverse dynamics is performed on the remaining viable configurations.	//SECTION 5
	 * 	- total energy is calculated for each solution							//SECTION 6
	 * 	- the minimum energy cost is saved and added to energyList.				//SECTION 7
	 * 	- the selected solution is temporary saved if is needed.				//SECTION 8
	 *
	 *	@param inputPath: the input Path to interpolate and calculate inv kinematics and dynamics
	 *
	 *	@return A boolean that is TRUE if a solution has been found, FALSE otherwise.
	 */


	private boolean solver(Path inputPath){


		//		|-----------|
		//		| SECTION 1 |
		//		|-----------|
		Trajectory pathTrajectory;
		int trajectoryLength;
		if(inputPath.interpolatedPath == null){
			pathTrajectory = inputPath.interpolate();
			trajectoryLength = pathTrajectory.points.size();
		} else{
			pathTrajectory = inputPath.getTrajectory();
			trajectoryLength = pathTrajectory.points.size();
		}

		//		System.out.println("Interpolation over. Points size: " + trajectoryLength + " Time instants : " + pathTrajectory.timeInstants.size());
		//		System.out.println("Number of points in Ext Torques " + pathTrajectory.extTorques.size());
		//		for(int l = 0 ; l < pathTrajectory.timeInstants.size(); l++){
		//			System.out.println(pathTrajectory.timeInstants.get(l));
		//		}

		//		for(int l = 0 ; l < pathTrajectory.points.size(); l++){
		//			Matrix.displayVector(pathTrajectory.points.get(l).getPosition());
		//		}

		thetaM = new double[6][8][trajectoryLength];

		long time1 = System.nanoTime();
		int error = 0;

		//		|-------------|
		//		| SECTION 2&3 |
		//		|-------------|

		for(int j = 0; j < trajectoryLength && error == 0; j++){

			error = onlineInvKinematics(pathTrajectory.points.get(j), j);
		}

		long time2 = System.nanoTime();

		//		|-----------|
		//		| SECTION 4 |
		//		|-----------|

		long timeConfigStart = System.nanoTime();

		int[] Solution_vector = new int[thetaM[0][0].length];
		ArrayList<int[]> Solutions = new ArrayList<int[]>();
		configIteration(thetaM, 0, 0, pathTrajectory.timeInstants,Solution_vector, Solutions);

		long timeConfigEnd = System.nanoTime();

		//		displayThetaM(134);


		int NumOfSolutions = Solutions.size();
		//		System.out.println("End of Path Planning. The number of solutions found is : " + NumOfSolutions);
		//		System.out.println("Time for the inverse Kinematics: " + (time2-time1)/1E9);

		//		//FOR DEBUG
		//		System.out.println("");
		//		for(int j = 0 ; j < Solutions.get(3).length; j++){
		//			System.out.println(Solutions.get(1)[j]);
		//		}
		//

		if(NumOfSolutions == 0)
		{
			/*
			 * If there are no feasible solutions we print an error log and we add ZERO energy to the list.
			 */
			//			System.out.println("Impossible to reach that target since there is no solution to the inverse kinematic problem.");
			EnergyPoint newPoint = new EnergyPoint(0d);
			this.energyList.add(newPoint);
			return false;
		}


		//		|-----------------|
		//		| SECTION 5,6,7,8 |
		//		|-----------------|
		long time3 = System.nanoTime();
		onlineInvDynamics(pathTrajectory, Solutions);
		long time4 = System.nanoTime();


		System.out.println("Time for the inverse Kinematics: " + (time2-time1)/1E9);
		System.out.println("Time for Config iteration: " + (timeConfigEnd-timeConfigStart)/1E9);
		System.out.println("Time for the inverse Dynamics: " + (time4-time3)/1E9); //MOST time consuming task
		System.out.println("Time for the whole procedure: " + (time4-time1)/1E9);

		return true;
	}
	
	
	
	/**
	 * Method that perform the inverse kinematics for the input position: the joint
	 * values are calculated to make the robot reach the specified target. 
	 * 
	 * @param nextPosition the position in which calculate the inverse kinematics
	 * @param point_i index to know in which position of ThetaM the result will be stored
	 * @return an integer that is -1 if an error occurred (the point is out of reach) or 0 otherwise.
	 */

	private int onlineInvKinematics(Target nextPosition, int point_i)
	{
		/*
		 * Need to get the coordinates of the WIRST center Pw
		 */

		//		double[] P = nextPosition.getPosition();
		double a6 = linkLength[5][0];
		double[] gw = {a6*nextPosition.getRotation()[0][2], a6*nextPosition.getRotation()[1][2], a6*nextPosition.getRotation()[2][2]};
		double[] Pw = Matrix.subtractVectors(nextPosition.getPosition(),gw);

		for(int i = 0; i < Pw.length; i++){Pw[i] = nextPosition.getPosition()[i] - gw[i];}


		double[] theta1 = solveTheta1(Pw);

		for(int i = 0; i < theta1.length; i++)
		{
			double[][] theta23 = solveTheta23(Pw, theta1[i], i);

			if(theta23[0][0] == 1000 || theta23[1][0] == 1000)
			{
				//				System.out.println("Trajectory impossible to track: theta_2 or theta_3 assume impossible values");
				//				System.out.println("hence, the whole sequence of points is impossible.");
				return -1; //will be used by the caller to stop the whole trajectory
			}

			for(int j = 0; j < theta23.length; j++)
			{
				double theta2 = theta23[0][j] - Math.PI/2;
				double theta3 = Math.PI/2 - theta23[1][j];
				if(theta3 < -Math.PI/2)
				{
					theta3 += 2*Math.PI; 
				}

				double[][] theta_456 = solveTheta456(nextPosition, theta1[i], theta2, theta3, point_i ); 



				for(int k = 0; k < theta_456[0].length; k++)
				{
					int index  = (int) ((i)*Math.pow(2, 2*(i)) + (j)*2 + k + 1);
					index--;


					thetaM[0][index][point_i] = theta1[i];
					thetaM[1][index][point_i] =	theta2;
					thetaM[2][index][point_i] =	theta3;
					thetaM[3][index][point_i] = theta_456[0][k];
					thetaM[4][index][point_i] =	-theta_456[1][k];
					thetaM[5][index][point_i] = theta_456[2][k];

				}

			}
		}

		//CHECK limitation
		thetaM = this.applyLimitation(thetaM, point_i);

		//		if(point_i == 0){displayThetaM( point_i );}
		//		if(point_i == 1){displayThetaM( point_i );}
		//		if(point_i == 2){displayThetaM( point_i );}
		//		if(point_i == (thetaM[0][0].length - 2)){displayThetaM( point_i );}
		//		if(point_i == (thetaM[0][0].length - 1)){displayThetaM( point_i );}
		//So far thetaM seems correct

		return 0;
	}


	private void onlineInvDynamics(Trajectory trajectory, ArrayList<int[]> kinSolutions){

		/*
		 * Inverse dynamics and energy calculation is performed within this method.
		 * The minimum energy among the different solutions is stored into energyList.
		 */

		int N = trajectory.points.size();
		double minEnergy = Double.MAX_VALUE;
		int[] optKinSolution = kinSolutions.get(0);
		double[][] optDynSolution = new double[N][6];

		long time1 = 0;
		long time2 = 0;

		for(int[] solution : kinSolutions)
		{
			int firstConf = solution[0];
			double[] prevTheta = new double[6]; 
			double[] prevDTheta = {0, 0 , 0, 0 , 0, 0};
			double energy = 0;
			double[][] dynSolutionVector = new double[N][6];

			for(int i = 0; i < prevTheta.length; i++){prevTheta[i] = thetaM[i][firstConf][0];} //initialization of the vector PREV_THETA

			//			System.out.println("Theta vector (prev_theta) : ");
			//			Matrix.displayVector(prev_theta);


			for(int i = 1 ; i < N; i++) //i should start from 1----- i from 1 to 134
			{
				//MAIN LOOP where the dynamic analysis is performed
				int index = solution[i]; //there was a +1 - WHY?
				double[] theta = new double[6];
				for(int j = 0; j < 6; j++){ theta[j] = thetaM[j][index][i];} //there was a +1 - WHY?
				double[] dtheta = new double[6];
				double[] ddtheta = new double[6];

				//				if(i == 71){
				//					System.out.println("Theta vector for i=" + i + ": ");
				//					Matrix.displayVector(theta);
				//					System.out.println("Prev_Theta vector for i=" + i + ": ");
				//					Matrix.displayVector(prev_theta);
				//				}

				double dTime = (trajectory.timeInstants.get(i) - trajectory.timeInstants.get(i-1));
				double invDTime = 1/dTime;

				for(int j = 0; j < 6; j++){ 
					dtheta[j] = (theta[j] - prevTheta[j])*invDTime;
					ddtheta[j] = (dtheta[j] - prevDTheta[j])*invDTime;
				}



				//	FROM MATLAB   T_ex=Ti(:,cont7);
				//                f_ex=fi(:,cont7);

				//The following two lines will have to be substitued with the proper external forces
				//and torques when the interpolation has been changed accordingly

				double[] T_ext = {0,0,0};
				double[] f_ext = {0,0,0};
				double[] dynamicSol = new double[6];

				time1 = System.nanoTime();
				dynamicSol = dynamicAnalysis(theta, dtheta, ddtheta, T_ext, f_ext);
				//				dynamicSol = dynamicAnalysisMOD(theta, dtheta, ddtheta, T_ext, f_ext);
				time2 = System.nanoTime();

				for(int k = 0; k < dynamicSol.length; k++){
					dynSolutionVector[i][k] = dynamicSol[k];
				}

				//				if(i == 1 || i == 2){
				//					System.out.println("Torque joint 1 : " + dynamicSol[0]);
				//					System.out.println("Torque joint 2 : " + dynamicSol[1]);
				//					Matrix.displayVector(rc[0]);
				//				}


				double power = powerCalculation(dynamicSol, dtheta);


				//				System.out.println("");
				//				System.out.println("Power : " + power);

				energy += power*dTime;



				//				if(i == 70){
				//					System.out.println("Theta vector for i=" + i + ": ");
				//					Matrix.displayVector(theta);
				//					System.out.println("dTheta vector for i=" + i + ": ");
				//					Matrix.displayVector(dtheta);
				//					System.out.println("ddTheta vector for i=" + i + ": ");
				//					Matrix.displayVector(ddtheta);
				//					System.out.println("Dynamic solution vector for i=" + i + ": ");
				//					Matrix.displayVector(dynamicSol);
				//					System.out.println("The power consumed for i = " + i + " is : " + power);
				//					System.out.println("The energy consumed for i = " + i + " is : " + energy);
				//					System.out.println("Time instant i : " + trajectory.timeInstants.get(i) + " and time instant i-1 : " + trajectory.timeInstants.get(i-1));
				//					
				//				}


				for(int j = 0; j < 6; j++){
					prevTheta[j] = theta[j];
					prevDTheta[j] = dtheta[j];
				}
			}

			/*
			 * I save all the data regarding the optimal trajecotry in case they are needed.
			 * For the 3Dmap only the energy consumption would be needed but here the planner
			 * has to give more information to the robot on how to actually follow the trajectory.
			 */
			System.out.println("Combination with energy: " + energy);

			if(energy < minEnergy){

				//				System.out.println("Solution with energy: " + energy);

				minEnergy = energy;
				optDynSolution = dynSolutionVector;
				optKinSolution = solution;
			}
		}
		System.out.println("The chosen solution was the one requiring " + minEnergy + " J. ");
		EnergyPoint newPoint = new EnergyPoint(minEnergy);
		this.energyList.add(newPoint);

//		System.out.println("Last sample of inverse dynamics took : " + (time2 - time1)/1E9 );
		

//		this.writeSolutionToExcel(optKinSolution, optDynSolution, trajectory);
	}

	private double[] dynamicAnalysis(double[] theta, double[] dtheta, double[] ddtheta, double[] T, double[] f){

		/*
		 * Here we apply the Newton-Euler method to calculate the inverse dynamics of the system.
		 */



		/*
		 * alpha_1 		-> angular acceleration
		 * acc_endL 	-> acceleration of the end of the link
		 * acc_centerL 	-> acceleration of the center of the link 
		 */
		double[] dynamicSolution = new double[6];
		double[] w1 = {0,0,0}, alpha1 = {0,0,0};
		double[] accCenterL1 = {0,0,0};
		double[] accEndL1 = {0,0,0};
		double[] z0 = {0,0,1}, g0 = {0,0,-9.81};
		double[] gi;
		double[][] w = new double[6][3], alpha = new double[6][3], dw = new double[6][3], accEndL = new double[6][3];
		double[][] accCenterL = new double[6][3];

		double[] dw1 = {0,0,0};
		double[][] roti;
		double[][] r0i;

		//Forward Recursion
		for(int i = 0; i < 6; i++)
		{
			roti = (robot.hToFrom(i+1, i, theta)).getRotation();
			roti = Matrix.transpose(roti);
			//			double[] temp = Matrix.addMatrices(w1, Matrix.multiplyScalMatr(dtheta[i], z0)); // (w_0 + z0*dtheta(i))

			//FORMULA LIKE IN MATLAB
			w[i] = Matrix.addMatrices(Matrix.multiplyMatrixVector(roti, w1), Matrix.multiplyScalMatr(dtheta[i], z0));
			//

			double[] temp1 = Matrix.crossProduct(w[i], Matrix.multiplyScalMatr(dtheta[i], z0));

			for(int j = 0; j < 3; j++)
			{
				//				alpha[i][j] = alpha1[j] + z0[j] * ddtheta[i] + temp1[j];

				dw[i][j] = dw1[j] + z0[j]*ddtheta[i] + temp1[j];
			}

			//FORMULA LIKE IN MATLAB
			alpha[i] = Matrix.addMatrices(Matrix.multiplyMatrixVector(roti, alpha1), Matrix.addMatrices(temp1, Matrix.multiplyScalMatr(ddtheta[i], z0)));
			//

			//			
			//			ae(:,i)=((rot_i')*ae_1)+cross(dw(:,i),r(:,i))+cross(w(:,i),c1_2);
			//			ac(:,i)=((rot_i)'*ae_1)+cross(dw(:,i),rc(:,i))+cross(w(:,i),c1_1); 

			double[] c11 = Matrix.crossProduct(w[i], rc[i]);
			double[] c12 = Matrix.crossProduct(w[i], r[i]);

			//THESE FORMULAS ARE THE SAME AS MATLAB

			accCenterL[i] = Matrix.multiplyMatrixVector(roti, accEndL1);
			double[] temp3 = Matrix.addMatrices(Matrix.crossProduct(dw[i], rc[i]), Matrix.crossProduct(w[i], c11));
			accCenterL[i] = Matrix.addMatrices(accCenterL[i], temp3);
			accEndL[i] = Matrix.multiplyMatrixVector(roti, accEndL1);
			double[] temp4 = Matrix.addMatrices(Matrix.crossProduct(dw[i], r[i]), Matrix.crossProduct(w[i], c12));
			accEndL[i] = Matrix.addMatrices(accEndL[i], temp4);


			//FORMULAS AS ON THE BOOK
			/*
			    acc_endL[i] = Matrix.multiplyMatrixVector(rot_i, acc_endL1);
				double[] temp4 = Matrix.addMatrices(Matrix.crossProduct(dw[i], r[i]), Matrix.crossProduct(w[i], c12));
				acc_endL[i] = Matrix.addMatrices(acc_endL[i], temp4);

				//MAYBE in this formula RC[i] has to be changed to rici[i]
			 	acc_centerL[i] = Matrix.addMatrices(Matrix.addMatrices(acc_endL[i], Matrix.crossProduct(dw[i], rici[i])), Matrix.crossProduct(w[i], Matrix.crossProduct(w[i], rici[i])));
			 */




			/*
			 * Update the values for next iteration
			 */
			for(int j = 0; j < 3; j++){
				alpha1[j] = alpha[i][j];
				w1[j] = w[i][j];
				dw1[j] = dw[i][j];
				accEndL1[j] = accEndL[i][j];
				accCenterL1[j] = accCenterL[i][j];
			}
		}






		//Backward recursion


		/*
		 R0_i=translation_for_num(0,i,1,theta);
		 rot_i=translation_for_num(i,i+1,1,theta);

		 f(:,7)=fi;
		 Tr(:,7)=Ti;

		 r_ici is involved in cross products: care for dimensions!
		 */
		double[] force_i = f; //is the external force applied at the end effector. It will be used to calculate the next force in NE method.
		double[] force_next = new double[f.length];
		double[] torque_i = T;//is the external torque applied at the end effector. It will be used to calculate the next torque in NE method.
		double[] torque_next = new double[T.length];

		double[] temp1, temp2, temp3, temp4;

		for(int i = 5; i >= 0 ; i--) //in matlab the loop is 6:1
		{


			roti = (robot.hToFrom(i+2, i+1, theta)).getRotation(); //the indixes of the transf matrices are the same as in matlab though
			r0i = (robot.hToFrom(i+1, 0, theta)).getRotation();

			r0i = Matrix.transpose(r0i);

			gi = Matrix.multiplyMatrixVector(r0i, g0);

			//force_next = Mass(i)*a_center(:,i)-mass(i)*gi + rot_i*f(:,i+1);

			for(int j = 0; j < force_next.length; j++)
			{
				force_next[j] = mass[i]*(accCenterL[i][j] - gi[j]);

			}
			force_next = Matrix.addMatrices(force_next, Matrix.multiplyMatrixVector(roti, force_i));


			temp1 = Matrix.crossProduct(Matrix.multiplyMatrixVector(roti, force_i), rici[i]);
			temp2 = Matrix.crossProduct(force_next,rc[i]);
			temp3 = Matrix.crossProduct(w[i], Matrix.multiplyMatrixVector(inertiaTens.get(i), w[i]));
			temp4 = Matrix.multiplyMatrixVector(inertiaTens.get(i), alpha[i]);
			torque_next = Matrix.multiplyMatrixVector(roti, torque_i);

			//			if(i == 0){
			//				Matrix.displayVector(torque_next);
			//			}

			for(int j = 0; j < torque_next.length; j++)
			{
				torque_next[j] += temp1[j] - temp2[j] + temp3[j] + temp4[j];
			}

			//			if(i == 0){
			//				Matrix.displayVector(temp1);
			//				Matrix.displayVector(temp2);
			//				Matrix.displayVector(temp3);
			//				Matrix.displayVector(temp4);
			//				Matrix.displayVector(torque_next);
			//			}

			//Trr(i)=Tr(:,i)'*[0;0;1];
			//This calculation (from MATLAB) basically extrapulates only the third element of the 1x3 vector
			dynamicSolution[i] = torque_next[2];

			//UPDATE of Torque and Force
			for(int j = 0; j < force_i.length; j++){
				force_i[j] = force_next[j];
				torque_i[j] = torque_next[j];
			}

		}



		return dynamicSolution;
	}

	private double[] dynamicAnalysisMOD(double[] theta, double[] dtheta, double[] ddtheta, double[] T, double[] f){

		/*
		 * Here we apply the Newton-Euler method to calculate the inverse dynamics of the system.
		 */



		/*
		 * alpha_1 		-> angular acceleration
		 * acc_endL 	-> acceleration of the end of the link
		 * acc_centerL 	-> acceleration of the center of the link 
		 */
		double[] dynamicSolution = new double[6];
		double[] w1 = {0,0,0}, alpha1 = {0,0,0};
		double[] accCenterL1 = {0,0,0};
		double[] accEndL1 = {0,0,0};
		double[] z0 = {0,0,1}, g0 = {0,0,-9.81};
		double[] gi;
		double[][] w = new double[6][3], alpha = new double[6][3], dw = new double[6][3], accEndL = new double[6][3];
		double[][] accCenterL = new double[6][3];

		double[] dw1 = {0,0,0};
		double[][] roti;
		double[][] r0i;

		//Forward Recursion
		for(int i = 0; i < 6; i++)
		{
			roti = (robot.hToFrom(i+1, i, theta)).getRotation();
			roti = Matrix.transpose(roti);
			double[] temp = Matrix.addMatrices(w1, Matrix.multiplyScalMatr(dtheta[i], z0)); // (w_0 + z0*dtheta(i))
			w[i] = Matrix.multiplyMatrixVector(roti, temp); //R(i-1)_i' * (w_0 + z0*dtheta(i))


			double[] temp1 = Matrix.crossProduct(w[i], Matrix.multiplyScalMatr(dtheta[i], z0));

			for(int j = 0; j < 3; j++)
			{
				alpha[i][j] = alpha1[j] + z0[j] * ddtheta[i] + temp1[j];

				dw[i][j] = dw1[j] + z0[j]*ddtheta[i] + temp1[j];
			}
			alpha[i] = Matrix.multiplyMatrixVector(roti, alpha[i]);  


			//			
			//			ae(:,i)=((rot_i')*ae_1)+cross(dw(:,i),r(:,i))+cross(w(:,i),c1_2);
			//			ac(:,i)=((rot_i)'*ae_1)+cross(dw(:,i),rc(:,i))+cross(w(:,i),c1_1); 

			double[] c11 = Matrix.crossProduct(w[i], rici[i]);
			double[] c12 = Matrix.crossProduct(w[i], r[i]);


			//FORMULAS AS ON THE BOOK

			accEndL[i] = Matrix.multiplyMatrixVector(roti, accEndL1);
			double[] temp4 = Matrix.addMatrices(Matrix.crossProduct(alpha[i], r[i]), Matrix.crossProduct(w[i], c12));
			accEndL[i] = Matrix.addMatrices(accEndL[i], temp4);

			//MAYBE in this formula RC[i] has to be changed to rici[i]
			accCenterL[i] = Matrix.addMatrices(Matrix.addMatrices(accEndL[i], Matrix.crossProduct(alpha[i], rici[i])), Matrix.crossProduct(w[i], Matrix.crossProduct(w[i], rici[i])));





			/*
			 * Update the values for next iteration
			 */
			for(int j = 0; j < 3; j++){
				alpha1[j] = alpha[i][j];
				w1[j] = w[i][j];
				dw1[j] = dw[i][j];
				accEndL1[j] = accEndL[i][j];
				accCenterL1[j] = accCenterL[i][j];
			}
		}






		//Backward recursion


		/*
		 R0_i=translation_for_num(0,i,1,theta);
		 rot_i=translation_for_num(i,i+1,1,theta);

		 f(:,7)=fi;
		 Tr(:,7)=Ti;

		 r_ici is involved in cross products: care for dimensions!
		 */
		double[] force_i = f; //is the external force applied at the end effector. It will be used to calculate the next force in NE method.
		double[] force_next = new double[f.length];
		double[] torque_i = T;//is the external torque applied at the end effector. It will be used to calculate the next torque in NE method.
		double[] torque_next = new double[T.length];

		double[] temp1, temp2, temp3, temp4;

		for(int i = 5; i >= 0 ; i--) //in matlab the loop is 6:1
		{


			roti = (robot.hToFrom(i+2, i+1, theta)).getRotation(); //the indixes of the transf matrices are the same as in matlab though
			r0i = (robot.hToFrom(i+1, 0, theta)).getRotation();

			r0i = Matrix.transpose(r0i);

			gi = Matrix.multiplyMatrixVector(r0i, g0);

			//force_next = Mass(i)*a_center(:,i)-mass(i)*gi + rot_i*f(:,i+1);
			for(int j = 0; j < force_next.length; j++)
			{
				force_next[j] = mass[i]*accCenterL[i][j] - mass[i]*gi[j];
			}
			force_next = Matrix.addMatrices(force_next, Matrix.multiplyMatrixVector(roti, force_i));


			temp1 = Matrix.crossProduct(Matrix.multiplyMatrixVector(roti, force_i), rici[i]);
			temp2 = Matrix.crossProduct(force_next,rc[i]);
			temp3 = Matrix.crossProduct(w[i], Matrix.multiplyMatrixVector(inertiaTens.get(i), w[i]));
			temp4 = Matrix.multiplyMatrixVector(inertiaTens.get(i), alpha[i]);
			torque_next = Matrix.multiplyMatrixVector(roti, torque_i);

			//			if(i == 0){
			//				Matrix.displayVector(torque_next);
			//			}

			for(int j = 0; j < torque_next.length; j++)
			{
				torque_next[j] += temp1[j] - temp2[j] + temp3[j] + temp4[j];
			}

			//			if(i == 0){
			//				Matrix.displayVector(temp1);
			//				Matrix.displayVector(temp2);
			//				Matrix.displayVector(temp3);
			//				Matrix.displayVector(temp4);
			//				Matrix.displayVector(torque_next);
			//			}

			//Trr(i)=Tr(:,i)'*[0;0;1];
			//This calculation (from MATLAB) basically extrapulates only the third element of the 1x3 vector
			dynamicSolution[i] = torque_next[2];

			//UPDATE of Torque and Force
			for(int j = 0; j < force_i.length; j++){
				force_i[j] = force_next[j];
				torque_i[j] = torque_next[j];
			}

		}



		return dynamicSolution;
	}

	private double powerCalculation(double[] T, double[] dtheta){
		double P = 0;

		for(int i = 0; i < T.length; i++){
			P += Math.abs(T[i]*dtheta[i]);
		}

		return P;
	}


	private double[] solveTheta1(double[] Pw)
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

	private double[][] solveTheta23(double[] Pw, double theta_1 , int index)
	{

		double a2 = linkLength[1][0];
		double a4 = linkLength[3][0];
		double a5 = linkLength[4][0];

		double d1 = linkLength[0][1];
		double d3 = linkLength[2][1];

		//standard values for initialization if target is out of reach
		double[] theta_3 = {1000 , 0}; 
		double[] theta_2 = {1000, 0};

		double[][] theta_23 = new double[2][2];

		double rr;
		double a;

		if(index == 0)
		{
			rr = Math.sqrt(Math.pow((Pw[0] - a2*Math.cos(theta_1)), 2) + Math.pow((Pw[1] - a2*Math.sin(theta_1)), 2));
		}
		else
		{
			rr = - Math.sqrt(Math.pow((Pw[0] - a2*Math.cos(theta_1)), 2) + Math.pow((Pw[1] - a2*Math.sin(theta_1)), 2));
		}
		//a is the cosine of THETA 3
		a = ((-(d3*d3) - (a4 + a5)*(a4 + a5) + (Pw[2] - d1)*(Pw[2] - d1) + rr*rr ))/(2*d3*(a4 + a5));

		//		System.out.println(" Questa cazzo di 'rr' vale: " + rr );
		//		System.out.println(" Questa cazzo di 'a' vale: " + a );

		if(a < -1 || a > 1 )
		{
			theta_23[0] = theta_2;
			theta_23[1] = theta_3;
			//			System.out.println("Impossible to reach that target.");
			return theta_23;
		}
		else
		{
			double temp = Math.sqrt(1 - a*a);
			theta_3[0] = Math.atan2(temp, a);
			theta_3[1] = Math.atan2(-temp, a);

			theta_23[1] = theta_3;
		}
		for(int i = 0; i < theta_3.length; i++)
		{
			theta_2[i] = Math.atan2((Pw[2] - d1),rr) + Math.atan2((a4 + a5)*Math.sin(theta_3[i]), d3 + (a4 + a5)*Math.cos(theta_3[i])); 
		}

		theta_23[0] = theta_2;		
		return theta_23;

	}

	private double[][] solveTheta456(Target target_position, double theta_1, double theta_2, double theta_3, int point_i){

		/*
		 * TODO
		 * 1: get the target rotation
		 * 2: get from the robot the rotational matrix from 0 to TOOL
		 * 3: build the matrix from which the Euler's angles will be calculated
		 * 
		 * IMPORTANT: needed a way to read the joint angles of the previous configuration 
		 * HOW DO I GET THAT?
		 */



		double[][] R_zyz = new double[3][3];
		double[] theta_values = {theta_1, theta_2, theta_3, 0, 0 , 0};
		Target T0_3 = robot.hToFrom(3, 0, theta_values);
		double[][] R0_3 = T0_3.getRotation();
		double[][] R33 = {
				{0, 0, 1},
				{0, 1, 0},
				{-1, 0, 0}

		};
		double[][] R6tcp = {
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}

		};
		double[][] R0tcp = target_position.getRotation();


		//Rzyz=(R0_3*R3_3)'*R0_TCP*R6_TCP';
		R_zyz = Matrix.multiplyMatrices(R0tcp, Matrix.transpose(R6tcp));
		R_zyz = Matrix.multiplyMatrices(Matrix.transpose(Matrix.multiplyMatrices(R0_3, R33)), R_zyz);


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
			if(-0.0001 < theta_5[i] && theta_5[i] < 0.0001){

				//				System.out.println("There is a Singularity in : " + point_i);

				//theta_5 is 0 and we're in a singularity case

				if(point_i == 0) {
					//previous_theta4_TEST = 0;
					theta_4[i] = 0;
				}
				else { 
					//previous_theta4_TEST = theta4Optimal(theta_5[i], theta_1, theta_2, theta_3, point_i);
					theta_4[i] = theta4Optimal(theta_5[i], theta_1, theta_2, theta_3, point_i);
				}

				//theta_4[i] = previous_theta4_TEST;
				theta_46 = Math.atan2(R_zyz[1][0], R_zyz[0][0]);

				theta_6[i] = theta_46 - theta_4[i];

			}
			else
			{//theta_5 is NOT 0 and we can calculate the other angles
				temp = 1d/Math.sin(theta_5[i]);
				theta_4[i] = Math.atan2(R_zyz[0][2]*temp, -R_zyz[1][2]*temp);
				theta_6[i] = Math.atan2(R_zyz[2][0]*temp, R_zyz[2][1]*temp);
			}
		}

		double[][] theta_456 = new double[3][2];
		theta_456[0] = theta_4;
		theta_456[1] = theta_5;
		theta_456[2] = theta_6;

		//DEBUG LINE
		//		if(point_i == 0){
		//			System.out.println("R_zyz :");
		//			Matrix.displayMatrix(R_zyz);
		//			System.out.println("R0_3 : ");
		//			Matrix.displayMatrix(R0_3);
		//			System.out.println("Theta1: " + theta_1 + "  theta2: " + theta_2 + "  theta3: " + theta_3);
		//			System.out.println("Theta4: " + theta_4[0] + " or " + theta_4[1]);
		//			System.out.println("Theta5: " + theta_5[0] + " or " + theta_5[1]);
		//			System.out.println("Theta6: " + theta_6[0] + " or " + theta_6[1]);
		//			
		//		}



		return theta_456;
		/*
		 * REMINDER
		 * in case performance change significantly it could be possible to assign
		 * the theta values directly to the 3x2 matrix and then return it.
		 */
	}

	private double theta4Optimal(double theta5, double theta1, double theta2, double theta3, int index){
		/*
		 * Theta4 of the previous point (i-1) could be any of the 8 previous configurations.
		 * Therefore it is chosen the angle belonging to the configuration that covers less space 
		 * between i-1 and i in terms of theta-1,2,3.
		 * 
		 * We ALREADY know that theta5 is a SINGULARITY (very close to zero)
		 */
		double minSum = Double.MAX_VALUE;
		double temp = 0;
		int configIndex = 0;

		if (theta5 < 0){
			for(int i = 0; i < thetaM[0].length; i++){
				if(thetaM[0][i][index-1] == 0){
					continue;
				}
				if(thetaM[4][i][index-1] < 0){
					//theta_sing=[theta_sing,abs(thetaM(1:3,i,t-1)-[theta1;theta2;theta3])];
					temp = Math.abs(thetaM[0][i][index-1] - theta1) + Math.abs(thetaM[1][i][index-1] - theta2) + Math.abs(thetaM[2][i][index-1] - theta3);
					if(temp < minSum){
						minSum = temp;
						configIndex = i;
					}
				}


			}

		} else {
			if (theta5 > 0){

				for(int i = 0; i < thetaM[0].length; i++){
					if(thetaM[0][i][index-1] == 0){
						continue;
					}
					if(thetaM[4][i][index-1] > 0){
						//theta_sing=[theta_sing,abs(thetaM(1:3,i,t-1)-[theta1;theta2;theta3])];
						temp = Math.abs(thetaM[0][i][index-1] - theta1) + Math.abs(thetaM[1][i][index-1] - theta2) + Math.abs(thetaM[2][i][index-1] - theta3);
						if(temp < minSum){
							minSum = temp;
							configIndex = i;
						}
					}
				}

			} else {

				for(int i = 0; i < thetaM[0].length; i++){
					if(thetaM[0][i][index-1] == 0){
						continue;
					}

					temp = Math.abs(thetaM[0][i][index-1] - theta1) + Math.abs(thetaM[1][i][index-1] - theta2) + Math.abs(thetaM[2][i][index-1] - theta3);
					if(temp < minSum){
						minSum = temp;
						configIndex = i;
					}
				}

			}

		}

		return thetaM[3][configIndex][index-1];
	}

	private double[][][] applyLimitation(double[][][] ThetaM, int point_i)
	{

		for(int i = 0; i < ThetaM[0].length ; i++)
		{
			for(int j = 0; j < ThetaM.length ; j++)
			{
				if(ThetaM[j][i][point_i] < this.jointLimits[j][0] || ThetaM[j][i][point_i] > this.jointLimits[j][1])
				{
					//IF a joint exceeds its limits, the all configuration is compromised, hence 
					//it will be put to zero.


					//System.out.println("Il giunto " + j + " oltrepassa i suoi limiti");

					for(int k = 0; k < ThetaM.length ; k++ )
					{
						ThetaM[k][i][point_i] = 0;
					}

					break; //no need to check the other joints in this very configuration.
				}
			}
		}
		return ThetaM;
	}


	private void configIteration(double[][][] ThetaM, int pointi, int prevConfig, ArrayList<Double> trajectoryTimes, int[] currentSolution, ArrayList<int[]> solutions){
		/*
		 * This function should check if the path can still be followed even if it doesn't 
		 * exceed the joints limitations. The reason is because the path could be too fast for
		 * the robot, and in that case the sequence of configuration chosen has to be discarded.
		 */
		int N = ThetaM[0][0].length;
		int error = 0;


		for(int currentConfig = 0; currentConfig < ThetaM[0].length; currentConfig++)
		{
			if (pointi == 0)
			{
				error = 0;
			}
			else
			{
				if(ThetaM[0][currentConfig][pointi] == 0 && ThetaM[1][currentConfig][pointi] == 0)
				{
					//In this case the configuration has already been put to ZERO from Apply_limitation
					//so there's no reason to check it again
					error  = 1;
				}
				else{
					error = checkConfiguration(currentConfig, prevConfig, pointi, ThetaM, trajectoryTimes);
				}
			}


			if(error == 0) 
			{
				//if the configuration checked was OK and with no problems then we can add it
				//to the solution vector.
				currentSolution[pointi] = currentConfig;
				if (pointi < N-1 )//N-2 or N-1?? pointi assumes values 0:134=N-1
				{
					configIteration(ThetaM, pointi+1, currentConfig, trajectoryTimes, currentSolution, solutions);
				}
				else
				{
					int[] newSolution = new int[currentSolution.length];

					for(int k = 0; k < newSolution.length; k++){
						newSolution[k] = currentSolution[k];
					}

					solutions.add(newSolution);
					//					System.out.println("New solution found, the number of solutions is : " + solutions.size());

				}

			}

		}

	}

	private int checkConfiguration(int currentConfig, int prevConfig, int point, double[][][] ThetaM, ArrayList<Double> trajectoryTimes){
		int error = 0;
		//point_i will be for sure >= 1
		double dTime = trajectoryTimes.get(point) - trajectoryTimes.get(point - 1);
		double[] wMax = {this.speedLimits[0], this.speedLimits[1], this.speedLimits[2], this.speedLimits[3], this.speedLimits[4], this.speedLimits[5]};
		double[] thMax = new double[wMax.length];

		for(int i = 0; i< wMax.length; i++){ 
			wMax[i] = wMax[i] * Math.PI/180;
			thMax[i] = wMax[i] * dTime;
		}
		double[] thAssump = thMax;
		double deltaJoint;

		for(int j = 0 ; j < ThetaM.length ; j++)
		{
			deltaJoint = Math.abs(ThetaM[j][currentConfig][point] - ThetaM[j][prevConfig][point-1]);

			if((j <= 2 || j == 4) && (deltaJoint > thAssump[j]))
			{
				return error = 1;
			}

			if((j == 3 || j == 5) && (deltaJoint > thAssump[j]))
			{
				//				if(j == 5){
				//					System.out.println("theta 6 exceeds its limits in i = " + point_i);
				//				}

				double delta_jointp = Math.abs(ThetaM[j][currentConfig][point] + 2*Math.PI - ThetaM[j][prevConfig][point-1]);
				double delta_jointm = Math.abs(ThetaM[j][currentConfig][point] - 2*Math.PI - ThetaM[j][prevConfig][point-1]);

				if(delta_jointp < thAssump[j]) { ThetaM[j][currentConfig][point] += 2*Math.PI; }
				else
				{
					if(delta_jointm < thAssump[j]){ ThetaM[j][currentConfig][point] -= 2*Math.PI;}
					else
					{
						return error = 1;
					}
				}
			}	
		}	
		return error;
	}

	private void copyRobotParameters(){

		this.linkLength = robot.getLinkLength();
		this.jointLimits = robot.getJointLimits();
		this.speedLimits = robot.getSpeedLimits();
		this.mass = robot.getLinkMasses();
		this.inertiaTens = robot.getInertiaTens();
		this.rc = robot.getRc();
		this.rici = robot.getRici();
		this.r = robot.getR();

		this.acceleration = 1.0;
		this.velocity = 0.2;	

	}



	public void displayThetaM(int index){
		//System.out.println(thetaM.length + " " + thetaM[0].length + " " + thetaM[0][0].length);
		System.out.println("");
		System.out.println("ThetaM of index " + index);

		for(int i = 0; i < thetaM[0][0].length; i++){
			if(i != index){
				continue;
			}

			for(int j = 0; j < thetaM.length; j++){
				for(int k = 0; k < thetaM[0].length; k++){
					System.out.print(thetaM[j][k][i] + " ");
				}
				System.out.println(" ");
			}
		}
		System.out.println(" ");
	}

	public double[][][] getThetaM(){
		return this.thetaM;
	}

	public List<EnergyPoint> getEnergyList(){
		return this.energyList;
	}

	public void setAcceleration(double acceleration){
		this.acceleration = acceleration;
	}

	public void setVelocity(double velocity){
		this.velocity = velocity;
	}


}
