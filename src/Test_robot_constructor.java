import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.io.SAXReader;
import org.dom4j.*;
import org.jaxen.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class Test_robot_constructor {

	public static void main(String[] args) {
		String path = "C:\\Users\\Enrico\\Google Drive\\Tesi\\Prova Costruttore Robot\\ABB IRB140.txt"; //percorso del file con i dati
		//String path = "C:\\Users\\Enrico\\Desktop\\MatriceProva.txt";
		
		Robot newRobot = new Robot(path);
		int DOF = newRobot.getDOF();
		System.out.println("The degree of freedom of the robot are: " + DOF); //checked and working
		System.out.println("The links lengths are the following: ");
		double[][] linkLength = newRobot.getParameters("link length");
		for(int i=0; i<linkLength.length; i++) //checked and working
		{
			for(int j=0;j<linkLength[0].length;j++)
			{
				System.out.print(linkLength[i][j] + " ");
			}
			System.out.println("");
		}
		
		
		double[][] jointLimit = newRobot.getParameters("joint limits");
		for(int i=0; i<jointLimit.length; i++) //checked and working
		{
			for(int j=0;j<jointLimit[0].length;j++)
			{
				System.out.print(jointLimit[i][j]/Math.PI + " ");
			}
			System.out.println("");
		}
		
		
		double[][] COG = newRobot.getParameters("center of gravity");
		for(int i=0; i<COG.length; i++) //checked and working
		{
			for(int j=0;j<COG[0].length;j++)
			{
				System.out.print(COG[i][j] + " ");
			}
			System.out.println("");
		}
		
		/*
		 * Trying to test the path class
		 */
		double[][] matrix1 = {
				{0.0, 0.0, 1.0, 0.59627},
				{1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.66282},
				{0.0, 0.0, 0.0, 1.0}
		};
		double[][] matrix2 = {
				{0.0, 0.0, 1.0, 0.59627},
				{1.0, 0.0, 0.0, -0.21132},
				{0.0, 1.0, 0.0, 0.58041},
				{0.0, 0.0, 0.0, 1.0}
		};
//		double[][] matrix1 = {
//				{0.0, 0.0, 1.0, 1},
//				{1.0, 0.0, 0.0, 4},
//				{0.0, 1.0, 0.0, -2},
//				{0.0, 0.0, 0.0, 1.0}
//		};
//		double[][] matrix2 = {
//				{0.0, 0.0, 1.0, 0},
//				{1.0, 0.0, 0.0, 2},
//				{0.0, 1.0, 0.0, 11},
//				{0.0, 0.0, 0.0, 1.0}
//		};
		
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;
		
		Target testTarget1 = new Target(matrix1);
		Target testTarget2 = new Target(matrix2);
		Path test_path = new Path(testTarget1, testTarget2, t_sample, x_sample);
		
		test_path.setMaxAcc(acceleration);
		test_path.setMaxVel(velocity);
		test_path.interpolate();
		Trajectory test_output = test_path.getTrajectory();
		System.out.println("Number of points : " + test_output.points.size());
		
		/*
		for(Target i : test_output.Points)
		{
			double[] position_vector = i.getPosition();
			
			System.out.println("Punto : " + position_vector[0] + " " + position_vector[1] + " " + position_vector[2]);
		}
		*/
		
		//double[][][] ThetaM = new double[2][3][4];
		
		//System.out.println("La seconda dimensione di ThetaM è : " + ThetaM[0].length);
		
		/*
		//Testing Target methods
		
		double[][] matrix3 = {
				{1.0, Math.PI, 0.0, 4.6},
				{2.1, 1.0, 0.0, 4.0},
				{2.45, 0.0, 1.0, -1.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		Target test_target3 = new Target(matrix3);
		double[][] test = test_target3.getInvHomMatrix();
		/*
		double[][] test = {
				{1.0, 0.0, 0.0, 0},
				{0, 1.0, 0.0, 0},
				{0, 0.0, 1.0, 0},
				{0.0, 0.0, 0.0, 1}
		};
		*/
		/*
		double[][] test_ident = MultiplyMatrices(test,matrix3);
		
		System.out.println("Righe : " + test_ident.length + " Colonne : " + test_ident[0].length);
		
		for(int i=0; i< test_ident.length;i++)
		{
			for(int j=0; j<test_ident[0].length; j++)
			{
				System.out.print(test_ident[i][j] + " ");
			}
			System.out.println(" ");
		}
		*/
		
		
		
		/*---------------------
		 * XML TESTING PART
		 * --------------------
		 */
		try {
			 String FilePath = "C:\\Users\\Enrico\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML Costruttore Robot.xml";
 	         File inputFile = new File(FilePath);
	         SAXReader reader = new SAXReader();
	         Document document = reader.read( inputFile );
	         
	         
	         System.out.println("Root element : " + document.getRootElement().getName());
	         String rootName = document.getRootElement().getName();
	         
	         Element root = document.getRootElement();
	         
	         List<Node> nodes = document.selectNodes(rootName + "/link");
	         
	         
	         for(Node node : nodes)
	         {
	        	// int linkNumber = Integer.parseInt(node.selectSingleNode("number").getText());
	        	 int linkNumber = (int) Double.parseDouble(node.selectSingleNode("number").getText());
	        	 double linkMass = Double.parseDouble(node.selectSingleNode("mass").getText());
	        	 System.out.println("\nCurrent Element : " + node.getName() + " " + linkNumber);
	        	 System.out.println("Mass = " + linkMass );
	        	 
	         }
	         
	         int DoF =(int) Double.parseDouble(document.selectSingleNode(rootName + "/DOF").getText());
	         System.out.println("Degrees of Freedom = " + DoF);
	         
		}
		catch (DocumentException e) {
	         e.printStackTrace();
	      }
		
		
		/*---------------------
		 * XML ROBOT CONSTRUCTOR TESTING PART
		 * --------------------
		 */
		 String FilePathXML = "C:\\Users\\Enrico\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML Costruttore Robot.xml";
		
		Robot newRobotXML = new Robot(FilePathXML, true);
		
		DOF = newRobotXML.getDOF();
		System.out.println("The degree of freedom of the robot are: " + DOF); //checked and working
		System.out.println("The links lengths are the following: ");
		linkLength = newRobotXML.getParameters("link length");
		for(int i=0; i<linkLength.length; i++) //checked and working
		{
			for(int j=0;j<linkLength[0].length;j++)
			{
				System.out.print(linkLength[i][j] + " ");
			}
			System.out.println("");
		}
		
		
		jointLimit = newRobotXML.getJointLimits();
		for(int i=0; i<jointLimit.length; i++) //checked and working
		{
			for(int j=0;j<jointLimit[0].length;j++)
			{
				System.out.print(jointLimit[i][j]/Math.PI + " ");
			}
			System.out.println("");
		}
		
		
		COG = newRobotXML.getCoGMatrix();
		for(int i=0; i<COG.length; i++) //checked and working
		{
			for(int j=0;j<COG[0].length;j++)
			{
				System.out.print(COG[i][j] + " ");
			}
			System.out.println("");
		}
		
//		double[][] R33 = {
//				{0, 0, 1},
//				{0, 1, 0},
//				{-1, 0, 0}
//				
//		};
//		double[][] R6tcp = {
//				{1, 0, 0},
//				{0, 1, 0},
//				{0, 0, 1}
//				
//		};
//		double[][] R0tcp = test_output.points.get(3).getRotation();
//		double[][] R_zyz = Matrix.multiplyMatrices(R0tcp, Matrix.transpose(R6tcp));
//		R_zyz = Matrix.transpose(Matrix.multiplyMatrices(R_zyz, R33));
//		
//		
//		System.out.println("Righe : " + R_zyz.length + " Colonne : " + R_zyz[0].length);
//		
//		for(int i=0; i< R_zyz.length;i++)
//		{
//			for(int j=0; j<R_zyz[0].length; j++)
//			{
//				System.out.print(R_zyz[i][j] + " ");
//			}
//			System.out.println(" ");
//		}
		
		Target[] targets = new Target[2];
		targets[0] = testTarget1;
		targets[1] = testTarget2;
		
		OnlinePlanner dynSolution = new OnlinePlanner(targets, newRobotXML);
		
		dynSolution.run();
		
		
//		Path newPath = new Path(targets[0], targets[1], t_sample, x_sample);
//		newPath.setMaxAcc(acceleration);
//		newPath.setMaxVel(velocity);
//		newPath.interpolate();
//		
//		for(int i = 0; i < newPath.getTrajectory().points.size(); i++){
//			Matrix.displayVector(newPath.getTrajectory().points.get(i).getPosition());
//		}
//		System.out.println(" ");
//		for(int i = 0; i < newPath.getTrajectory().timeInstants.size(); i++){
//			System.out.println(newPath.getTrajectory().timeInstants.get(i));
//		}
//		
//		System.out.println("Number of points : " + newPath.getTrajectory().points.size() + " and time istants: " + newPath.getTrajectory().timeInstants.size());
//		
//		System.out.println("");
//		System.out.println("");
//		double theta1Test = -3.1416;
//		double theta2Test = 0.6655;
//		double theta3Test = 2.4033;
//		double[] jointValues = {theta1Test,theta2Test,theta3Test,0,0,0};
//		
//		double[][] R0_3 = newRobotXML.Hto_from(3, 0, jointValues).getRotation();
//		
//		Matrix.displayMatrix(R0_3);
//		double[] thetaValues = {0,0,0,0,0,0}; 
//		
//		System.out.println(" ");
//		Matrix.displayMatrix(newRobotXML.Hto_from(4, 3,thetaValues).getHomMatrix());
//		
//		double[] joint_values = {0,0,0,0,0,0};
//		linkLength = newRobotXML.getLinkLength();
//		
//		double[][] H0_1 = {
//				{Math.cos(joint_values[0]), -Math.sin(joint_values[0]), 0, 0},
//				{Math.sin(joint_values[0]), Math.cos(joint_values[0]), 0, 0},
//				{0, 0, 1, linkLength[0][1]},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H0_1);
//		System.out.println("");
//		double[][] H1_2 = {
//				{Math.cos(joint_values[1]), -Math.sin(joint_values[1]), 0, linkLength[1][0]},
//				{0, 0, -1, 0},
//				{Math.sin(joint_values[1]), Math.cos(joint_values[1]), 0, 0},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H1_2);
//		System.out.println("");
//		double[][] H2_3 = {
//				{Math.cos(joint_values[2]), -Math.sin(joint_values[2]), 0, 0},
//				{Math.sin(joint_values[2]), Math.cos(joint_values[2]), 0, linkLength[2][1]},
//				{0, 0, 1, 0},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H2_3);
//		System.out.println("");
//		double[][] H3_4 = {
//				{0, 0, 1, linkLength[3][0] + linkLength[4][0]},
//				{Math.sin(joint_values[3]), Math.cos(joint_values[3]), 0, 0},
//				{-Math.cos(joint_values[3]), Math.sin(joint_values[3]), 0, 0},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H3_4);
//		System.out.println("");
//		double[][] H4_5 = {
//				{0, 0, -1, 0},
//				{Math.sin(joint_values[4]), Math.cos(joint_values[4]), 0, 0},
//				{Math.cos(joint_values[4]), -Math.sin(joint_values[4]), 0, 0},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H4_5);
//		System.out.println("");
//		double[][] H5_6 = {
//				{0, 0, 1, 0},
//				{Math.sin(joint_values[5]), Math.cos(joint_values[5]), 0, 0},
//				{-Math.cos(joint_values[5]), Math.sin(joint_values[5]), 0, 0},
//				{0, 0, 0, 1}
//		};
//		
//		Matrix.displayMatrix(H5_6);
//		System.out.println("");
//		double[][] H6_endEff = {
//				{1, 0, 0, 0},
//				{0, 1, 0, 0},
//				{0, 0, 1, 0},
//				{0, 0, 0, 1}
//		};
		
		
//		Matrix.displayMatrix(H6_endEff);
		
		
		/*
		 * TESTING THE DYNAMIC ANALYSIS
		 */
		
		double[] prev_theta = {0.0, -0.23355842765405876, 0.12967262418466752, 3.141592653589793, -0.10388580346939125, 3.141592653589793};
		double[] prev_dtheta = {0, 0 , 0, 0 , 0, 0};
		double[] dtheta = {-0.00876825075922701, -5.284036864505737E-4, -0.0041616514068643085, -0.08451675027405514, -0.004693740419386971, 0.08406068565993685};
		double[] ddtheta = {-0.8768250759227009, -0.05284036864505737, -0.41616514068643085, -8.451675027405514, -0.46937404193869714, 8.406068565993685};
		
		double[] w_1 = {0,0,0}, alpha_1 = {0,0,0};
		double[] acc_centerL1 = {0,0,0};
		double[] acc_endL1 = {0,0,0};
		double[] z0 = {0,0,1}, g0 = {0,0,-9.81};
		double[] gi;
		double[][] w = new double[6][3], alpha = new double[6][3], dw = new double[6][3], acc_endL = new double[6][3];
		double[][] acc_centerL = new double[6][3];
		
 		double[] dw_1 = {0,0,0};
		double[][] rot_i;
		double[][] R0_i;
		
		double[][] rici = newRobotXML.getRici();
		double[][] r = newRobotXML.getR();
		double[][] rc = newRobotXML.getRc();
		double[] mass = newRobotXML.getLinkMasses();
		List<double[][]> inertiaTens = newRobotXML.getInertiaTens();
		
		double[][][] thetaM = dynSolution.getThetaM();
		
		double[] dynamic_solution = {0,0,0,0,0,0};
		
		double[] theta = {-8.768250759227009E-5, -0.23356371169092327, 0.12963100767059887, 3.1407474860870526, -0.10393274087358512, 3.1424332604463925}; 
		
		for(int i = 0; i < theta.length; i++){
//			theta[i] = thetaM[]
		}
		
		for(int i = 0; i < 6; i++)
		{
			rot_i = (newRobotXML.hToFrom(i+1, i, theta)).getRotation();
			rot_i = Matrix.transpose(rot_i);
			
//			System.out.println("Rot_i transpose :");
//			Matrix.displayMatrix(rot_i);
			
//			double[] temp = Matrix.addMatrices(w_1, Matrix.multiplyScalMatr(dtheta[i], z0)); // (w_0 + z0*dtheta(i))
//			w[i] = Matrix.multiplyMatrixVector(rot_i, temp); //R0_1' * (w_0 + z0*dtheta(i))
			double[] temp = Matrix.multiplyScalMatr(dtheta[i], z0);
			
//			System.out.println("Temporary vector (z0*dtheta) :");
//			Matrix.displayVector(temp);
			
			//FORMULA LIKE IN MATLAB
			w[i] = Matrix.addMatrices(Matrix.multiplyMatrixVector(rot_i, w_1), Matrix.multiplyScalMatr(dtheta[i], z0));
//			Matrix.displayVector(w[i]);
			//
			
			double[] temp1 = Matrix.crossProduct(w[i], Matrix.multiplyScalMatr(dtheta[i], z0));
			
//			System.out.println("Cross product vector:");
//			Matrix.displayVector(temp1);
			
			for(int j = 0; j < 3; j++)
			{
				alpha[i][j] = alpha_1[j] + z0[j] * ddtheta[i] + temp1[j];
				
				dw[i][j] = dw_1[j] + z0[j]*ddtheta[i] + temp1[j];
			}
			alpha[i] = Matrix.multiplyMatrixVector(rot_i, alpha[i]);  
			
			//FORMULA LIKE IN MATLAB
			alpha[i] = Matrix.addMatrices(Matrix.multiplyMatrixVector(rot_i, alpha_1), Matrix.addMatrices(temp1, Matrix.multiplyScalMatr(ddtheta[i], z0)));
			//
//			System.out.println("Alpha vector:");
//			Matrix.displayVector(alpha[i]);
			
			
			double[] c11 = Matrix.crossProduct(w[i], rc[i]);
			double[] c12 = Matrix.crossProduct(w[i], r[i]);
			
			acc_centerL[i] = Matrix.multiplyMatrixVector(rot_i, acc_endL1);
			double[] temp3 = Matrix.addMatrices(Matrix.crossProduct(dw[i], rc[i]), Matrix.crossProduct(w[i], c11));
			acc_centerL[i] = Matrix.addMatrices(acc_centerL[i], temp3);
			acc_endL[i] = Matrix.multiplyMatrixVector(rot_i, acc_endL1);
			double[] temp4 = Matrix.addMatrices(Matrix.crossProduct(dw[i], r[i]), Matrix.crossProduct(w[i], c12));
			acc_endL[i] = Matrix.addMatrices(acc_endL[i], temp4);
			
			
//			System.out.println("cross product c11:");
//			Matrix.displayVector(c11);
//			System.out.println("cross product c12:");
//			Matrix.displayVector(c12);
			
//			System.out.println("dw vector:");
//			Matrix.displayVector(dw[i]);
//			System.out.println("Acceleration center link:");
//			Matrix.displayVector(acc_centerL[i]);
//			System.out.println("Acceleration end link:");
//			Matrix.displayVector(acc_endL[i]);
			
			/*
			 * Update the values for next iteration
			 */
		
			for(int j = 0; j < 3; j++){
				alpha_1[j] = alpha[i][j];
				w_1[j] = w[i][j];
				dw_1[j] = dw[i][j];
				acc_endL1[j] = acc_endL[i][j];
				acc_centerL1[j] = acc_centerL[i][j];
				}
		}
		
		//Backward recursion
		double[] T_ext = {0,0,0};
		double[] f_ext = {0,0,0};
		
				/*
				 R0_i=translation_for_num(0,i,1,theta);
				 rot_i=translation_for_num(i,i+1,1,theta);
				 
				 f(:,7)=fi;
				 Tr(:,7)=Ti;
				 
				 r_ici is involved in cross products: care for dimensions!
				 */
				double[] force_i = f_ext; //is the external force applied at the end effector. It will be used to calculate the next force in NE method.
				double[] force_next = new double[f_ext.length];
				double[] torque_i = T_ext;//is the external torque applied at the end effector. It will be used to calculate the next torque in NE method.
				double[] torque_next = new double[T_ext.length];
				
				double[] temp1, temp2, temp3, temp4;
				
				for(int i = 5; i >= 0 ; i--) //in matlab the loop is 6:1
				{
					
					
					rot_i = (newRobotXML.hToFrom(i+2, i+1, theta)).getRotation(); //the indixes of the transf matrices are the same as in matlab though
					R0_i = (newRobotXML.hToFrom(i+1, 0, theta)).getRotation();
					
					
//					System.out.println("Rot_i matrix :");
//					Matrix.displayMatrix(rot_i);
//					System.out.println("R0_i matrix :");
//					Matrix.displayMatrix(R0_i);
					
					
					
					R0_i = Matrix.transpose(R0_i);
					
					gi = Matrix.multiplyMatrixVector(R0_i, g0);
				
					//force_next = Mass(i)*a_center(:,i)-mass(i)*gi + rot_i*f(:,i+1);
					for(int j = 0; j < force_next.length; j++)
					{
						force_next[j] = mass[i]*acc_centerL[i][j] - mass[i]*gi[j];
					}
					force_next = Matrix.addMatrices(force_next, Matrix.multiplyMatrixVector(rot_i, force_i));
					
//					System.out.println("force_next vector:");
//					Matrix.displayVector(force_next);
					
//					Tr(:,i)=rot_i*Tr(:,i+1)+cross(rot_i*f(:,i+1),rici(:,i))-cross(f(:,i),rc(:,i))+cross(w(:,i),Iner(:,:,i)*w(:,i))+Iner(:,:,i)*alpha(:,i);

					
					temp1 = Matrix.crossProduct(Matrix.multiplyMatrixVector(rot_i, force_i), rici[i]);
					temp2 = Matrix.crossProduct(force_next,rc[i]);
					temp3 = Matrix.crossProduct(w[i], Matrix.multiplyMatrixVector(inertiaTens.get(i), w[i]));
					temp4 = Matrix.multiplyMatrixVector(inertiaTens.get(i), alpha[i]);
					torque_next = Matrix.multiplyMatrixVector(rot_i, torque_i);
					
//					System.out.println("temp1:");
//					Matrix.displayVector(temp1);
//					System.out.println("temp2:");
//					Matrix.displayVector(temp2);
//					System.out.println("Part of the cross product:");
//					Matrix.displayVector(Matrix.multiplyMatrixVector(inertiaTens.get(i), w[i]));
//					System.out.println("w[i] vector:");
//					Matrix.displayVector(w[i]);
//					System.out.println("inertia tensor "+ (i+1) +":");
//					Matrix.displayMatrix(inertiaTens.get(i));
//					System.out.println("temp3:");
//					Matrix.displayVector(temp3);
//					System.out.println("temp4:");
//					Matrix.displayVector(temp4);
					
					
					
					for(int j = 0; j < torque_next.length; j++)
					{
						torque_next[j] += temp1[j] - temp2[j] + temp3[j] + temp4[j];
					}
					
//					System.out.println("torque_next vector:");
//					Matrix.displayVector(torque_next);
					
					//Trr(i)=Tr(:,i)'*[0;0;1];
					//This calculation (from MATLAB) basically extrapulates only the third element of the 1x3 vector
					dynamic_solution[i] = torque_next[2];
					
					//UPDATE of Torque and Force
					for(int j = 0; j < force_i.length; j++){
						force_i[j] = force_next[j];
						torque_i[j] = torque_next[j];
					}
				}
//				System.out.println("Solution vector:");
//				Matrix.displayVector(dynamic_solution);
		
	}
	
	public static double[][] MultiplyMatrices(double[][] left, double[][] right){

		/*
		 * Function implemented to multiply 2 matrices (useful in Hfrom_to)
		 */

		double[][] M = new double[left.length][right[0].length];
		
		for(int row = 0; row < M.length; row++)
			for(int column = 0; column < M[0].length; column++)
				M[row][column] = 0;
		

		if(left[0].length == right.length){

			for(int row = 0; row < M.length; row++)
			{
				for(int column = 0; column < M[0].length; column++)
				{
					for(int i = 0; i < left[0].length; i++)
					{
					M[row][column] += left[row][i] * right[i][column];
					}
				}
		
			}

		}else{
			System.out.println("The dimension of the matrices do not agree. [m*n n*p = m*p]");
		}

		return M;
	}
	
	
}
