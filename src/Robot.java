import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.dom4j.io.SAXReader;
import org.dom4j.*;
import org.jaxen.*;
import java.io.File;

/*
 * TODO it defines all the parameters of the robot
 *  - use XML standard to define the constructors and the data
 * 
 */

public class Robot {

	private String model;
	private int dof;
	private double[][] linkLength;
	private double[] linkMasses;
	private List<double[][]> inertiaTens;
	private double[][] jointLimits;
	private double[][] cogMatrix; //center of gravity matrix (SolidWorks)
	private double[][] rc; //it has to be structured as 6x3 matrix to work down in the dynamic
	private double[][] rici; //6x3 matrix as well. The reason is because in loops 1x3 vectors are needed and will be taken by having r_c[i] instructions
	Target location;

	public static void main(String[] args){

	}

	/* 
	 * you insert the model, the function looks for it in the cataloge.
	 * If it finds it, it initializes the robot with the data found
	 * in the cataloge, otherwise it gives you 2 choices: you can either
	 * get the data from a file or insert them manually. (?)
	 */



	//Part of CONSTRUCTORS (if more than one are needed)

	public Robot(String model_path){

		String section_limiter = "**--";
		String comment_limiter = "%%";
		String s;

		int section = 0;

		try{
			FileReader f;
			f = new FileReader(model_path);
			BufferedReader b;
			b = new BufferedReader(f);

			while(true){
				s = b.readLine();
				if(s == null)
				{
					break;
				}

				if(s.startsWith(comment_limiter) || s.isEmpty())
				{
					continue;
				}
				else
				{
					if(section == 0)
					{
						this.model = s;
					}
					if(s.startsWith(section_limiter))
					{

						notifyActualSection(section);
						System.out.println(s);
						copyRobotParameters(b, section);
						section++;
					}

				}

			}
		}
		catch (IOException e){
			System.out.println("Error in opening the file");
		}

	}
	
	public Robot(String xml_FilePath, boolean xml_true){
		
		xmlBuilder(xml_FilePath);
		jointDistances();
	}
	

	//END part of CONSTRUCTORS




	//Part of GET methods



	public int getDOF(){
		return this.dof;
	}

	public double[][] getLinkLength() {
		return linkLength;
	}

	public double[][] getJointLimits() {
		return jointLimits;
	}

	public double[][] getParameters(String arg)
	{

		/*
		 * Function to get the parameters of the robot.
		 * Consider using only some keywords no matter the overall string given in input maybe
		 * I don't think it would be better having a single method for every different parameter
		 * just because they are all of the same type;
		 */

		double[][] defaultRes = new double[0][0];

		if(arg.equalsIgnoreCase("link length") || arg.equalsIgnoreCase("link_length"))
		{
			return this.linkLength;
		}
		if(arg.equalsIgnoreCase("joint limits") || arg.equalsIgnoreCase("joint_limits"))
		{
			return this.jointLimits;
		}
		if(arg.equalsIgnoreCase("center of gravity"))
		{
			return this.cogMatrix;
		}

		return defaultRes; 
	}
	
	public List<double[][]> getInertiaTens(){
		return this.inertiaTens;
	}
	
	public double[][] getSingleInertiaTens(int index){
		return this.inertiaTens.get(index);
	}

	public double[] getLinkMasses(){
		return this.linkMasses;
	}
	
	public double[][] getRc() {
		return rc;
	}

	public double[][] getRici() {
		return rici;
	}
	//END of GET methods

	
	//Part of SET methods (if needed)


	//END of SET methods

	//OTHER PUBLIC methods
	


	public Target Hto_from(int from, int to, double[] joint_values){

		/* It returns the homogeneous matrix (H) representing position and orientation of
		 * frame "from" with respect to frame "to" 
		 * 
		 * joint_values = values of the joint coordinates
		 */

		List<double[][]> H_indexed = new ArrayList<double[][]>();
		Target T;
		double[][] I = new double[4][4];
		double[][] H = new double[4][4];
		int max;
		int min;

		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++)
				if(i == j)
					I[i][j] = 1;
				else I[i][j] = 0;

		H = I;

		if(from == to){
			T = new Target(H);
			return T;
		}


		double[][] H0_1 = {
				{Math.cos(joint_values[0]), -Math.sin(joint_values[0]), 0, 0},
				{Math.sin(joint_values[0]), Math.cos(joint_values[0]), 0, 0},
				{0, 0, 1, this.linkLength[0][1]},
				{0, 0, 0, 1}
		};
		double[][] H1_2 = {
				{Math.cos(joint_values[1]), -Math.sin(joint_values[1]), 0, this.linkLength[1][0]},
				{0, 0, -1, 0},
				{Math.sin(joint_values[1]), Math.cos(joint_values[1]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H2_3 = {
				{Math.cos(joint_values[2]), -Math.sin(joint_values[2]), 0, 0},
				{Math.sin(joint_values[2]), Math.cos(joint_values[2]), 0, 0},
				{0, 0, 1, this.linkLength[2][1]},
				{0, 0, 0, 1}
		};
		double[][] H3_4 = {
				{0, 0, 1, this.linkLength[3][0] + this.linkLength[4][0]},
				{Math.sin(joint_values[3]), Math.cos(joint_values[3]), 0, 0},
				{-Math.cos(joint_values[3]), Math.sin(joint_values[3]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H4_5 = {
				{0, 0, -1, 0},
				{Math.sin(joint_values[4]), Math.cos(joint_values[4]), 0, 0},
				{Math.cos(joint_values[4]), -Math.sin(joint_values[4]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H5_6 = {
				{0, 0, 1, 0},
				{Math.sin(joint_values[5]), Math.cos(joint_values[5]), 0, 0},
				{-Math.cos(joint_values[5]), Math.sin(joint_values[5]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H6_endEff = {
				{1, 0, 0, 0},
				{0, 1, 0, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1}
		};

		H_indexed.add(H0_1);
		H_indexed.add(H1_2);
		H_indexed.add(H2_3);
		H_indexed.add(H3_4);
		H_indexed.add(H4_5);
		H_indexed.add(H5_6);
		H_indexed.add(H6_endEff);


		if(from > to){
			max = from;
			min = to;
		}else{
			max = to;
			min = from;
		}

		for(int i = max - 1; i >= min; i--)
			H = MultiplyMatrices(H_indexed.get(i), H);

		T = new Target(H);

		if(from < to){
			H = T.getInvHomMatrix();
			T.setHomMatrix(H);
		}

		return T;

	}
	


	//END of PUBLIC methods

	//
	//Part of PRIVATE methods that are used/called by inner methods
	//

	private void initializeArray(double[] array){

		for(int i = 0; i < array.length; i++)
			array[i] = 0;

	}
	
	private void initializeMatrix(double[][] m){

		for(int i = 0; i < m.length; i++)
			for(int j = 0; j < m[0].length; j++)
				m[i][j] = 0;

	}
	
	private double[][] MultiplyScalMatr(double scalar, double[][] m){

		for(int i = 0; i < m.length; i++)
			for(int j = 0; j < m[0].length; j++)
				m[i][j] = scalar*m[i][j];

		return m;

	}
	
	private double[][] MultiplyMatrices(double[][] left, double[][] right){

		/*
		 * Function implemented to multiply 2 matrices (useful in Hto_from)
		 */

		double[][] M = new double[left.length][right[0].length];

		for(int row = 0; row < M.length; row++)
			for(int column = 0; column < M[0].length; column++)
				M[row][column] = 0;

		if(left[0].length == right.length){

			for(int row = 0; row < M.length; row++)
				for(int column = 0; column < M[0].length; column++)
					for(int i = 0; i < left[0].length; i++)
						M[row][column] += left[row][i] * right[i][column];

		}else{
			System.out.println("The dimension of the matrices do not agree. [m*n n*p = m*p]");
		}

		return M;

	}
	
	private double[][] SubtractMatrices(double[][] m1, double[][] m2){

		/*
		 * It returns m = m1 - m2.
		 * If dimensions do not match, it returns null matrix of dim = m1.
		 */

		double[][] m = new double[m1.length][m1[0].length];

		if(m2.length == m1.length && m2[0].length == m1[0].length){

			for(int i = 0; i < m1.length; i++)
				for(int j = 0; j < m1[0].length; j++)
					m[i][j] = m1[i][j] - m2[i][j];

		}else{
			System.out.println("Matrices dimension do not match.");
			initializeMatrix(m);
		}

		return m;

	}
	
	private void copyRobotParameters (BufferedReader b, int section) throws IOException
	{
		String s;
		String comment_limiter = "%%";
		String section_limiter = "**--";
		
		double[][] temp_matr = new double[3][3];
		
		
		int row_number = 0; //it's used to remember which row of the matrix we are copying
		//Section where I still switch but I instantiate the attributes of the robot
		//to be filled with data from the file

		switch(section){
		case 2: 
			this.linkLength = new double[this.dof][2];
			break;
		case 3:
			this.linkMasses = new double[this.dof];
			break;
		case 4:
			this.jointLimits = new double[this.dof][2];
			break;
		case 5:
			this.cogMatrix = new double[3][this.dof];
			break;
		case 6:
			this.inertiaTens = new ArrayList<double[][]>();
			

		default:
			break;
		}



		while(true){
			s = b.readLine();
			if(s == null)
			{
				break;
			}
			if(s.startsWith(comment_limiter) || s.isEmpty())
			{
				continue;
			}
			if(s.startsWith(section_limiter))
			{
				break;
			}
			if(section == 0)
			{
				this.model = s;
				System.out.println(this.model);
				continue;
			}
			String[] tokens = s.split(" ");
			//System.out.println("the number of tokens in this string is " + tokens.length);
			System.out.println(" ");
			for(int i = 0; i<tokens.length; i++ )
			{
				
				System.out.print(tokens[i] + " ");
			}
			
			switch(section){
			case 1: 
				this.dof = Integer.parseInt(tokens[0]);
				break;
			case 2:
				this.linkLength[row_number][0] = Double.parseDouble(tokens[0]);
				this.linkLength[row_number][1] = Double.parseDouble(tokens[1]);
				row_number++;
				break;
			case 3:
				for(int i = 0; i<tokens.length; i++ )
				{
					this.linkMasses[i] = Double.parseDouble(tokens[i]);
				}
				break;
			case 4:
				//here the values are read in degrees but will be converted into RAD
				this.jointLimits[row_number][0] = Double.parseDouble(tokens[0]) * Math.PI / 180;
				this.jointLimits[row_number][1] = Double.parseDouble(tokens[1]) * Math.PI / 180;
				row_number++;
				break;
			case 5:
				/*
				 * REMINDER
				 * The center of gravity matrix is a 3xDOF matrix with ( x )
				 * 													   ( y )
				 * 													   ( z )
				 * in every column for each joint of the robot.
				 */
				for(int i = 0; i<tokens.length; i++ )
				{

					this.cogMatrix[row_number][i] = Double.parseDouble(tokens[i]);

				}
				row_number++;
				break;
			case 6:
				temp_matr[row_number][0] = Double.parseDouble(tokens[0]);
				temp_matr[row_number][1] = Double.parseDouble(tokens[1]);
				temp_matr[row_number][2] = Double.parseDouble(tokens[2]);
				row_number++;
				
				if(row_number == 3)
				{
					row_number = 0;
					this.inertiaTens.add(temp_matr);
				}
				
				break;
			case 7:
				break;

			default:
				break;
			}
		}
	}

	private void notifyActualSection(int section)
	{
		String sec;

		switch(section){
		case 0:
			sec = "Model of the Robot";
			break;
		case 1: 
			sec = "Degrees of Freedom";
			break;
		case 2: 
			sec = "Move to the wirst: joints lengths";
			break;
		case 3:
			sec = "Link masses ";
			break;
		case 4: 
			sec = "Joints working range: joints limits";
			break;
		case 5:
			sec = "Center of Gravity ";
			break;
		case 6:
			sec = "Inertia Tensors ";
			break;

		default:
			sec = " ";
			break;
		}
		System.out.println("The current section being compiled is: " + sec);
	}

	

	private void cog_wrtFrame0(double[][] cg){

		/*
		 * Input: cg(cog_Matrix) with respect to global frames (SolidWorks)
		 * Output: It gives you cog_Matrix (cg) values with respect to frame 0
		 */


		double[] joint_values = new double[this.dof];
		double[][] cg_i = new double[3][1];
		double[][] cg0 = new double[3][this.dof];
		double[][] cg0_i = new double[3][1];
		Target T0_i = new Target();
		double temp;

		/*
		 * I divide cg by 1000 in order to get m instead of mm
		 */
		for(int i = 0; i < cg.length; i++)
			for(int j = 0; j < cg[0].length; j++)
				cg[i][j] /= 1000;

		for(int i = 0; i < this.dof; i++){

			cg_i[0][0] = cg[0][i];
			cg_i[1][0] = cg[1][i];
			cg_i[2][0] = cg[2][i];

			initializeArray(joint_values);

			T0_i = Hto_from(i+1, 0, joint_values);
			cg0_i = MultiplyMatrices(T0_i.getInvRotation(), cg_i);

			cg0[0][i] = cg0_i[0][0];
			cg0[1][i] = cg0_i[1][0];
			cg0[2][i] = cg0_i[2][0];

		}

	}

	private List<double[][]> inTensi_wrtFramei(List<double[][]> inTens_sw){

		/*
		 * Input: Number of the joint of which you want the inertia tensors
		 * Output: Inertia tensors (3x3 matrix)
		 */

		int linkNum = 0;
		double prop = 98000/46249; /* ??????????????????? */
		double[] joint_values = new double[this.dof];
		double[] massesInMm = new double[this.linkMasses.length];
		double[][] BUFFER = new double[3][3];
		double[][] inTensi_i = new double[3][3];
		double[][] inTens0_i = new double[3][3];
		List<double[][]> inTens_prop = new ArrayList<double[][]>();
		List<double[][]> inTens = new ArrayList<double[][]>();
		Target T = new Target();
		
		/*
		 * from link_masses in m to link_masses in mm
		 */
		for(int i = 0; i < massesInMm.length; i++)
			massesInMm[i] = this.linkMasses[i]*1000;

		cog_wrtFrame0(this.cogMatrix); 
		
		/*
		 * Values from cog_wrtFrame0 are in m so I have to multiply them by 1000
		 */
		for(int i = 0; i < cogMatrix.length; i++)
			for(int j = 0; j < cogMatrix[0].length; j++)
				cogMatrix[i][j] *= 1000;
		
		inTens_prop = inTens_sw;
		/* I multiply inTens_00 * prop */
		for(double[][] link:inTens_prop)
			for(int i = 0; i < 3; i++)
				for(int j = 0; j < 3; j++)
					link[i][j] *= prop;
		
		for(double[][] link:inTens_sw){

			link[0][0] += massesInMm[linkNum]*(Math.pow(this.cogMatrix[1][linkNum], 2) + Math.pow(this.cogMatrix[2][linkNum], 2));
			link[1][1] += massesInMm[linkNum]*(Math.pow(this.cogMatrix[0][linkNum], 2) + Math.pow(this.cogMatrix[2][linkNum], 2));
			link[2][2] += massesInMm[linkNum]*(Math.pow(this.cogMatrix[0][linkNum], 2) + Math.pow(this.cogMatrix[1][linkNum], 2));

			link[0][1] -= massesInMm[linkNum]*(this.cogMatrix[0][linkNum]*this.cogMatrix[1][linkNum]);
			link[1][0] = link[0][1];
			link[0][2] -= massesInMm[linkNum]*(this.cogMatrix[0][linkNum]*this.cogMatrix[2][linkNum]);
			link[2][0] = link[0][2];
			link[1][2] -= massesInMm[linkNum]*(this.cogMatrix[1][linkNum]*this.cogMatrix[2][linkNum]);
			link[2][1] = link[1][2];

			inTens0_i = SubtractMatrices(inTens_prop.get(linkNum), link);

			initializeArray(joint_values);

			T = Hto_from(linkNum, 0, joint_values);
			/* I have to perform the following calculation: Ii_i = R0_i^(-1)*I0_i*R0_i
			 * to do so, I use a buffer 
			 */
			BUFFER = MultiplyMatrices(T.getInvRotation(), inTens0_i);
			inTensi_i = MultiplyMatrices(BUFFER, T.getRotation());

			/* I need inTens to be in kg*m^2(so far they were in g*mm^2)
			 *  then I multiply every element of the list for 10^-9 */			
			inTensi_i = Matrix.multiplyScalMatr(10E-9, inTensi_i);
			inTens.add(inTensi_i);

			linkNum++;

		}

		return inTens;

	}

	private void xmlBuilder(String FilePath){
		
		try{
			 File inputFile = new File(FilePath);
	         SAXReader reader = new SAXReader();
	         Document document = reader.read( inputFile );
	         
	         
	         String rootName = document.getRootElement().getName();
	         int DoF =(int) Double.parseDouble(document.selectSingleNode(rootName + "/DOF").getText());
	         
	         initializeVariables(DoF); //Initialization of all Robot variables according to the number of DOF
	         
	         List<Node> links = document.selectNodes(rootName + "/link");
	         
	         for(Node link : links)
	         {
	        	 /*
	     		 * In this method the link in input is read according to the XML format and
	     		 * its data are retrieved and used to fill Robot's matrices.
	     		 * 
	     		 * The variable initialization might be verbose but makes clearer the reading of the
	     		 * XML file as he proceeds. All the variables are first instanced and then copied into
	     		 * the Robot's private variables.
	     		 */
	        	 
	        	 int linkNumber = (int) Double.parseDouble(link.selectSingleNode("number").getText()) - 1;
	        	 Node node_length = document.selectSingleNode(rootName + "/link/length");
	        	 double a = Double.parseDouble(link.selectSingleNode("length/a").getText());
	        	 double d = Double.parseDouble(link.selectSingleNode("length/d").getText());
	        	 double mass = Double.parseDouble(link.selectSingleNode("mass").getText());
	        	 double max_limit = Double.parseDouble(link.selectSingleNode("limit/max").getText());
	        	 double min_limit = Double.parseDouble(link.selectSingleNode("limit/min").getText());
	        	 double[] CoG = new double[3];
	        	 CoG[0] = Double.parseDouble(link.selectSingleNode("CoG/x").getText());
	        	 CoG[1] = Double.parseDouble(link.selectSingleNode("CoG/y").getText());
	        	 CoG[2] = Double.parseDouble(link.selectSingleNode("CoG/z").getText());
	        	 
	        	 double[][] inertia_link = new double[3][3];
	        	 inertia_link[0][0] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_11").getText());
	        	 inertia_link[1][1] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_22").getText());
	        	 inertia_link[2][2] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_33").getText());
	        	 inertia_link[0][1] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_12").getText());
	        	 inertia_link[0][2] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_13").getText());
	        	 inertia_link[1][2] = Double.parseDouble(link.selectSingleNode("InertiaTensor/elem_23").getText());
	        	 inertia_link[2][1] = inertia_link[1][2]; //the InertiaTensor matrix is symmetric
	        	 inertia_link[2][0] = inertia_link[0][2];
	        	 inertia_link[1][0] = inertia_link[0][1];
	        	 
	        	 this.linkLength[linkNumber][0] = a;
	        	 this.linkLength[linkNumber][1] = d;
	        	 this.linkMasses[linkNumber] = mass;
	        	 this.jointLimits[linkNumber][0] = min_limit * Math.PI / 180;
	        	 this.jointLimits[linkNumber][1] = max_limit * Math.PI / 180;
	        	 
	        	 for(int i=0; i < CoG.length; i++){this.cogMatrix[i][linkNumber] = CoG[i];}
	        	 
	        	 this.inertiaTens.add(linkNumber, inertia_link);
	        	 
	         }
	         
		}
		catch(DocumentException e){
			System.out.println("Check that the FilePath name is correct please. Unable to locate the file ");
			e.printStackTrace();
		}
	}
	
	private void jointDistances(){
		double[] theta = new double[this.dof];
		Arrays.fill(theta, 0);
		
		for(int i=0; i < this.dof; i++){
			Target Tnum = this.Hto_from(i, 0, theta);
			double[][] temp = Tnum.getInvHomMatrix();
			double[] cogLink = {this.cogMatrix[0][i],this.cogMatrix[1][i],this.cogMatrix[2][i],1}; //Vector to be considered 4x1
			/*
			 * The whole calculation could be optimized but is not really worth it
			 */
			double[] rc4Link = Matrix.multiplyMatrixVector(temp, cogLink);
			double[] rcLink = {rc4Link[0],rc4Link[1],rc4Link[2]};
			
			this.rc[i] = rcLink;
		}
		
		//r=[a1,0,0;0,d2,0;(a2+a3),0,0;0,0,0;0,0,0; 0,0,0] in MATLAB notation
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
		
		double[][] r = {
				{linkLength[1][0],0,0},
				{0,linkLength[2][1],0},
				{linkLength[3][0]+linkLength[4][0],0,0},
				{0,0,0},
				{0,0,0},
				{0,0,0}
		};
		this.rici = Matrix.subtractMatrices(rc, r);
		
	}
	
	private void initializeVariables(int DOF){
		
		this.dof = DOF;
		this.linkLength = new double[DOF][2];
		this.linkMasses = new double[DOF];
		this.jointLimits = new double[DOF][2];
		this.cogMatrix = new double[3][DOF];
		this.inertiaTens = new ArrayList<double[][]>();
		this.rc = new double[DOF][3];
		this.rici = new double[DOF][3];
		
	}
	
	
	// END of PRIVATE methods


}
