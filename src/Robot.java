import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Collections;

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
	
	public double[][] getCoGMatrix(){
		return this.cogMatrix;
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
	


	public Target hToFrom(int from, int to, double[] joint_values){

		/* It returns the homogeneous matrix (H) representing position and orientation of
		 * frame "from" with respect to frame "to" 
		 * 
		 * joint_values = values of the joint coordinates
		 */

		List<double[][]> hIndexed = new ArrayList<double[][]>();
		Target T;
		double[][] I = new double[4][4];
		double[][] H = new double[4][4];
		int max;
		int min;

		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				if(i == j){
					I[i][j] = 1;
				}else {
					I[i][j] = 0;
				}
				H[i][j] = I[i][j];
			}
		}

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
				{Math.sin(joint_values[2]), Math.cos(joint_values[2]), 0, this.linkLength[2][1]},
				{0, 0, 1, 0},
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

		hIndexed.add(H0_1); //index 0 in the list
		hIndexed.add(H1_2); //index 1
		hIndexed.add(H2_3); //index 2 etc.
		hIndexed.add(H3_4);
		hIndexed.add(H4_5);
		hIndexed.add(H5_6);
		hIndexed.add(H6_endEff);


		if(from > to){
			max = from;
			min = to;
		}else{
			max = to;
			min = from;
		}

		for(int i = max - 1; i >= min; i--){
			H = Matrix.multiplyMatrices(hIndexed.get(i), H); //original line
//			H = Matrix.multiplyMatrices(H, hIndexed.get(i));
		}

		T = new Target(H);

		if(from < to){
			H = T.getInvHomMatrix();
//			T.setHomMatrix(H);
			T = new Target(H);
		}

		return T;

	}
	


	//END of PUBLIC methods

	//
	//Part of PRIVATE methods that are used/called by inner methods
	//

	
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

	private void inTensiWRTFramei(){

		/*
		 * Input: Number of the joint of which you want the inertia tensors
		 * Output: Inertia tensors (3x3 matrix)
		 */

		int linkNum = 0;
		double prop = 98000.0/46249.0 ; /* ??????????????????? */
		
//		System.out.println("Prop: " + prop);
		
		double[] jointValues = new double[this.dof];
		double[] massesInGr = new double[this.linkMasses.length];
		double[][] BUFFER = new double[3][3];
		double[][] inTensi = new double[3][3];
		double[][] inTens0 = new double[3][3];
		Target T = new Target();
		
		/*
		 * from link_masses in m to link_masses in mm
		 */
		for(int i = 0; i < massesInGr.length; i++){
			massesInGr[i] = this.linkMasses[i]*1000;
		}
		
		//cog_wrtFrame0(this.cogMatrix); 
	
		
	
		
		for(double[][] link:this.inertiaTens){
			
			double[][] linkProp = new double[link.length][link[0].length];
			for(int i = 0; i < linkProp.length; i++){
				for(int j = 0; j < linkProp[0].length; j++){
					linkProp[i][j] = link[i][j]*prop;
				}
			}
//			
//			System.out.println("Il tensore d'inerzia "+ (linkNum+1) + "(originale) vale :");
//			Matrix.displayMatrix(link);
//			System.out.println(" ");
			
			link[0][0] += massesInGr[linkNum]*(Math.pow(this.cogMatrix[1][linkNum], 2) + Math.pow(this.cogMatrix[2][linkNum], 2));			
			link[1][1] += massesInGr[linkNum]*(Math.pow(this.cogMatrix[0][linkNum], 2) + Math.pow(this.cogMatrix[2][linkNum], 2));
			link[2][2] += massesInGr[linkNum]*(Math.pow(this.cogMatrix[0][linkNum], 2) + Math.pow(this.cogMatrix[1][linkNum], 2));

			link[0][1] -= massesInGr[linkNum]*(this.cogMatrix[0][linkNum]*this.cogMatrix[1][linkNum]);
			link[1][0] = link[0][1];
			link[0][2] -= massesInGr[linkNum]*(this.cogMatrix[0][linkNum]*this.cogMatrix[2][linkNum]);
			link[2][0] = link[0][2];
			link[1][2] -= massesInGr[linkNum]*(this.cogMatrix[1][linkNum]*this.cogMatrix[2][linkNum]);
			link[2][1] = link[1][2];

//			System.out.println("Mass " + (linkNum+1)+ ": " + massesInGr[linkNum]);
//			System.out.println("CoGMatrix [0]["+linkNum+"]  : " + this.cogMatrix[0][linkNum]);
//			System.out.println("CoGMatrix [1]["+linkNum+"]  : " + this.cogMatrix[1][linkNum]);
//			System.out.println("CoGMatrix [2]["+linkNum+"]  : " + this.cogMatrix[2][linkNum]);
//			System.out.println("Il tensore d'inerzia "+ (linkNum+1) + " (modificato) vale :");
//			Matrix.displayMatrix(link);
//			System.out.println(" ");
//			System.out.println("Mentre quello motliplicato per prop vale: ");
//			Matrix.displayMatrix(linkProp);
//			System.out.println(" ");
//			
			
			inTens0 = Matrix.subtractMatrices(linkProp, link);
			
//			
//			System.out.println("Igc"+(linkNum+1)+" :");
//			Matrix.displayMatrix(inTens0);
//			System.out.println(" ");
			
			Arrays.fill(jointValues, 0);

			T = hToFrom(linkNum+1, 0, jointValues);
			/* I have to perform the following calculation: Ii_i = R0_i^(-1)*I0_i*R0_i
			 * to do so, I use a buffer 
			 */
			BUFFER = Matrix.multiplyMatrices(T.getInvRotation(), inTens0);
			inTensi = Matrix.multiplyMatrices(BUFFER, T.getRotation());

//			System.out.println("R0_" + (linkNum+1) +" vale :");
//			Matrix.displayMatrix(T.getRotation());
//			System.out.println(" ");
//			System.out.println("R0_"+ (linkNum+1) +"^-1 vale :");
//			Matrix.displayMatrix(T.getInvRotation());
//			System.out.println(" ");
			
			/* I need inTens to be in kg*m^2(so far they were in g*mm^2)
			 *  then I multiply every element of the list by 10^-9 */			
			inTensi = Matrix.multiplyScalMatr(1E-9, inTensi); //it seems like it should be E-10....why??
			
//			
//			System.out.println(" ");
//			System.out.println("Tensore d'inerzia i_i modificato di indice  " + (linkNum+1));
//			Matrix.displayMatrix(inTensi);
//			System.out.println(" ");
			
			link = inTensi;
			for(int i = 0; i < link.length; i++){
				for(int j = 0; j < link[0].length; j++){
					link[i][j] = inTensi[i][j];
				}
			}
			
			double[][] TEST = new double[3][3];
			Matrix.fill(TEST, 0);
			link = TEST;
			
//			System.out.println(" ");
//			System.out.println("Tensore d'inerzia modificato di indice  " + linkNum);
//			Matrix.displayMatrix(inTensi_i);
//			System.out.println(" ");

			linkNum++;

		}
		
		
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
	         
	        this.inTensiWRTFramei();
	        this.jointDistances();
//	        System.out.println(" ");
//			Matrix.displayMatrix(this.inertiaTens.get(0));
//			System.out.println(" ");
	        
			System.out.println("Rc matrix:");
			Matrix.displayMatrix(this.rc);
			System.out.println(" ");
			
			System.out.println("Rici matrix: ");
			Matrix.displayMatrix(this.rici);
			System.out.println(" ");
	        
		}
		catch(DocumentException e){
			System.out.println("Check that the FilePath name is correct please. Unable to locate the file ");
			e.printStackTrace();
		}
	}
	
	private void jointDistances(){
		
		
		/*
		 * Change the CoG matrix unit from mm to m
		 */
		
		for(int i = 0; i < cogMatrix.length; i++)
			for(int j = 0; j < cogMatrix[0].length; j++)
				cogMatrix[i][j] /= 1000;
		
		
		double[] theta = new double[this.dof];
		Arrays.fill(theta, 0);
		
		for(int i=0; i < this.dof; i++){
			Target Tnum = this.hToFrom(i+1, 0, theta);
//			
//			System.out.println("Position Vector. Index: " + (i+1));
//			Matrix.displayVector(Tnum.getPosition());
//			System.out.println(" ");
//			
//			System.out.println("Rotation matrix. Index: " + (i+1));
//			Matrix.displayMatrix(Tnum.getRotation());
//			System.out.println(" ");
			
			System.out.println("Homogenous matrix. Index: " + (i+1));
			Matrix.displayMatrix(Tnum.getHomMatrix());
			System.out.println(" ");
			
			
			double[][] temp = Tnum.getInvHomMatrix();
			
//			System.out.println("Inverse homogenous matrix. Index: " + (i+1));
//			Matrix.displayMatrix(temp);
//			System.out.println(" ");
			
			double[] cogLink = {this.cogMatrix[0][i],this.cogMatrix[1][i],this.cogMatrix[2][i],1}; //Vector to be considered 4x1
			/*
			 * The whole calculation could be optimized but is not really worth it
			 */
			double[] rc4Link = Matrix.multiplyMatrixVector(temp, cogLink);
			double[] rcLink = {rc4Link[0],rc4Link[1],rc4Link[2]};
			
			System.out.println("Rc4Link. Index: " + (i+1));
			Matrix.displayVector(rc4Link);
			System.out.println(" ");
			
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
