import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Math;
import java.util.ArrayList;
import java.util.List;

/*
 * TODO it defines all the parameters of the robot
 * 
 */

public class Robot {

	private String model;
	private int dof;
	private double[][] link_length;
	private double[][] link_masses;
	private double[][] inertia_matr;
	private double[][] inertia_tens;
	private double[][] joint_limits;
	private double[][] cog_Matrix; //center of gravity matrix (SolidWorks)
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

						NotifyActualSection(section);
						System.out.println(s);
						CopyRobotParameters(b, section);
						section++;
					}

				}

			}
		}
		catch (IOException e){
			System.out.println("Error in opening the file");
		}

	}

	//END part of CONSTRUCTORS




	//Part of GET methods



	public int getDOF(){
		return this.dof;
	}

	public double[][] getParameters(String arg)
	{

		/*
		 * Function to get the parameters of the robot.
		 * Consider using only some key words no matter the overall string given in input maybe
		 * I don't think it would be better having a single method for every different parameter
		 * just because they are all of the same type;
		 */

		double[][] defaultRes = new double[0][0];

		if(arg.equalsIgnoreCase("link length") || arg.equalsIgnoreCase("link_length"))
		{
			return this.link_length;
		}
		if(arg.equalsIgnoreCase("link masses") || arg.equalsIgnoreCase("link_masses"))
		{
			return this.link_masses;
		}
		if(arg.equalsIgnoreCase("joint limits") || arg.equalsIgnoreCase("joint_limits"))
		{
			return this.joint_limits;
		}
		if(arg.equalsIgnoreCase("center of gravity"))
		{
			return this.cog_Matrix;
		}

		return defaultRes; 
	}

	//END of GET methods




	//Part of SET methods (if needed)


	//END of SET methods

	//OTHER PUBLIC methods



	//END of PUBLIC methods

	//
	//Part of PRIVATE methods that are used/called by inner methods
	//


	private void CopyRobotParameters (BufferedReader b, int section) throws IOException
	{
		String s;
		String comment_limiter = "%%";
		String section_limiter = "**--";

		int row_number = 0; //it's used to remember which row of the matrix we are copying
		//Section where I still switch but I instantiate the attributes of the robot
		//to be filled with data from the file

		switch(section){
		case 2: 
			this.link_length = new double[this.dof][2];
			break;
		case 3:
			this.link_masses = new double[this.dof][1];
			break;
		case 4:
			this.joint_limits = new double[this.dof][2];
			break;
		case 5:
			this.cog_Matrix = new double[3][this.dof];
			break;

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
			for(int i = 0; i<tokens.length; i++ )
			{
				System.out.println(tokens[i]);
			}
			switch(section){
			case 1: 
				this.dof = Integer.parseInt(tokens[0]);
				break;
			case 2:
				this.link_length[row_number][0] = Double.parseDouble(tokens[0]);
				this.link_length[row_number][1] = Double.parseDouble(tokens[1]);
				row_number++;
				break;
			case 3:
				for(int i = 0; i<tokens.length; i++ )
				{
					this.link_masses[i][0] = Double.parseDouble(tokens[i]);
				}
				break;
			case 4:
				//here the values are read in degrees but will be converted into RAD
				this.joint_limits[row_number][0] = Double.parseDouble(tokens[0]) * Math.PI / 180;
				this.joint_limits[row_number][1] = Double.parseDouble(tokens[1]) * Math.PI / 180;
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

					this.cog_Matrix[row_number][i] = Double.parseDouble(tokens[i]);

				}
				row_number++;
				break;
			case 6:
				break;
			case 7:
				break;

			default:
				break;
			}
		}
	}

	private void NotifyActualSection(int section)
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

		default:
			sec = " ";
			break;
		}
		System.out.println("The current section being compiled is: " + sec);
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

	private Target Hto_from(int from, int to, double[] joint_values){

		/* It returns the homogeneous matrix (H) representing position and orientation of
		 * frame "from" with respect to frame "to" 
		 * 
		 * joint_values = values of the joint coordinates
		 */

		List<double[][]> H_indexed = new ArrayList<double[][]>(this.dof);
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
				{Math.cos(joint_values[1]), -Math.sin(joint_values[1]), 0, 0},
				{Math.sin(joint_values[1]), Math.cos(joint_values[1]), 0, 0},
				{0, 0, 1, this.link_length[1][2]},
				{0, 0, 0, 1}
		};
		double[][] H1_2 = {
				{Math.cos(joint_values[2]), -Math.sin(joint_values[2]), 0, this.link_length[2][1]},
				{0, 0, -1, 0},
				{Math.sin(joint_values[2]), Math.cos(joint_values[2]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H2_3 = {
				{Math.cos(joint_values[3]), -Math.sin(joint_values[3]), 0, 0},
				{Math.sin(joint_values[3]), Math.cos(joint_values[3]), 0, 0},
				{0, 0, 1, this.link_length[3][2]},
				{0, 0, 0, 1}
		};
		double[][] H3_4 = {
				{0, 0, 1, this.link_length[4][1] + this.link_length[5][1]},
				{Math.sin(joint_values[4]), Math.cos(joint_values[4]), 0, 0},
				{-Math.cos(joint_values[4]), Math.sin(joint_values[4]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H4_5 = {
				{0, 0, -1, 0},
				{Math.sin(joint_values[5]), Math.cos(joint_values[5]), 0, 0},
				{Math.cos(joint_values[5]), -Math.sin(joint_values[5]), 0, 0},
				{0, 0, 0, 1}
		};
		double[][] H5_6 = {
				{0, 0, 1, 0},
				{Math.sin(joint_values[6]), Math.cos(joint_values[6]), 0, 0},
				{-Math.cos(joint_values[6]), Math.sin(joint_values[6]), 0, 0},
				{0, 0, 0, 1}
		};

		H_indexed.add(H0_1);
		H_indexed.add(H1_2);
		H_indexed.add(H2_3);
		H_indexed.add(H3_4);
		H_indexed.add(H4_5);
		H_indexed.add(H5_6);


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

	private void initializeArray(double[] array){

		for(int i = 0; i < array.length; i++)
			array[i] = 0;

	}

	private void joint_distances(double[][] cg){

		/*
		 * Input: cg(cog_Matrix) with respect to global frames (SolidWorks)
		 * Output: cg coordinates with respect to frame 0
		 */


		double[] joint_values = new double[this.dof];
		double[][] cg_i = new double[3][1];
		double[][] cg0 = new double[3][this.dof];
		double[][] cg0_i = new double[3][1];
		Target T0_i = new Target();

		/*
		 * I divide cg by 1000 in order to get mm instead of m
		 */
		for(int i = 0; i < cg.length; i++)
			for(int j = 0; j < cg[0].length; j++)
				cg[i][j] /= 1000;

		for(int i = 0; i < this.dof; i++){

			cg_i[0][0] = cg[0][i];
			cg_i[1][0] = cg[1][i];
			cg_i[2][0] = cg[2][i];

			initializeArray(joint_values);

			/* From i to 0 because I want the inverse matrix */
			T0_i = Hto_from(0, i, joint_values);
			cg0_i = MultiplyMatrices(T0_i.getInvRotation(), cg_i);

			cg0[0][i] = cg0_i[0][0];
			cg0[1][i] = cg0_i[1][0];
			cg0[2][i] = cg0_i[2][0];

		}

	}

	// END of PRIVATE methods


}
