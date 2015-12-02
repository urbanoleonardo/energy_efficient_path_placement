
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
				{1.0, 0.0, 0.0, 2.3},
				{0.0, 1.0, 0.0, -3.5},
				{0.0, 0.0, 1.0, 1.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		double[][] matrix2 = {
				{1.0, 0.0, 0.0, 4.6},
				{0.0, 1.0, 0.0, 4.0},
				{0.0, 0.0, 1.0, -1.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double x_sample = 0.1;
		double t_sample = 0.1;
		double acceleration = 2.0;
		double velocity = 6.0;
		
		Target test_target1 = new Target(matrix1);
		Target test_target2 = new Target(matrix2);
		Path test_path = new Path(test_target1, test_target2, x_sample, t_sample);
		
		test_path.SetMaxAcc(acceleration);
		test_path.SetMaxVel(velocity);
		test_path.Interpolate();
		Trajectory test_output = test_path.GetTrajectory();
		System.out.println("Numero di punti : " + test_output.Points.size());
		
		
		for(Target i : test_output.Points)
		{
			double[] position_vector = i.getPosition();
			
			System.out.println("Punto : " + position_vector[0] + " " + position_vector[1] + " " + position_vector[2]);
		}
		
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
