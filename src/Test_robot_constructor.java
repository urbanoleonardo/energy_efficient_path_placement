
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
		for(int i=0; i<jointLimit.length; i++)
		{
			for(int j=0;j<jointLimit[0].length;j++)
			{
				System.out.print(jointLimit[i][j] + " ");
			}
			System.out.println("");
		}
	}

}
