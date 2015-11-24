
public class Test_robot_constructor {

	public static void main(String[] args) {
		String path = "C:\\Users\\Enrico\\Google Drive\\Tesi\\Prova Costruttore Robot\\ABB IRB140.txt"; //percorso del file con i dati
		//String path = "C:\\Users\\Enrico\\Desktop\\MatriceProva.txt";
		
		Robot newRobot = new Robot(path);
		int DOF = newRobot.getDOF();
		System.out.println("The degree of freedom of the robot are: " + DOF);

	}

}
