

public class Target {

	private double[][] hom_matrix;
	private double[] position;
	private double[][] rotation;


	public static void main(String[] args){

	}
	
	
	
	
	//Part of CONSTRUCTORS

	public Target(double[][] T){

		hom_matrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		if(T.length == 4 || T[0].length == 4){		
			for(int i = 0; i < T.length; i++)
				for (int j = 0; j < T[i].length; j++){
					hom_matrix[i][j] = T[i][j];
					if(i != 4 && j == 4)
						position[i] = T[i][j];
					if(i != 4 && j != 4)
						rotation[i][j] = T[i][j];
				}
		}else{
			System.out.println("Wrong dimensions. Homogeneus matrix must be 4x4.");
			return;
		}

	}

	public Target(double[] pos, double[][] rot){

		hom_matrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		hom_matrix[3][3] = 1;

		if(rotation.length == 3){
			for(int i = 0; i < pos.length;i++){
				hom_matrix[i][3] = pos[i];
				position[i] = pos[i];
			}
		}else{
			System.out.println("Wrong dimensions. Position vector must be a vector of length = 3.");
			return;
		}

		if(rot.length == 3 || rot[0].length == 3){
			for(int i = 0; i < 3; i++)
				for(int j = 0; j < 3; j++){
					hom_matrix[i][j] = rot[i][j];
					rotation[i][j] = rot[i][j];
				}
		}else{
			System.out.println("Wrong dimensions. Rotation matrix must be 3x3.");
			return;
		}
	}

	//END of CONSTRUCTORS
	
	
	//Part of GET methods
	
	public double[][] getHomMatrix(){

		return hom_matrix;

	}

	public double[] getPosition(){

		return position;

	}

	public double[][] getRotation(){

		return rotation;

	}
	
	//END of GET methods
	
	//OTHER PUBLIC methods
	
	
	
	
	
	//END of PUBLIC methods
	
	
	
	//PRIVATE methods
	
	
	

}
