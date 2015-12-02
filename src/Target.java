

public class Target {

	private double[][] hom_matrix;
	private double[] position;
	private double[][] rotation;


	public static void main(String[] args){

	}

	//Part of CONSTRUCTORS

	public Target(){

		hom_matrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		double[][] I = new double[4][4];
		I = identityMatr();

		for(int i = 0; i < I.length; i++)
			for (int j = 0; j < I[i].length; j++){
				hom_matrix[i][j] = I[i][j];
				if(i != 3 && j == 3)
					position[i] = I[i][j];
				if(i != 3 && j != 3)
					rotation[i][j] = I[i][j];
			}

		return;

	}

	public Target(double[][] T){

		hom_matrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		if(T.length == 4 && T[0].length == 4){		
			for(int i = 0; i < T.length; i++)
				for (int j = 0; j < T[i].length; j++){
					hom_matrix[i][j] = T[i][j];
					if(i != 3 && j == 3)
						position[i] = T[i][j];
					if(i != 3 && j != 3)
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

		hom_matrix = identityMatr();

		if(pos.length == 3){
			for(int i = 0; i < pos.length; i++){
				hom_matrix[i][3] = pos[i];
				position[i] = pos[i];
			}
		}else{
			System.out.println("Wrong dimensions. Position vector must be a vector of length = 3.");
			return;
		}

		if(rot.length == 3 && rot[0].length == 3){
			for(int i = 0; i < rot.length; i++)
				for(int j = 0; j < rot[0].length; j++){
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
	
public double[][] getInvHomMatrix(){
		
		/*
		 * It returns the inverse matrix using Gaussian Elimination Algorithm
		 * Input: 4x4 matrix
		 */

		double temp;
		double[][] INV = new double[4][4];
		INV = identityMatr();

		for(int k = 0; k < this.getHomMatrix().length; k++){

			temp = this.getHomMatrix()[k][k];	

			for(int j = 0; j < this.getHomMatrix()[0].length; j++){

				this.getHomMatrix()[k][j] /= temp;									
				INV[k][j] /= temp;		

			}

			for(int i=0; i < this.getHomMatrix().length; i++){

				temp = this.getHomMatrix()[i][k];	

				for(int j = 0; j < this.getHomMatrix()[0].length; j++){

					if(i == k)
						break;

					this.getHomMatrix()[i][j] -= this.getHomMatrix()[k][j]*temp;						
					INV[i][j] -= INV[k][j]*temp;

				}
			}
		}

		return INV;

}

	public double[] getPosition(){

		return position;

	}

	public double[][] getRotation(){

		return rotation;

	}
	
	public double[][] getInvRotation(){
		
		/*
		 * R is orthogonal so its inverse matrix is equal to its transpose matrix
		 */

		double[][] INV = new double[4][4];

		for(int i = 0; i < this.getHomMatrix().length; i++)
			for(int j = 0; j < this.getHomMatrix()[0].length; j++)
				INV[i][j] = this.getHomMatrix()[j][i];

		return INV;

	}

	//END of GET methods

	//SET methods
	
	
	//END of SET methods
	
	
	//OTHER PUBLIC methods

	public void TranslateTarget (double[] vector)
	{
		for(int i = 0 ; i < this.position.length ; i++)
		{
		this.position[i] += vector[i];
		this.hom_matrix[i][3] = vector[i];
		}
	}
	
	
	
	
	//END of PUBLIC methods


	//PRIVATE methods

	private double[][] identityMatr(){

		double[][] I = new double[4][4];

		for(int i = 0; i < I.length; i++)
			for(int j = 0; j < I[0].length; j++)
				if(i == j)
					I[i][j] = 1;
				else I[i][j] = 0;

		return I;

	}

	//END of PRIVATE methods

}
