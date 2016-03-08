import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

public class Target {

	private double[][] homMatrix;
	private double[] position;
	private double[][] rotation;


	public static void main(String[] args){

	}

	//Part of CONSTRUCTORS

	public Target(){

		homMatrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		double[][] I = new double[4][4];
		I = identityMatr();

		for(int i = 0; i < I.length; i++)
			for (int j = 0; j < I[i].length; j++){
				homMatrix[i][j] = I[i][j];
				if(i != 3 && j == 3)
					position[i] = I[i][j];
				if(i != 3 && j != 3)
					rotation[i][j] = I[i][j];
			}
	}

	public Target(double[][] T){

		homMatrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		if(T.length == 4 && T[0].length == 4){		
			for(int i = 0; i < T.length; i++)
				for (int j = 0; j < T[0].length; j++){
					homMatrix[i][j] = T[i][j];
					if(i != 3 && j == 3)
						position[i] = T[i][j];
					if(i != 3 && j != 3)
						rotation[i][j] = T[i][j];
				}
		}else{
			System.out.println("Wrong dimensions. Homogeneus matrix must be 4x4.");
		}

	}

	public Target(double[] pos, double[][] rot){

		homMatrix = new double[4][4];
		position = new double[3];
		rotation = new double[3][3];

		homMatrix = identityMatr();

		if(pos.length == 3){
			for(int i = 0; i < pos.length; i++){
				homMatrix[i][3] = pos[i];
				position[i] = pos[i];
			}
		}else{
			System.out.println("Wrong dimensions. Position vector must be a vector of length = 3.");
			return;
		}

		if(rot.length == 3 && rot[0].length == 3){
			for(int i = 0; i < rot.length; i++)
				for(int j = 0; j < rot[0].length; j++){
					homMatrix[i][j] = rot[i][j];
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

		return homMatrix;

	}
	
public double[][] getInvHomMatrix(){
		
		/*
		 * Inverse of an Affine Matrix
		 */
		
		double[][] rotMatrix = new double[3][3];
		Matrix.fill(rotMatrix, this.rotation);
		
		rotMatrix = Matrix.transpose(rotMatrix);
		
		double[] posVector = new double[3];
		for(int i = 0; i < posVector.length; i++){
			posVector[i] = - position[i];
		}
		posVector = Matrix.multiplyMatrixVector(rotMatrix, posVector);

		Target result = new Target(posVector,rotMatrix);

		return result.getHomMatrix();

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

		double[][] INV = new double[3][3];

		for(int i = 0; i < this.getRotation().length; i++)
			for(int j = 0; j < this.getRotation()[0].length; j++)
				INV[i][j] = this.getRotation()[j][i];

		return INV;

	}

	//END of GET methods

	//SET methods
	
	public void setHomMatrix(double[][] H){
		
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++)
				this.homMatrix[i][j] = H[i][j];
		
	}
	
	//END of SET methods
	
	
	//OTHER PUBLIC methods

	public void translateTarget (double[] vector)
	{
		
		BigDecimal[] buffer1 = new BigDecimal[vector.length];
		BigDecimal[] buffer2 = new BigDecimal[vector.length];
		BigDecimal[] res = new BigDecimal[vector.length];
		MathContext mc1 = new MathContext(2, RoundingMode.HALF_DOWN);
		MathContext mc2 = new MathContext(6, RoundingMode.HALF_DOWN);
		
//		for(int i = 0 ; i < this.position.length ; i++)
//		{
//			
//		this.position[i] += vector[i];
//		this.homMatrix[i][3] += vector[i];
//		
//		}
		
		for(int i = 0 ; i < this.position.length ; i++)
		{
						
					buffer1[i] = new BigDecimal(vector[i], mc1);				
					buffer2[i] = new BigDecimal(this.position[i], mc2);
					
					res[i] = buffer2[i].add(buffer1[i], mc2);
						
						this.position[i] = res[i].doubleValue();
						this.homMatrix[i][3] = res[i].doubleValue();
						
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
