import java.lang.Math;
import java.util.ArrayList;
import java.util.List;

public class Matrix {

	public static double[][] multiplyMatrices(double[][] left, double[][] right){

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
	
	public static double multiplyMatrices(double[] left, double[][] right){

		/*
		 * Function implemented to multiply 2 matrices (useful in Hto_from)
		 */

		double M = 0;


		if(left.length == right.length){

			for(int i = 0; i < left.length; i++){
						M += left[i] * right[0][i];
			}
		}else{
			System.out.println("The dimension of the matrices do not agree. [m*n n*p = m*p]");
		}

		return M;

	}
	
	public static double[] multiplyMatrixVector(double[][] matrix, double[] vector)
	{
		double[] M_vector = new double[matrix.length];
		
		
			for(int k = 0; k < matrix.length; k++){ M_vector[k] = 0;}
			for(int l = 0; l < M_vector.length; l++)
					for(int i = 0; i < matrix[0].length; i++)
						M_vector[l] += matrix[l][i] * vector[i];
			
			return M_vector;
		
	}
	
	public static double[][] multiplyMatrices(double[][] left, double[] right){

		

		double[][] M = new double[left.length][right.length];
		

		for(int row = 0; row < M.length; row++)
			for(int column = 0; column < M[0].length; column++)
				M[row][column] = 0;

		
		
		if(left[0].length == 1){

			for(int row = 0; row < M.length; row++)
				for(int column = 0; column < M[0].length; column++)
					for(int i = 0; i < left[0].length; i++)
						M[row][column] += left[row][i] * right[i];
			
			return M;
		}else{
			System.out.println("The dimension of the matrices do not agree. [m*n n*p = m*p]");
		}

		return M;
	}
	
	public static double[][] transpose(double[][] matrix){
		double[][] matrix_transpose = new double[matrix[0].length][matrix.length];
		
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[0].length; j++)
			{
				matrix_transpose[j][i] = matrix[i][j];
			}
		}
		
		return matrix_transpose;
		
	}
	
	public static double[][] transpose(double[] vector){
		double[][] vector_transpose = new double[1][vector.length];
		
		for(int i = 0; i < vector.length; i++)
		{
			
			vector_transpose[0][i] = vector[i];
			
		}
		
		return vector_transpose;
		
	}
	
	public static double[][] multiplyScalMatr(double scalar, double[][] m){
		double[][] result = new double[m.length][m[0].length];
		
		for(int i = 0; i < m.length; i++)
			for(int j = 0; j < m[0].length; j++)
				result[i][j] = scalar*m[i][j];

		return result;

	}

	public static double[] multiplyScalMatr(double scalar, double[] m){
		
		double[] result = new double[m.length];
		
		for(int i = 0; i < m.length; i++)
				result[i] = scalar*m[i];

		return result;

	}
	
	public static double[][] subtractMatrices(double[][] m1, double[][] m2){

		/*
		 * It returns m = m1 - m2.
		 * 
		 */

		double[][] m = new double[m1.length][m1[0].length];

		if(m2.length == m1.length && m2[0].length == m1[0].length){

			for(int i = 0; i < m1.length; i++)
				for(int j = 0; j < m1[0].length; j++)
					m[i][j] = m1[i][j] - m2[i][j];

		}else{
			System.out.println("Matrices dimension do not match.");
		}

		return m;

	}
	
	public static double[] subtractVectors(double[] m1, double[] m2){
		double[] result = new double[m1.length];
		
		for(int i = 0; i < m1.length; i++){
			result[i] = m1[i] - m2[i];
		}
		
		return result;
		
	}

	
	
	public static double[] crossProduct(double[] A, double[] B)
	{
		double[] result = {0,0,0};
		
		if(A.length != 3 || B.length != 3)
		{
			return result;
		}
		else
		{
			result[0] = A[1]*B[2] - A[2]*B[1];
			result[1] = A[2]*B[0] - A[0]*B[2];
			result[2] = A[0]*B[1] - A[1]*B[0];
			
		}
		return result;	
	}
	
	public static double[][] addMatrices(double[][] A, double[][] B)
	{
		int row = A.length;
		int column = A[0].length;
		double[][] C = new double[row][column];
		for(int i = 0; i < row; i++)
			for(int j = 0; j < column; j++)
				C[i][j] = A[i][j] + B[i][j];
		
		return C;
	}
	
	public static double[] addMatrices(double[] A, double[] B)
	{
		int row = A.length;
		double[] C = new double[row];
		for(int i = 0; i < row; i++)
			C[i] = A[i] + B[i];
		
		return C;
	}
	
	public static void displayMatrix(double[][] matrix){
		for(int i = 0; i < matrix.length; i++){			
			for(int j = 0; j < matrix[0].length; j++){
			
					System.out.print(matrix[i][j] + " ");
				}
			System.out.println(" ");
		}
	}
	
	public static void displayVector(double[] matrix){
//		System.out.println(" ");
		for(int i = 0; i < matrix.length; i++){
			System.out.print(matrix[i] + " ");
		}
		System.out.println(" ");
	}
	
	public static void fill(double[][] matrix, double[][] values){
		for(int i = 0; i < matrix.length; i++){			
			for(int j = 0; j < matrix[0].length; j++){
				matrix[i][j] = values[i][j];
				}
		}
	}
	
	public static void fill(double[][] matrix, double value){
		for(int i = 0; i < matrix.length; i++){			
			for(int j = 0; j < matrix[0].length; j++){
				matrix[i][j] = value;
				}
		}
	}
	
	public static int[] getIndMaxValue3D(double[][][] matrix){
		
		double max = 0;
		int[] index = new int[3];
		
		for(int i = 0; i < matrix.length; i++)
			for(int j = 0; j < matrix[0].length; j++)
				for(int k = 0; k < matrix[0][0].length; k++)
					
					if(matrix[i][j][k] > max){
						if(matrix[i][j][k] == 0)
							continue;
						else{
						max = matrix[i][j][k];
						index[0] = i;
						index[1] = j;
						index[2] = k;
						}
					}
		
		return index;
		
	}
	
	public static int[] getIndMinValue3D(double[][][] matrix){
		
		double min = Double.MAX_VALUE;
		int[] index = new int[3];
		
		for(int i = 0; i < matrix.length; i++)
			for(int j = 0; j < matrix[0].length; j++)
				for(int k = 0; k < matrix[0][0].length; k++)
					
					if(matrix[i][j][k] < min){
						min = matrix[i][j][k];
						index[0] = i;
						index[1] = j;
						index[2] = k;
					}
		
		return index;
		
	}
	
}
