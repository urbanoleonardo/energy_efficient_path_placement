/*
 * It defines the working envelope.
 */


public class Working_Envelope {

	/*
	 * 3x1 array with x, y, z minimum values [m]
	 */
	private float[][] min_values;
	/*
	 * 3x1 array with x, y, z maximum values [m]
	 */
	private float[][] max_values;

	private float resolution;

	public static void main(String[] args){

	}

	/*
	 * Constructor considering only the max values of the volume. In this case you
	 * consider null min_values array and a defaul resolution = 0.05 */
	public void Working_Envelope(float[][] max){

		min_values = new float[3][1];
		max_values = new float[3][1];

		for(int i = 0; i < min_values.length; i++)
			for(int j = 0; j < min_values[0].length; j++){
				min_values[i][j] = 0;
				max_values[i][j] = max[i][j];
			}

		resolution = (float) 0.05;

	}

	/*
	 * Constructor considering also min_values.
	 */
	public void Working_Envelope(float[][] min, float[][] max, float res){

		min_values = new float[3][1];
		max_values = new float[3][1];

		for(int i = 0; i < min_values.length; i++)
			for(int j = 0; j < min_values[0].length; j++){
				min_values[i][j] = min[i][j];
				max_values[i][j] = max[i][j];
			}

		resolution = res;

	}

	public float[][] getMinValues(){

		return min_values;

	}

	public float[][] getMaxValues(){

		return max_values;

	}

	public float getResolution(){

		return resolution;

	}

	public void setResolution(float res){

		resolution = res;

	}


}
