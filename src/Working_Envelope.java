/*
 * It defines the working envelope.
 */


public class Working_Envelope {

	/*
	 * array with x, y, z minimum values [m]
	 */
	private double[] min_values;
	/*
	 * array with x, y, z maximum values [m]
	 */
	private double[] max_values;

	private double resolution;

	public static void main(String[] args){

	}

	/*
	 * Constructor considering only the max values of the volume. In this case you
	 * consider null min_values array and a default resolution = 0.05 
	 */
	
	public Working_Envelope(double[] max){

		min_values = new double[3];
		max_values = new double[3];

		for(int i = 0; i < min_values.length; i++){
				min_values[i] = 0;
				max_values[i] = max[i];
			}

		resolution = (double) 0.05;

	}
	
	/*
	 * Constructor considering also min_values.
	 */
	public Working_Envelope(double[] min, double[] max){

		min_values = new double[3];
		max_values = new double[3];

		for(int i = 0; i < min_values.length; i++){
				min_values[i] = min[i];
				max_values[i] = max[i];
			}

		resolution = (double) 0.05;

	}
	
	/*
	 * Constructor considering also resolution.
	 */
	public Working_Envelope(double[] min, double[] max, double res){

		min_values = new double[3];
		max_values = new double[3];

		for(int i = 0; i < min_values.length; i++){
				min_values[i] = min[i];
				max_values[i] = max[i];
			}

		resolution = res;

	}

	public double[] getMinValues(){

		return min_values;

	}

	public double[] getMaxValues(){

		return max_values;

	}

	public double getResolution(){

		return resolution;

	}

	public void setResolution(double res){

		resolution = res;

	}


}
