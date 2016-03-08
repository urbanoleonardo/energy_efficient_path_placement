/*
 * It defines the working envelope.
 */
import java.math.*;

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

		resolution = 5/100.0;

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

		resolution = 5/100.0;

	}

	/*
	 * Constructor considering also resolution.
	 */
	public Working_Envelope(double[] min, double[] max, double r){

		min_values = new double[3];
		max_values = new double[3];

		String s = null;
		int order = 0;
		int n = 0;

		for(int i = 0; i < min_values.length; i++){
			min_values[i] = min[i];
			max_values[i] = max[i];
		}

		if(r < 1.0){
			s = Double.toString(r);
//			System.out.println("String s = ");
//			System.out.println(s);
			order = s.length() - 2;
//			System.out.println("order = ");
//			System.out.println(order);
			n = (int) (r*Math.pow(10, order));
//			System.out.println("n = ");
//			System.out.println(n);
			resolution = n/Math.pow(10, order);
//			System.out.println("Resolution = ");
//			System.out.println(resolution);
		}else
			resolution = r;

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

	public int[] getSize(){

		/*
		 * To round up to the next decimal I have used the formula a/b = (a + b - 1)/b (more efficient..)
		 */

		int[] size = new int[3];
		double[] s = new double[3];

		for(int i = 0; i < 3; i++)
			//			size[i] =  (int) ((this.max_values[i] - this.min_values[i] + this.resolution - 1)/this.resolution);
			size[i] =  (int) (Math.ceil((this.max_values[i] - this.min_values[i])/this.resolution) + 1);

		return size;

	}
	
//	public double getIncrement(int i){
//		
////		String s = null;
////		int order = 0;
////		int n = 0;
//		double increment;
////		
////		s = Double.toString(resolution);
////		order = s.length() - 2;
////		n = (int) (resolution*Math.pow(10, order));
//		
//		increment = i*resolution;
//		
//		return increment;
//		
//	}
	
	public double[][][][] getWorkingEnvelopeMatrix(){
		
		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);
		int[] size = this.getSize();
		double[][][][] m = new double[size[0]][size[1]][size[2]][3];
		
		for(int x = 0; x < m.length; x++)
			for(int y = 0; y < m[0].length; y++)
				for(int z = 0; z < m[0][0].length; z++){
						
					buffer[0] = new BigDecimal(min_values[0] + x*resolution, mt);
					buffer[1] = new BigDecimal(min_values[1] + y*resolution, mt);
					buffer[2] = new BigDecimal(min_values[2] + z*resolution, mt);
						
						m[x][y][z][0] = buffer[0].doubleValue();
						m[x][y][z][1] = buffer[1].doubleValue();
						m[x][y][z][2] = buffer[2].doubleValue();
						
				}
		
		return m;
		
	}

	public void setResolution(double res){

		resolution = res;

	}


}
