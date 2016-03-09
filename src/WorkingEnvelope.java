import java.math.*;

public class WorkingEnvelope {

	/*
	 * It defines the working envelope, that is the volume in which you are
	 * going to shift the initial position of the trajectory
	 */

	private double[] minValues;

	private double[] maxValues;

	private double resolution;

	// Part of CONSTRUCTORS

	/*
	 * Constructor considering only the max values of the volume. In this case
	 * you consider null minValues array and a default resolution = 0.05 [m]
	 */

	public WorkingEnvelope (double[] max) {

		minValues = new double[3];

		maxValues = new double[3];

		for (int i = 0; i < minValues.length; i++) {

			minValues[i] = 0;

			maxValues[i] = max[i];

		}

		resolution = 5 / 100.0;

	}

	/*
	 * Constructor in which you define also minValues.
	 */

	public WorkingEnvelope (double[] min, double[] max) {

		minValues = new double[3];
		maxValues = new double[3];

		for (int i = 0; i < minValues.length; i++) {

			minValues[i] = min[i];
			maxValues[i] = max[i];

		}

		resolution = 5 / 100.0;

	}

	/*
	 * Constructor in which you define also the resolution.
	 */

	public WorkingEnvelope(double[] min, double[] max, double r) {

		minValues = new double[3];
		maxValues = new double[3];

		int order = 0;
		int n = 0;
		
		String s = null;

		for (int i = 0; i < minValues.length; i++) {

			minValues[i] = min[i];
			maxValues[i] = max[i];

		}

		/*
		 * Since we are using doubles, we have to be aware that numbers are not
		 * always what they are supposed to be EX: 0.05 could be 0.0499999999 In
		 * order to avoid problems with array's indexes when we are going to
		 * scan the working envelope etc., I have to ensure that 0.05 is
		 * actually 0.5! To do so, I make the following calculations
		 */

		if (r < 1.0) {

			s = Double.toString(r);
			
			/*
			 * I first store in "order" the number of decimals
			 */
			order = s.length() - 2;
			
			/*
			 * then I store in "n" the number that has to be scaled of
			 * 10^(-order)
			 */
			n = (int) (r * Math.pow(10, order));
			
			resolution = n / Math.pow(10, order);

		} else

			resolution = r;

	}

	// END of CONSTRUCTORS
	
	// MAIN method
	
	public static void main(String[] args) {
		
	}
	
	// END of MAIN method

	// Part of GET methods

	public double[] getMinValues() {

		return minValues;

	}

	public double[] getMaxValues() {

		return maxValues;

	}

	public double getResolution() {

		return resolution;

	}

	public int[] getSize() {

		int[] size = new int[3];

		for (int i = 0; i < 3; i++)

			size[i] = (int) (Math.ceil((this.maxValues[i] - this.minValues[i]) / this.resolution) + 1);

		return size;

	}

	public double[][][][] getWorkingEnvelopeMatrix() {

		/*
		 * It creates a 3D matrix of arrays in which there will be stored the
		 * values of x, y, z for each point of the working envelope
		 */

		int[] size = this.getSize();
		double[][][][] m = new double[size[0]][size[1]][size[2]][3];
		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		for (int x = 0; x < m.length; x++)
			for (int y = 0; y < m[0].length; y++)
				for (int z = 0; z < m[0][0].length; z++){

						/*
						 * In order to ensure a certain precision (2 decimals)
						 * for each double stored in the matrix, I am going to
						 * use a BigDecimal buffer. BigDecimals gives you the
						 * possibility to set the resolution of the input number
						 * as you prefer. I then get the doubleValue of the
						 * buffer and I store it in the matrix.
						 * 
						 */
						
						buffer[0] = new BigDecimal(minValues[0] + x * resolution, mt);
						buffer[1] = new BigDecimal(minValues[1] + y * resolution, mt);
						buffer[2] = new BigDecimal(minValues[2] + z * resolution, mt);
						
						m[x][y][z][0] = buffer[0].doubleValue();
						m[x][y][z][1] = buffer[1].doubleValue();
						m[x][y][z][2] = buffer[2].doubleValue();

					}

		return m;

	}

	// END of GET methods

	// SET methods

	public void setResolution(double res) {

		resolution = res;

	}

	// END of SET methods

	// OTHER PUBLIC methods

	// END of OTHER PUBLIC methods

}
