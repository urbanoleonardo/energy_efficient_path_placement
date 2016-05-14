
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.LinkedList;
import java.util.List;

public class WorkingEnvelope {

	/*
	 * It defines the working envelope, that is the volume in which you are
	 * going to shift the initial position of the trajectory
	 */

	private final double xMin;
	private final double yMin;
	private final double zMin;
	
	private final double xMax;
	private final double yMax;
	private final double zMax;

	private double resolution;
	
	private List<Point3D> weList;

	// Part of CONSTRUCTORS

	/*
	 * Constructor considering only the max values of the volume. In this case
	 * you consider null minValues array and a default resolution = 0.05 [m]
	 */

	public WorkingEnvelope (double[] max) {

		if(max == null || max.length != 3){
			
			throw new IllegalArgumentException();
			
		}
		
		xMin = 0;
		yMin = 0;
		zMin = 0;
		
		xMax = max[0];
		yMax = max[1];
		zMax = max[2];

		resolution = 5 / 100.0;
		
//		weList = buildWeList();

	}

	/*
	 * Constructor in which you define also minValues.
	 */

	public WorkingEnvelope (double[] min, double[] max) {

		xMin = min[0];
		yMin = min[1];
		zMin = min[2];
		
		xMax = max[0];
		yMax = max[1];
		zMax = max[2];

		resolution = 5 / 100.0;
		
//		weList = buildWeList();

	}

	/*
	 * Constructor in which you define also the resolution.
	 */

	public WorkingEnvelope(double[] min, double[] max, double r) {

		xMin = min[0];
		yMin = min[1];
		zMin = min[2];
		
		xMax = max[0];
		yMax = max[1];
		zMax = max[2];

		int order = 0;
		int n = 0;
		
		String s = null;

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
//			order = s.length() - 2;
			order = s.length() - s.indexOf('.') - 1;
			
			/*
			 * then I store in "n" the number that has to be scaled of
			 * 10^(-order)
			 */
			n = (int) (r * Math.pow(10, order));
			
			resolution = n / Math.pow(10, order);

		} else

			resolution = r;
		
//		weList = buildWeList();

	}

	// END of CONSTRUCTORS
	
	// MAIN method
	
	public static void main(String[] args) {
		
	}
	
	// END of MAIN method

	// Part of GET methods

	public double[] getMinValues() {

		return new double[] {xMin, yMin, zMin};

	}

	public double[] getMaxValues() {

		return new double[] {xMax, yMax, zMax};

	}

	public double getResolution() {

		return resolution;

	}
	
	public List<Point3D> getWeList(){
		
		return weList;
		
	}

	public int[] getSize() {

		int[] size = new int[3];
		
		size[0] = (int) (Math.ceil((this.xMax - this.xMin) / this.resolution) + 1);
		size[1] = (int) (Math.ceil((this.yMax - this.yMin) / this.resolution) + 1);
		size[2] = (int) (Math.ceil((this.zMax - this.zMin) / this.resolution) + 1);

		return size;

	}

	private List<Point3D> buildWeList() {

		/*
		 * It creates a 3D matrix of arrays in which there will be stored the
		 * values of x, y, z for each point of the working envelope
		 */

		int[] size = this.getSize();
		
		List<Point3D> weList = new LinkedList<Point3D>();
		
		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		for (int x = 0; x < size[0]; x++)
			for (int y = 0; y < size[1]; y++)
				for (int z = 0; z < size[2]; z++){

						/*
						 * In order to ensure a certain precision (2 decimals)
						 * for each double stored in the matrix, I am going to
						 * use a BigDecimal buffer. BigDecimals gives you the
						 * possibility to set the resolution of the input number
						 * as you prefer. I then get the doubleValue of the
						 * buffer and I store it in the matrix.
						 * 
						 */
						
						buffer[0] = new BigDecimal(xMin + x * resolution, mt);
						buffer[1] = new BigDecimal(yMin + y * resolution, mt);
						buffer[2] = new BigDecimal(zMin + z * resolution, mt);
						
						Point3D p = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());
						weList.add(p);
						
					}
		
		return weList;

	}

	// END of GET methods

	// SET methods

	public void setResolution(double res) {

		resolution = res;
		
		weList = buildWeList();

	}

	// END of SET methods

	// OTHER PUBLIC methods

	// END of OTHER PUBLIC methods

}

