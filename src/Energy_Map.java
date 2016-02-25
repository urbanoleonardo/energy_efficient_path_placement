import java.util.Scanner;

public class Energy_Map {

	public static void main(String[] args){
		
		Scanner input = new Scanner(System.in);
		
		Working_Envelope we;
		Path p;
		Robot r;
		
		double[] min_values = new double[3];
		double[] max_values = new double[3];
		double resolution;
		
		double[] in_pos = new double[3];
		double[][] in_rot = new double[3][3];
		Target initial_position;
		double[] fin_pos = new double[3];
		double[][] fin_rot = new double[3][3];
		Target final_position;
		double t_sample;
		double x_sample;
		
		String model;
		
		double[] current = new double[3];
		
		/*
		 * Getting working envelope data from user
		 */
		
		System.out.println("WORKING ENVELOPE");
		
		for(int i = 0; i < min_values.length; i++){
		System.out.println("Minimum values[" + i + "]: ");
		min_values[i] = input.nextDouble();
		}
		
		for(int i = 0; i < max_values.length; i++){
		System.out.println("Maximum values[" + i + "]: ");
		max_values[i] = input.nextDouble();
		}
		
		System.out.println("Resolution:");
		resolution = input.nextDouble();
		
		we = new Working_Envelope(min_values, max_values, resolution);
		
		/*
		 * Getting path data from user and path interpolation
		 */
		
		System.out.println("PATH");
		
		System.out.println("Initial position");
		System.out.println("position: ");
		for(int i = 0; i < in_pos.length; i++){
			System.out.println("position[" + i + "] = ");
			in_pos[i] = input.nextDouble();
		}
		System.out.println("rotation:");
		for(int i = 0; i < in_rot.length; i++)
			for(int j = 0; j < in_rot[0].length; j++){			
				System.out.println("rotation[" + i + "]" + "[" + j + "] = ");
				in_rot[i][j] = input.nextDouble();			
		}
		
		initial_position = new Target(in_pos, in_rot);
		
		System.out.println("Final position");
		System.out.println("position:");
		for(int i = 0; i < fin_pos.length; i++){
			System.out.println("position[" + i + "] = ");
			fin_pos[i] = input.nextDouble();
		}
		System.out.println("rotation:");
		for(int i = 0; i < fin_rot.length; i++)
			for(int j = 0; j < fin_rot[0].length; j++){			
				System.out.println("rotation[" + i + "]" + "[" + j + "] = ");
				fin_rot[i][j] = input.nextDouble();			
		}
		
		final_position = new Target(fin_pos, fin_rot);
		
		System.out.println("Time resolution: ");
		t_sample = input.nextDouble();
		
		System.out.println("Space resolution: ");
		x_sample = input.nextDouble();
		
		p = new Path(initial_position, final_position, t_sample, x_sample);
		
		p.interpolate();
		
		/*
		 * Getting robot data from user
		 */
		
		System.out.println("ROBOT");
		
		System.out.println("Model: ");
		model = input.nextLine();
		
		r = new Robot(model);
		
		/*
		 * 3 loops in order to scan every single point of the Working Envelope
		 */
		
		for(double x = we.getMinValues()[1]; x <= we.getMaxValues()[1]; x += we.getResolution())
			for(double y = we.getMinValues()[2]; y <= we.getMaxValues()[2]; y += we.getResolution())
				for(double z = we.getMinValues()[3]; z <= we.getMaxValues()[3]; z += we.getResolution()){
					
					current[1] = x;
					current[2] = y;
					current[3] = z;
					
					p.translateTrajectory(current);
					
				}
	
	
	}
	
	//Part of CONSTRUCTORS
	//END part of CONSTRUCTORS

	//Part of GET methods
	//END of GET methods

	//Part of SET methods
	//END of SET methods

	//OTHER PUBLIC methods
	//END of PUBLIC methods
	
	//Part of PRIVATE methods
	//END of PRIVATE methods
	
}
