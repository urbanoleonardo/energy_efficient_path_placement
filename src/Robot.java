import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/*
 * TODO it defines all the parameters of the robot
 * 
 */

public class Robot {

	private String model;
	private int dof;
	private float[][] link_length;
	private float[][] link_masses;
	private float[][] inertia_matr;
	private float[][] inertia_tens;
	private float[][] joint_limits;
	Target location;
	
	public static void main(String[] args){
		
	}
	
	/* 
	 * you insert the model, the function looks for it in the cataloge.
	 * If it finds it, it initializes the robot with the data found
	 * in the cataloge, otherwise it gives you 2 choices: you can either
	 * get the data from a file or insert them manually. (?)
	 */
	public Robot(String model_path){
		
		String section_limiter = "**--";
		String comment_limiter = "%%";
		String s;
	    
		int section = 0;
		
		try{
		FileReader f;
		f = new FileReader(model_path);
		BufferedReader b;
		b = new BufferedReader(f);
		
		while(true){
			s = b.readLine();
			if(s == null)
			{
				break;
			}
			
			if(s.startsWith(comment_limiter) || s.isEmpty())
			{
				continue;
			}
			else
			{
				if(section == 0)
				{
					this.model = s;
				}
				if(s.startsWith(section_limiter))
				{
					
					NotifyActualSection(section);
					System.out.println(s);
					CopyRobotParameters(b, section);
					section++;
				}
				
			}
			
		}
		}
		catch (IOException e){
			System.out.println("Error in opening the file");
		}
		
	}
	
	
	private void NotifyActualSection(int section)
	{
		String sec;
		
		switch(section){
		case 0:
			sec = "Model of the Robot";
			break;
		case 1: 
			sec = "Degrees of Freedom";
			break;
		case 2: 
			sec = "Move to the wirst";
			break;
		case 3: 
			sec = "Joints working range";
			break;
		default:
			sec = " ";
			break;
		}
		System.out.println("The current section being compiled is: " + sec);
	}
	private void CopyRobotParameters (BufferedReader b, int section) throws IOException
	{
		String s;
		String comment_limiter = "%%";
		String section_limiter = "**--";
		
		int row_number = 0; //it's used to remember which row of the matrix we are copying
		//Section where I still switch but I instantiate the attributes of the robot
		//to be filled with data from the file
		
		switch(section){
		case 2: 
			this.link_length = new float[this.dof][2];
			 
			break;
		case 3:
			this.link_masses = new float[this.dof][1];
			break;
		case 4:
			this.joint_limits = new float[this.dof][2];
			break;
		
		default:
			break;
		}
		
		
		
		while(true){
		s = b.readLine();
		if(s == null)
		{
			break;
		}
		if(s.startsWith(comment_limiter) || s.isEmpty())
		{
			continue;
		}
		if(s.startsWith(section_limiter))
		{
			break;
		}
		if(section == 0)
		{
			this.model = s;
			System.out.println(this.model);
			continue;
		}
		String[] tokens = s.split(" ");
		//System.out.println("the number of tokens in this string is " + tokens.length);
		for(int i = 0; i<tokens.length; i++ )
		{
			System.out.println(tokens[i]);
		}
		switch(section){
		case 1: 
			this.dof = Integer.parseInt(tokens[0]);
			break;
		case 2:
			this.link_length[row_number][0] = Float.parseFloat(tokens[0]);
			this.link_length[row_number][1] = Float.parseFloat(tokens[1]);
			row_number++;
			break;
		
		default:
			break;
		}
		}
	}
	public int getDOF(){
		return this.dof;
	}
	
	public float[][] inertiaTens(){
		float[][] ajeje = new float[4][4];
		return ajeje;
	}
	
}
