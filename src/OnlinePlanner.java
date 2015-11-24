/*
TO DO 
decide how to get the current position of the robot, so the starting point that is
needed for the inverse kinematics. 
*/
public class OnlinePlanner implements Runnable{

	private Target currPosition;
	private int targetsLength;
	private Target[] targets;
	private Robot robot;
	public OnlinePlanner(Target[] targets, Robot robot)
	{
		
		targetsLength = copytargets(targets);
		this.robot = robot;
	}
	
	public void run(){
		if (targetsLength == 1)
	{
	
	}
	else
	{
		
	}
		
		
	}
	private int copytargets(Target[] targetsTOstore)
	{
		int length = targetsTOstore.length;
		this.targets = new Target[length];
		for(int i = 0; i<length; i++)
		{
			this.targets[i] = targetsTOstore[i];
		}
		return length;
	}
	private void Online_Kinematics(Target Current_position, Target Next_position)
	{
		
	}
	private void Online_Dynamics()
	{
		
	}
}
