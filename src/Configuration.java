
/*
 * TODO it will be necessary to address the problem that arises when the robot has more or less
 * than 6 degrees of freedom. The value that is hard coded here should be somehow picked from 
 * the class ROBOT that has DOF as attribute.
 * 
 */
public class Configuration {
private float[] joints_values;

public Configuration(int DOF)
{
	this.joints_values = new float[DOF];
	
}
public Configuration(float[] joints_values, int DOF)
{
	int length = joints_values.length;
	if(length > DOF)
	{
		System.out.println("More joints values than the actual degrees of freedom specified have been given.");
		
	}
	for(int i = 0; i < length; i++)
	{
		this.joints_values[i] = joints_values[i];
	}
}

private void Limitation() //Does it have to be in ROBOT class or not?
{
	/*
	 * TODO
	 */
}

}
