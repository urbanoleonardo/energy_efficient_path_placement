import java.util.ArrayList;
import java.util.List;

/*
 * TODO it will be necessary to address the problem that arises when the robot has more or less
 * than 6 degrees of freedom. The value that is hard coded here should be somehow picked from 
 * the class ROBOT that has DOF as attribute.
 * 
 */
public class Configuration {
	private List<Double> jointValues;

	public Configuration()
	{
		this.jointValues = new ArrayList<Double>();

	}
	public Configuration(float[] joints_values)
	{
		int length = joints_values.length;

		for(int i = 0; i < length; i++)
		{
			this.jointValues.add((double) joints_values[i]);
		}
	}

	public Configuration(double[] joints_values)
	{
		int length = joints_values.length;

		for(int i = 0; i < length; i++)
		{
			this.jointValues.add(joints_values[i]);
		}
	}

	public Configuration(List<Double> joints_values)
	{
		this.jointValues.addAll(joints_values);

	}

	public List<Double> getJointValues(){
		return this.jointValues;
	}
	
	public double getJointValue(int jointNumber){
		return this.jointValues.get(jointNumber-1);
	}
	
	public void setJointValue(int jointNumber, double jointValue){
		this.jointValues.set(jointNumber-1, jointValue);
	}
	
	private void Limitation() //Does it have to be in ROBOT class or not?
	{
		/*
		 * TODO
		 */
	}

}
