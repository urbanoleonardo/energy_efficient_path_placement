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
		this.jointValues = new ArrayList<Double>();
		int length = joints_values.length;

		for(int i = 0; i < length; i++)
		{
			this.jointValues.add((double) joints_values[i]);
		}
	}

	public Configuration(double[] joints_values)
	{
		this.jointValues = new ArrayList<Double>();
		int length = joints_values.length;

		for(int i = 0; i < length; i++)
		{
			this.jointValues.add(joints_values[i]);
		}
	}

	public Configuration(List<Double> jointValues)
	{
		this.jointValues = new ArrayList<Double>();
		this.jointValues.addAll(jointValues);

	}

	public List<Double> getJointValues(){
		return this.jointValues;
	}

	public int size(){
		return this.jointValues.size();
	}

	public double getJointValue(int jointNumber){
		return this.jointValues.get(jointNumber-1);
	}

	public void setJointValue(int jointNumber, double jointValue){

		if(jointNumber > this.jointValues.size()){
			this.jointValues.add(jointValue);
		} else{
			this.jointValues.set(jointNumber-1, jointValue);
		}

	}

	public void setJointValue(List<Double> jointValues){

		this.jointValues.addAll(jointValues);

	}

	public void setJointValue(double[] jointValues){
		this.jointValues = new ArrayList<Double>();
		
		for(int i = 0; i < jointValues.length; i++){
			this.jointValues.add(jointValues[i]);
		}

	}
}
