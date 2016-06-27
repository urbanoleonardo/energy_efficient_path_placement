import org.jgap.Chromosome;
import org.jgap.FitnessFunction;
import org.jgap.IChromosome;

public class MinimizingMakeChangeFitnessFunction extends FitnessFunction {

	private final Target target1;
	private final Target target2;
	private final Robot robot;

	public MinimizingMakeChangeFitnessFunction (Target initialTarget, Target finalTarget, Robot robot){

		target1 = new Target (initialTarget.getHomMatrix());
		target2 = new Target (finalTarget.getHomMatrix());
		this.robot = robot;
	}


	public double evaluate( IChromosome a_subject )
	{
		int delta = Integer.MAX_VALUE;
		
		double x = getPositionDouble(a_subject, 0);
		double y = getPositionDouble(a_subject, 1);
		double z = getPositionDouble(a_subject, 2);

		double[] currPos = {x,y,z};

		Target newStart = new Target(this.target1.getHomMatrix());
		Target newEnd = new Target(this.target2.getHomMatrix());

		newStart.translateTarget(currPos);
		newEnd.translateTarget(currPos);

//		OnlinePlannerCallable thread = new OnlinePlannerCallable(newStart, newEnd, robot, new Point3D(x,y,z));
		PlannerList thread = new PlannerList(newStart, newEnd, robot, new Point3D(x,y,z));

		EnergyPoint result = thread.call();
		
//		System.out.println("Energy of this gene : " + result.getEnergy());


		double fitness = result.getEnergy() == 0 ? delta : result.getEnergy();

		if(result.getEnergy() != 0){
			GATest.eCloud.add(result);
		}
		
		return (delta - fitness);
	}

	public static double getPositionDouble(IChromosome a_potentialSolution, int a_position){

		Integer spaceCoordinate = (Integer) a_potentialSolution.getGene(a_position).getAllele();

		double spaceCoordinateD = (double) spaceCoordinate;

		return spaceCoordinateD/100;
	}

	public static double getEnergyAtGene (IChromosome a_potentialSolution, Target target1, Target target2, Robot robot){

		double x = getPositionDouble(a_potentialSolution, 0);
		double y = getPositionDouble(a_potentialSolution, 1);
		double z = getPositionDouble(a_potentialSolution, 2);

		double[] currPos = {x,y,z};

		Target newStart = new Target(target1.getHomMatrix());
		Target newEnd = new Target(target2.getHomMatrix());

		newStart.translateTarget(currPos);
		newEnd.translateTarget(currPos);

//		OnlinePlannerCallable thread = new OnlinePlannerCallable(newStart, newEnd, robot, new Point3D(x,y,z));
		PlannerList thread = new PlannerList(newStart, newEnd, robot, new Point3D(x,y,z));
		
		EnergyPoint result = thread.call();
		
		
		return result.getEnergy();
	}
}
