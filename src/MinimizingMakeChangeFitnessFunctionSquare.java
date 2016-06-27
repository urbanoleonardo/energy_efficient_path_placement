import org.jgap.Chromosome;
import org.jgap.FitnessFunction;
import org.jgap.IChromosome;

public class MinimizingMakeChangeFitnessFunctionSquare extends FitnessFunction{

	private final Target point1;
	private final Target point2;
	private final Target point3;
	private final Target point4;
	private final Robot robot;

	public MinimizingMakeChangeFitnessFunctionSquare(SquarePath path, Robot robot){

		point1 = path.getPathPositions()[0];
		point2 = path.getPathPositions()[1];
		point3 = path.getPathPositions()[2];
		point4 = path.getPathPositions()[3];
		
		this.robot = robot;
	}


	public double evaluate( IChromosome a_subject )
	{
		int delta = Integer.MAX_VALUE;
		
		double x = getPositionDouble(a_subject, 0);
		double y = getPositionDouble(a_subject, 1);
		double z = getPositionDouble(a_subject, 2);

		double[] currPos = {x,y,z};

		Target point1 = new Target(this.point1.getHomMatrix());
		Target point2 = new Target(this.point2.getHomMatrix());
		Target point3 = new Target(this.point3.getHomMatrix());
		Target point4 = new Target(this.point4.getHomMatrix());
		
		point1.translateTarget(currPos);
		point2.translateTarget(currPos);
		point3.translateTarget(currPos);
		point4.translateTarget(currPos);
		
		Path newSquarePath = new SquarePath(point1,point2,point3,point4);
		newSquarePath.setMaxAcc(1.0);
		newSquarePath.setMaxVel(0.2);


		
		
		
//		OnlinePlannerCallable thread = new OnlinePlannerCallable(newStart, newEnd, robot, new Point3D(x,y,z));
		PlannerList thread = new PlannerList(newSquarePath, robot, new Point3D(x,y,z));

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

	public static double getEnergyAtGene (IChromosome a_potentialSolution, SquarePath path , Robot robot){

		double x = getPositionDouble(a_potentialSolution, 0);
		double y = getPositionDouble(a_potentialSolution, 1);
		double z = getPositionDouble(a_potentialSolution, 2);

		double[] currPos = {x,y,z};

		Target point1 = new Target (path.getPathPositions()[0].getHomMatrix());
		Target point2 = new Target (path.getPathPositions()[1].getHomMatrix());
		Target point3 = new Target (path.getPathPositions()[2].getHomMatrix());
		Target point4 = new Target (path.getPathPositions()[3].getHomMatrix());
		
		point1.translateTarget(currPos);
		point2.translateTarget(currPos);
		point3.translateTarget(currPos);
		point4.translateTarget(currPos);
		
		Path newSquarePath = new SquarePath(point1,point2,point3,point4); 
		newSquarePath.setMaxAcc(1.0);
		newSquarePath.setMaxVel(0.2);

//		OnlinePlannerCallable thread = new OnlinePlannerCallable(newStart, newEnd, robot, new Point3D(x,y,z));
		PlannerList thread = new PlannerList(newSquarePath, robot, new Point3D(x,y,z));
		
		EnergyPoint result = thread.call();
		
		
		return result.getEnergy();
	}
}
