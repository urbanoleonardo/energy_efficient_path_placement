import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jgap.Chromosome;
import org.jgap.Configuration;
import org.jgap.FitnessFunction;
import org.jgap.Gene;
import org.jgap.GeneticOperator;
import org.jgap.Genotype;
import org.jgap.IChromosome;
import org.jgap.IChromosomePool;
import org.jgap.InvalidConfigurationException;
import org.jgap.Population;
import org.jgap.audit.*;
import org.jgap.impl.ChromosomePool;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.GreedyCrossover;
import org.jgap.impl.IntegerGene;

import com.sun.j3d.utils.applet.MainFrame;

public class GATestOLD {

	final static int MAX_ALLOWED_EVOLUTIONS = 500;

	public static List<EnergyPoint> eCloud = new LinkedList<EnergyPoint>();
	
	public static void main(String[] args) throws Exception {
		/*
		 * This method is running a monitor to early stop the evolution.
		 * Moreover is using chromosomePool filled with the results from the grid, 
		 * that are also used to form the initial population.
		 */
//		GATest2();
//		GATestSquare();
//		hybridGATest1();
		hybridGATest2();
//		hybridGATestSquare();
//		discreteTestSquare();
		
	}
	

	public static void hybridGATest1() throws Exception {
		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";

		double[] min = { -1.0, -1.0, 0.0 };
		double[] max = { 1.0, 1.0, 1.0 };

		int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
		int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
		
		WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);


//		double[][] h1 = { { 0.0, 0.0, 1.0, 0.59627 }, 
//				{ 1.0, 0.0, 0.0, 0.0 }, 
//				{ 0.0, 1.0, 0.0, 0.66282 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
//		double[][] h2 = { { 0.0, 0.0, 1.0, 0.59627 },
//				{ 1.0, 0.0, 0.0, -0.21132 }, 
//				{ 0.0, 1.0, 0.0, 0.58041 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
		
		double[][] h1 = { 	
				{ 1.0, 0.0, 0.0, 0.0 }, 
				{ 0.0, 1.0, 0.0, 0.2 }, 
				{ 0.0, 0.0, -1.0, 0.0 },
				{ 0.0, 0.0, 0.0, 1.0 } };
		double[][] h2 = 	{ 	
				{ 1.0, 0.0, 0.0, 0.0 },
				{ 0.0, 1.0, 0.0, -0.2 }, 
				{ 0.0, 0.0, -1.0, 0.0 },
				{ 0.0, 0.0, 0.0, 1.0 } };


		Target testTarget1 = new Target(h1);
		Target testTarget2 = new Target(h2);
		Target[] targets = new Target[2];
		double[] backToCow = { - testTarget1.getPosition()[0],  - testTarget1.getPosition()[1],  - testTarget1.getPosition()[2]};
		testTarget1.translateTarget(backToCow);
		testTarget2.translateTarget(backToCow);
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;
		Path p = new Path(testTarget1, testTarget2, t_sample, x_sample);
		p.setMaxAcc(acceleration);
		p.setMaxVel(velocity);

		Robot robot = new Robot(xmlFilePath, true);

		List<EnergyPoint> initialEnergyGrid = createEnergyCloud(robot, p , workingEnvelope);
		eCloud.addAll(initialEnergyGrid);
		
		Configuration conf = new DefaultConfiguration();
		conf.setPreservFittestIndividual(true);
		conf.setKeepPopulationSizeConstant(false);
		
		
		FitnessFunction myFunc = new MinimizingMakeChangeFitnessFunction(testTarget1, testTarget2, robot);

		conf.setFitnessFunction( myFunc );
		
		Gene[] sampleGenes = new Gene[ 3 ];

		sampleGenes[0] = new IntegerGene(conf, minInt[0], maxInt[0] );  
		sampleGenes[1] = new IntegerGene(conf, minInt[1], maxInt[1] );  
		sampleGenes[2] = new IntegerGene(conf, minInt[2], maxInt[2] );  

		Chromosome sampleChromosome = new Chromosome(conf, sampleGenes );
		conf.setSampleChromosome( sampleChromosome );

		int populationSize = 5000;
//		int populationSize = 10000;

		conf.setPopulationSize( populationSize );

		
		ChromosomePool pool = new ChromosomePool();
		conf.setChromosomePool(pool);
		
		initializePool(pool, conf, sampleChromosome, initialEnergyGrid);
		Population pop = initializePopulation(conf, sampleChromosome, initialEnergyGrid);
		
		


		System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
		System.out.println("The population size is : " + pop.size());
		

//		Genotype population = Genotype.randomInitialGenotype( conf );
		Genotype population = new Genotype(conf, pop);
		population.fillPopulation(populationSize - pop.size());
		
		long startTime = System.currentTimeMillis();
		
		List monitors = new Vector();
		monitors.add(new TimedMonitor(30));
		monitors.add(new FitnessImprovementMonitor(3, 3, 0.1d));
		IEvolutionMonitor monitor = new ChainedMonitors(monitors, 2);

		//		for( int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++ )
		//		  {
		//		      population.evolve();
		//		  }
		//		
		//		IChromosome bestSolutionSoFar = population.getFittestChromosome();

		List<String> messages = population.evolve(monitor);  
		if (messages.size() > 0) {
			for (String msg : messages) {
				System.out.println("Message from monitor: " + msg+"\n");
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
	    System.out.println("Total evolution time: " + (endTime - startTime) + " ms");
		
		IChromosome bestSolutionSoFar = population.getFittestChromosome();
		System.out.println( "The best solution contained the following: " );

		System.out.println( MinimizingMakeChangeFitnessFunction.getEnergyAtGene(bestSolutionSoFar, testTarget1, testTarget2, robot) + " J totally required." );

		System.out.println("x = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 0 ));
		System.out.println("y = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 1 ));
		System.out.println("z = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 2 ));

		Map3D demo = new Map3D(xmlFilePath, eCloud);
		

	}
	public static void discreteTestSquare() throws Exception {
		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST SW.xml";

		double[] min = { -2.0, -2.0, 0.0 };
		double[] max = { 2.0, 2.0, 2.0 };

		int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
		int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
		
		WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);
		workingEnvelope.setResolution(0.1);
		
		

		double[][] vertex1 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{-1.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex2 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{-1.0, 0.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex3 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, -0.4},
				{-1.0, 0.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex4 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, -0.4},
				{-1.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		
		Target point1 = new Target(vertex1);
		Target point2 = new Target(vertex2);
		Target point3 = new Target(vertex3);
		Target point4 = new Target(vertex4);
		
		
	
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;

		Path squarePath = new SquarePath(point1,point2,point3,point4);
		squarePath.setMaxAcc(acceleration);
		squarePath.setMaxVel(velocity);
		squarePath.interpolate();
		
		Robot robot = new Robot(xmlFilePath, true);
		
		
		List<EnergyPoint> initialEnergyGrid = createEnergyCloudGeneral(robot, squarePath , workingEnvelope);
		eCloud.addAll(initialEnergyGrid);
		
		Map3D demo = new Map3D(xmlFilePath, eCloud);
		EnergyPoint minimum = findMinimum(initialEnergyGrid);
		System.out.println("The minimum consumption is : " + minimum.getEnergy());
		System.out.println("Located in ->  x = " + minimum.getX() + "  y = " + minimum.getY() + "  z = "+ minimum.getZ());
	}
	
	public static void GATest2() throws Exception {
//		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";
//		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST SW.xml";
		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST.xml";

		double[] min = { -2.0, -2.0, 0.0 };
		double[] max = { 2.0, 2.0, 2.0 };

		int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
		int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
		
		WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);
		workingEnvelope.setResolution(0.1);
		
//		double[][] h1 = { 	{ 1.0, 0.0, 0.0, 0.0 }, 
//				{ 0.0, 1.0, 0.0, 0.2 }, 
//				{ 0.0, 0.0, -1.0, 0.0 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
//		double[][] h2 = 	{ 	{ 1.0, 0.0, 0.0, 0.0 },
//				{ 0.0, 1.0, 0.0, -0.2 }, 
//				{ 0.0, 0.0, -1.0, 0.0 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
		double[][] h1 = { 
				{ 0.0, 0.0, 1.0, 0.0 }, 
				{ 0.0, 1.0, 0.0, 0.0 }, 
				{ -1.0, 0.0, 0.0, 0.0 },
				{ 0.0, 0.0, 0.0, 1.0 } };
		
		double[][] h2 = { 
				{ 0.0, 0.0, 1.0, 0.0 }, 
				{ 0.0, 1.0, 0.0, -0.4 }, 
				{ -1.0, 0.0, 0.0, 0.0 },
				{ 0.0, 0.0, 0.0, 1.0 } };
		
//		double[][] h2 = { 	
//				{ 1.0, 0.0, 0.0, -0.2 },
//				{ 0.0, 1.0, 0.0, -0.4 }, 
//				{ 0.0, 0.0, 1.0, 0.2 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
		
		

		Target testTarget1 = new Target(h1);double[][] vertex1 = {
				{0.0, 0.0, 1.0, 0.0},
				{1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex2 = {
				{0.0, 0.0, 1.0, 0.0},
				{1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex3 = {
				{0.0, 0.0, 1.0, 0.0},
				{1.0, 0.0, 0.0, -0.4},
				{0.0, 1.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex4 = {
				{0.0, 0.0, 1.0, 0.0},
				{1.0, 0.0, 0.0, -0.4},
				{0.0, 1.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		
		Target point1 = new Target(vertex1);
		Target point2 = new Target(vertex2);
		Target point3 = new Target(vertex3);
		Target point4 = new Target(vertex4);
		
		
		
		Target testTarget2 = new Target(h2);
		Target[] targets = new Target[2];
		double[] backToCow = { - testTarget1.getPosition()[0],  - testTarget1.getPosition()[1],  - testTarget1.getPosition()[2]};
		testTarget1.translateTarget(backToCow);
		testTarget2.translateTarget(backToCow);
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;
		Path p = new Path(testTarget1, testTarget2, t_sample, x_sample);
		p.setMaxAcc(acceleration);
		p.setMaxVel(velocity);

		Path squarePath = new SquarePath(point1,point2,point3,point4);
		squarePath.setMaxAcc(acceleration);
		squarePath.setMaxVel(velocity);
		squarePath.interpolate();
		
		Robot robot = new Robot(xmlFilePath, true);

//		List<EnergyPoint> initialEnergyGrid = createEnergyCloud(robot, p , workingEnvelope);
		
//		EnergyPoint minimum = findMinimum(initialEnergyGrid);
//		System.out.println("minimum with discrete -> x: " + minimum.getX() + " y: " + minimum.getY() + " z: " + minimum.getZ());
		//eCloud.addAll(initialEnergyGrid);
		Configuration conf = new DefaultConfiguration();
		
		conf.setPreservFittestIndividual(true);
		conf.setKeepPopulationSizeConstant(false);
		
		
		FitnessFunction myFunc = new MinimizingMakeChangeFitnessFunction(testTarget1, testTarget2, robot);

		conf.setFitnessFunction( myFunc );
		
		Gene[] sampleGenes = new Gene[ 3 ];
		
		sampleGenes[0] = new IntegerGene(conf, minInt[0], maxInt[0] );  
		sampleGenes[1] = new IntegerGene(conf, minInt[1], maxInt[1] );  
		sampleGenes[2] = new IntegerGene(conf, minInt[2], maxInt[2] );  
		
		Chromosome sampleChromosome = new Chromosome(conf, sampleGenes );
		conf.setSampleChromosome( sampleChromosome );

		int populationSize = 500;
//		int populationSize = 10000;

		conf.setPopulationSize( populationSize );
		//Population pop = initializePopulation(conf, sampleChromosome, initialEnergyGrid);
		
		Genotype population = Genotype.randomInitialGenotype(conf);
		//Genotype population = new Genotype(conf, pop);
		//population.fillPopulation(populationSize - pop.size());

		
		
		long startTime = System.currentTimeMillis();
		
		List monitors = new Vector();
		monitors.add(new TimedMonitor(120));
		monitors.add(new FitnessImprovementMonitor(1, 3, 0.01d));
		IEvolutionMonitor monitor = new ChainedMonitors(monitors, 2);

//				for( int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++ )
//				  {
//				      population.evolve();
//				  }
				
//				IChromosome bestSolutionSoFar = population.getFittestChromosome();

		List<String> messages = population.evolve(monitor);  
		if (messages.size() > 0) {
			for (String msg : messages) {
				System.out.println("Message from monitor: " + msg+"\n");
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
	    System.out.println("Total evolution time: " + (endTime - startTime) + " ms");
		
		IChromosome bestSolutionSoFar = population.getFittestChromosome();
		System.out.println( "The best solution contained the following: " );

		System.out.println( MinimizingMakeChangeFitnessFunction.getEnergyAtGene(bestSolutionSoFar, testTarget1, testTarget2, robot) + " J totally required." );

		System.out.println("x = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 0 ));
		System.out.println("y = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 1 ));
		System.out.println("z = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 2 ));

		Map3D demo = new Map3D(xmlFilePath, eCloud);
	}
	
	public static void GATestSquare() throws Exception {
		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST SW.xml";

		double[] min = { -2.0, -2.0, 0.0 };
		double[] max = { 2.0, 2.0, 2.0 };

		int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
		int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
		
		WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);
		
		

		double[][] vertex1 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{-1.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex2 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{-1.0, 0.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex3 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, -0.4},
				{-1.0, 0.0, 0.0, 0.4},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		double[][] vertex4 = {
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 1.0, 0.0, -0.4},
				{-1.0, 0.0, 0.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
		
		
//		double[][] vertex1 = {
//				{0.0, 0.0, 1.0, 0.59627},
//				{1.0, 0.0, 0.0, -0.21132},
//				{0.0, 1.0, 0.0, 0.38041},
//				{0.0, 0.0, 0.0, 1.0}
//		};
//		
//		double[][] vertex2 = {
//				{0.0, 0.0, 1.0, 0.49627},
//				{1.0, 0.0, 0.0, -0.21132},
//				{0.0, 1.0, 0.0, 0.38041},
//				{0.0, 0.0, 0.0, 1.0}
//		};
//		
//		double[][] vertex3 = {
//				{0.0, 0.0, 1.0, 0.49627},
//				{1.0, 0.0, 0.0, -0.31132},
//				{0.0, 1.0, 0.0, 0.38041},
//				{0.0, 0.0, 0.0, 1.0}
//		};
//		
//		double[][] vertex4 = {
//				{0.0, 0.0, 1.0, 0.59627},
//				{1.0, 0.0, 0.0, -0.31132},
//				{0.0, 1.0, 0.0, 0.38041},
//				{0.0, 0.0, 0.0, 1.0}
//		};
		
		Target point1 = new Target(vertex1);
		Target point2 = new Target(vertex2);
		Target point3 = new Target(vertex3);
		Target point4 = new Target(vertex4);
		
		
	
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;

		Path squarePath = new SquarePath(point1,point2,point3,point4);
		squarePath.setMaxAcc(acceleration);
		squarePath.setMaxVel(velocity);
		squarePath.interpolate();
		
		Robot robot = new Robot(xmlFilePath, true);
		
		Configuration conf = new DefaultConfiguration();
		
		conf.setPreservFittestIndividual(true);
		conf.setKeepPopulationSizeConstant(false);
		
		
		FitnessFunction myFunc = new MinimizingMakeChangeFitnessFunctionSquare((SquarePath) squarePath, robot);

		conf.setFitnessFunction( myFunc );
		
		Gene[] sampleGenes = new Gene[ 3 ];
		
		sampleGenes[0] = new IntegerGene(conf, minInt[0], maxInt[0] );  
		sampleGenes[1] = new IntegerGene(conf, minInt[1], maxInt[1] );  
		sampleGenes[2] = new IntegerGene(conf, minInt[2], maxInt[2] );  
		
		Chromosome sampleChromosome = new Chromosome(conf, sampleGenes );
		conf.setSampleChromosome( sampleChromosome );

		int populationSize = 5000;
//		int populationSize = 10000;

		conf.setPopulationSize( populationSize );
		Genotype population = Genotype.randomInitialGenotype(conf);

		
		long startTime = System.currentTimeMillis();
		
		List monitors = new Vector();
		monitors.add(new TimedMonitor(10));
		monitors.add(new FitnessImprovementMonitor(3, 3, 20d));
		IEvolutionMonitor monitor = new ChainedMonitors(monitors, 2);

				for( int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++ )
				  {
				      population.evolve();
				  }
				
				

//		List<String> messages = population.evolve(monitor);  
//		if (messages.size() > 0) {
//			for (String msg : messages) {
//				System.out.println("Message from monitor: " + msg+"\n");
//			}
//		}
				
		IChromosome bestSolutionSoFar = population.getFittestChromosome();
		long endTime = System.currentTimeMillis();
		System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
	    System.out.println("Total evolution time: " + (endTime - startTime) + " ms");
		
//		IChromosome bestSolutionSoFar = population.getFittestChromosome();
		System.out.println( "The best solution contained the following: " );

		System.out.println( MinimizingMakeChangeFitnessFunctionSquare.getEnergyAtGene(bestSolutionSoFar,(SquarePath) squarePath, robot) + " J totally required." );

		System.out.println("x = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 0 ));
		System.out.println("y = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 1 ));
		System.out.println("z = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 2 ));

		
		Set<EnergyPoint> hs = new HashSet<>();
		hs.addAll(eCloud);
		eCloud.clear();
		eCloud.addAll(hs);
		
		Map3D demo = new Map3D(xmlFilePath, eCloud);
	}
	
	public static void hybridGATest2() throws Exception {
	//		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";
//			String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST SW.xml";
			String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST.xml";
	
			double[] min = { -2.0, -2.0, 0.0 };
			double[] max = { 2.0, 2.0, 2.0 };
	
			int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
			int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
			
			WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);
			workingEnvelope.setResolution(0.1);
			
			double[][] h1 = { 
					{ 0.0, 0.0, 1.0, 0.0 }, 
					{ 0.0, 1.0, 0.0, 0.0 }, 
					{ -1.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 1.0 } };
			double[][] h2 = { 	
					{ 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 1.0, 0.0, -0.4 }, 
					{ -1.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 1.0 } };
			
			
	
			Target testTarget1 = new Target(h1);double[][] vertex1 = {
					{0.0, 0.0, 1.0, 0.0},
					{1.0, 0.0, 0.0, 0.0},
					{0.0, 1.0, 0.0, 0.0},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex2 = {
					{0.0, 0.0, 1.0, 0.0},
					{1.0, 0.0, 0.0, 0.0},
					{0.0, 1.0, 0.0, 0.4},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex3 = {
					{0.0, 0.0, 1.0, 0.0},
					{1.0, 0.0, 0.0, -0.4},
					{0.0, 1.0, 0.0, 0.4},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex4 = {
					{0.0, 0.0, 1.0, 0.0},
					{1.0, 0.0, 0.0, -0.4},
					{0.0, 1.0, 0.0, 0.0},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			
			Target point1 = new Target(vertex1);
			Target point2 = new Target(vertex2);
			Target point3 = new Target(vertex3);
			Target point4 = new Target(vertex4);
			
			
			
			Target testTarget2 = new Target(h2);
			Target[] targets = new Target[2];
			double[] backToCow = { - testTarget1.getPosition()[0],  - testTarget1.getPosition()[1],  - testTarget1.getPosition()[2]};
			testTarget1.translateTarget(backToCow);
			testTarget2.translateTarget(backToCow);
			double x_sample = 1E-3;
			double t_sample = 0.01;
			double acceleration = 1.0;
			double velocity = 0.2;
			Path p = new Path(testTarget1, testTarget2, t_sample, x_sample);
			p.setMaxAcc(acceleration);
			p.setMaxVel(velocity);
	
			Path squarePath = new SquarePath(point1,point2,point3,point4);
			squarePath.setMaxAcc(acceleration);
			squarePath.setMaxVel(velocity);
			squarePath.interpolate();
			
			Robot robot = new Robot(xmlFilePath, true);
	
			List<EnergyPoint> initialEnergyGrid = createEnergyCloud(robot, p , workingEnvelope);
			eCloud.addAll(initialEnergyGrid);
			Configuration conf = new DefaultConfiguration();
			
			System.out.println("The number of points displayed in eCloud at the moment is : " + eCloud.size());
			
			conf.setPreservFittestIndividual(true);
			conf.setKeepPopulationSizeConstant(false);
			
			
			FitnessFunction myFunc = new MinimizingMakeChangeFitnessFunction(testTarget1, testTarget2, robot);
	
			conf.setFitnessFunction( myFunc );
			
			Gene[] sampleGenes = new Gene[ 3 ];
			
			sampleGenes[0] = new IntegerGene(conf, minInt[0], maxInt[0] );  
			sampleGenes[1] = new IntegerGene(conf, minInt[1], maxInt[1] );  
			sampleGenes[2] = new IntegerGene(conf, minInt[2], maxInt[2] );  
			
			Chromosome sampleChromosome = new Chromosome(conf, sampleGenes );
			conf.setSampleChromosome( sampleChromosome );
	
			int populationSize = 5000;
	//		int populationSize = 10000;
	
			conf.setPopulationSize( populationSize );
			Population pop = initializePopulation(conf, sampleChromosome, initialEnergyGrid);
			
	//		Genotype population = Genotype.randomInitialGenotype(conf);
			Genotype population = new Genotype(conf, pop);
			population.fillPopulation(populationSize - pop.size());
	
			
			
			long startTime = System.currentTimeMillis();
			
			List monitors = new Vector();
			monitors.add(new TimedMonitor(30));
			monitors.add(new FitnessImprovementMonitor(1, 3, 0.1d));
			IEvolutionMonitor monitor = new ChainedMonitors(monitors, 2);
	
//					for( int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++ )
//					  {
//					      population.evolve();
//					  }
					
	//				IChromosome bestSolutionSoFar = population.getFittestChromosome();
	
			List<String> messages = population.evolve(monitor);  
			if (messages.size() > 0) {
				for (String msg : messages) {
					System.out.println("Message from monitor: " + msg+"\n");
				}
			}
	
			long endTime = System.currentTimeMillis();
			System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
		    System.out.println("Total evolution time: " + (endTime - startTime) + " ms");
			
			IChromosome bestSolutionSoFar = population.getFittestChromosome();
			System.out.println( "The best solution contained the following: " );
	
			System.out.println( MinimizingMakeChangeFitnessFunction.getEnergyAtGene(bestSolutionSoFar, testTarget1, testTarget2, robot) + " J totally required." );
	
			System.out.println("x = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 0 ));
			System.out.println("y = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 1 ));
			System.out.println("z = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 2 ));
	
			
	//		Set<EnergyPoint> hs = new HashSet<>();
	//		hs.addAll(eCloud);
	//		eCloud.clear();
	//		eCloud.addAll(hs);
			
			Map3D demo = new Map3D(xmlFilePath, eCloud);
		}


	public static void hybridGATestSquare() throws Exception {
			String xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST SW.xml";
	
			double[] min = { -2.0, -2.0, 0.0 };
			double[] max = { 2.0, 2.0, 2.0 };
	
			int[] minInt = {(int) min[0]*100, (int) min[1]*100, (int) min[2]*100};
			int[] maxInt = {(int) max[0]*100, (int) max[1]*100, (int) max[2]*100};
			
			WorkingEnvelope workingEnvelope = new WorkingEnvelope(min, max);
			workingEnvelope.setResolution(0.1);
			
			
	
			double[][] vertex1 = {
					{0.0, 0.0, 1.0, 0.0},
					{0.0, 1.0, 0.0, 0.0},
					{-1.0, 0.0, 0.0, 0.0},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex2 = {
					{0.0, 0.0, 1.0, 0.0},
					{0.0, 1.0, 0.0, 0.0},
					{-1.0, 0.0, 0.0, 0.4},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex3 = {
					{0.0, 0.0, 1.0, 0.0},
					{0.0, 1.0, 0.0, -0.4},
					{-1.0, 0.0, 0.0, 0.4},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			double[][] vertex4 = {
					{0.0, 0.0, 1.0, 0.0},
					{0.0, 1.0, 0.0, -0.4},
					{-1.0, 0.0, 0.0, 0.0},
					{0.0, 0.0, 0.0, 1.0}
			};
			
			
	//		double[][] vertex1 = {
	//				{0.0, 0.0, 1.0, 0.59627},
	//				{1.0, 0.0, 0.0, -0.21132},
	//				{0.0, 1.0, 0.0, 0.38041},
	//				{0.0, 0.0, 0.0, 1.0}
	//		};
	//		
	//		double[][] vertex2 = {
	//				{0.0, 0.0, 1.0, 0.49627},
	//				{1.0, 0.0, 0.0, -0.21132},
	//				{0.0, 1.0, 0.0, 0.38041},
	//				{0.0, 0.0, 0.0, 1.0}
	//		};
	//		
	//		double[][] vertex3 = {
	//				{0.0, 0.0, 1.0, 0.49627},
	//				{1.0, 0.0, 0.0, -0.31132},
	//				{0.0, 1.0, 0.0, 0.38041},
	//				{0.0, 0.0, 0.0, 1.0}
	//		};
	//		
	//		double[][] vertex4 = {
	//				{0.0, 0.0, 1.0, 0.59627},
	//				{1.0, 0.0, 0.0, -0.31132},
	//				{0.0, 1.0, 0.0, 0.38041},
	//				{0.0, 0.0, 0.0, 1.0}
	//		};
			
			Target point1 = new Target(vertex1);
			Target point2 = new Target(vertex2);
			Target point3 = new Target(vertex3);
			Target point4 = new Target(vertex4);
			
			
		
			double x_sample = 1E-3;
			double t_sample = 0.01;
			double acceleration = 1.0;
			double velocity = 0.2;
	
			Path squarePath = new SquarePath(point1,point2,point3,point4);
			squarePath.setMaxAcc(acceleration);
			squarePath.setMaxVel(velocity);
			squarePath.interpolate();
			
			Robot robot = new Robot(xmlFilePath, true);
			
			
			List<EnergyPoint> initialEnergyGrid = createEnergyCloudGeneral(robot, squarePath , workingEnvelope);
			eCloud.addAll(initialEnergyGrid);
			Configuration conf = new DefaultConfiguration();
			
			conf.setPreservFittestIndividual(true);
			conf.setKeepPopulationSizeConstant(false);
			
			
			FitnessFunction myFunc = new MinimizingMakeChangeFitnessFunctionSquare((SquarePath) squarePath, robot);
	
			conf.setFitnessFunction( myFunc );
			
			Gene[] sampleGenes = new Gene[ 3 ];
			
			sampleGenes[0] = new IntegerGene(conf, minInt[0], maxInt[0] );  
			sampleGenes[1] = new IntegerGene(conf, minInt[1], maxInt[1] );  
			sampleGenes[2] = new IntegerGene(conf, minInt[2], maxInt[2] );  
			
			Chromosome sampleChromosome = new Chromosome(conf, sampleGenes );
			conf.setSampleChromosome( sampleChromosome );
	
			int populationSize = 5000;
	//		int populationSize = 10000;
	
			conf.setPopulationSize( populationSize );
			Population pop = initializePopulation(conf, sampleChromosome, initialEnergyGrid);
	//		
			Genotype population = new Genotype(conf, pop);
	//		Genotype population = Genotype.randomInitialGenotype(conf);
			population.fillPopulation(populationSize - pop.size());
	
			
			long startTime = System.currentTimeMillis();
			
			List monitors = new Vector();
			monitors.add(new TimedMonitor(10));
			monitors.add(new FitnessImprovementMonitor(3, 3, 20d));
			IEvolutionMonitor monitor = new ChainedMonitors(monitors, 2);
	
					for( int i = 0; i < MAX_ALLOWED_EVOLUTIONS; i++ )
					  {
					      population.evolve();
					  }
					
					
	
	//		List<String> messages = population.evolve(monitor);  
	//		if (messages.size() > 0) {
	//			for (String msg : messages) {
	//				System.out.println("Message from monitor: " + msg+"\n");
	//			}
	//		}
					
			IChromosome bestSolutionSoFar = population.getFittestChromosome();
			long endTime = System.currentTimeMillis();
			System.out.println("The number of points displayed is eCloud size: " + eCloud.size());
		    System.out.println("Total evolution time: " + (endTime - startTime) + " ms");
			
	//		IChromosome bestSolutionSoFar = population.getFittestChromosome();
			System.out.println( "The best solution contained the following: " );
	
			System.out.println( MinimizingMakeChangeFitnessFunctionSquare.getEnergyAtGene(bestSolutionSoFar,(SquarePath) squarePath, robot) + " J totally required." );
	
			System.out.println("x = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 0 ));
			System.out.println("y = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 1 ));
			System.out.println("z = " + MinimizingMakeChangeFitnessFunction.getPositionDouble(bestSolutionSoFar, 2 ));
	
			
			Map3D demo = new Map3D(xmlFilePath, eCloud);
		}


	private static List<EnergyPoint> createEnergyCloud(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */
		ExecutorService execs = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		long time = System.currentTimeMillis();
		
		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];

		double weResolution = we.getResolution();
		double[] weMinValues = we.getMinValues();
		long[] weSize = we.getSize();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

		List<Future<EnergyPoint>> energyCloud = new LinkedList<Future<EnergyPoint>>();
		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();
		
		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		int[] index = new int[3];

		for (index[0] = 0; index[0] < weSize[0]; index[0]++){
			for (index[1] = 0; index[1] < weSize[1]; index[1]++){
				for (index[2] = 0; index[2] < weSize[2]; index[2]++){

					for(int j = 0; j < 3; j++){					
						buffer[j] = new BigDecimal(weMinValues[j] + index[j] * weResolution, mt);
					}
					Target[] curr = new Target[2];
					curr[0] = new Target(inPosH);
					curr[1] = new Target(finPosH);
					
					Point3D point = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());
					
					currPos = Matrix.copyVector(point.getPosition());
					curr[0].translateTarget(currPos);
					curr[1].translateTarget(currPos);

					PlannerList thread = new PlannerList(curr[0], curr[1], r, point);

					Future<EnergyPoint> result = execs.submit(thread);

					energyCloud.add(result);
				}
			}
		}
		
		execs.shutdown();
		
		for(Future<EnergyPoint> point : energyCloud){
			try {
				if(point.get().getEnergy() != 0.0){
					energyCloudFinal.add(point.get());
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				e.getCause();
			}
		}
		
		time = System.currentTimeMillis() - time;
		System.out.println("Time for calculating energy values: " + time/1000 + " sec");
		
		return energyCloudFinal;

	}
	
	private static List<EnergyPoint> createEnergyCloudGeneral(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */
		ExecutorService execs = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		long time = System.currentTimeMillis();
		
		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];
		
		Target point1 = null;
		Target point2 = null;
		Target point3 = null;
		Target point4 = null;
		boolean isSquare = false;
		
		double weResolution = we.getResolution();
		double[] weMinValues = we.getMinValues();
		long[] weSize = we.getSize();

		Target[] targets = p.getPathPositions();

		if(p instanceof SquarePath){
			point1 = new Target (targets[0].getHomMatrix());
			point2 = new Target (targets[1].getHomMatrix());
			point3 = new Target (targets[2].getHomMatrix());
			point4 = new Target (targets[3].getHomMatrix());
			isSquare = true;
			
		}else{
			inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
			finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());
		}
		
		

		List<Future<EnergyPoint>> energyCloud = new LinkedList<Future<EnergyPoint>>();
		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();
		
		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		int[] index = new int[3];

		for (index[0] = 0; index[0] < weSize[0]; index[0]++){
			for (index[1] = 0; index[1] < weSize[1]; index[1]++){
				for (index[2] = 0; index[2] < weSize[2]; index[2]++){

					for(int j = 0; j < 3; j++){					
						buffer[j] = new BigDecimal(weMinValues[j] + index[j] * weResolution, mt);
					}
					
					Target[] curr = null;
					Point3D point = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());
					currPos = Matrix.copyVector(point.getPosition());
					Path newPath = null;
					
					if(isSquare){
						Target point1new = new Target(point1.getHomMatrix());
						Target point2new = new Target(point2.getHomMatrix());
						Target point3new = new Target(point3.getHomMatrix());
						Target point4new = new Target(point4.getHomMatrix());
						
						
						point1new.translateTarget(currPos);
						point2new.translateTarget(currPos);
						point3new.translateTarget(currPos);
						point4new.translateTarget(currPos);
						
						newPath = new SquarePath(point1new,point2new,point3new,point4new);
						newPath.setMaxAcc(p.getMaxAcc());
						newPath.setMaxVel(p.getMaxVel());
						
					}else{
						curr = new Target[2];
						curr[0] = new Target(inPosH);
						curr[1] = new Target(finPosH);

						

						
						curr[0].translateTarget(currPos);
						curr[1].translateTarget(currPos);
						
						newPath = new Path(curr[0],curr[1]);
						newPath.setMaxAcc(p.getMaxAcc());
						newPath.setMaxVel(p.getMaxVel());

					}
					
					PlannerList thread = new PlannerList(newPath, r, point);

					Future<EnergyPoint> result = execs.submit(thread);

					energyCloud.add(result);
				}
			}
		}
		
		execs.shutdown();
		
		for(Future<EnergyPoint> point : energyCloud){
			try {
				if(point.get().getEnergy() != 0.0){
					energyCloudFinal.add(point.get());
				}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				e.getCause();
			}
		}
		
		time = System.currentTimeMillis() - time;
		System.out.println("Time for calculating energy values: " + time/1000 + " sec");
		
		return energyCloudFinal;

	}
	
	private static ChromosomePool initializePool(ChromosomePool pool, Configuration conf,Chromosome sampleChromosome,List<EnergyPoint> energyList) throws InvalidConfigurationException{
//		ChromosomePool pool = new ChromosomePool();
		
		for(EnergyPoint point : energyList){
			Gene[] newGenes = new Gene[ 3 ];
			
			newGenes = sampleChromosome.getGenes();
			
			newGenes[0].setAllele((int) point.getX() * 100);
			newGenes[1].setAllele((int) point.getY() * 100);
			newGenes[2].setAllele((int) point.getZ() * 100);
			
			IChromosome chrom = Chromosome.randomInitialChromosome(conf);
			chrom.setGenes(newGenes);
			
			pool.releaseChromosome(chrom);
		}
		
		return pool;
	}
	
	private static Population initializePopulation(Configuration conf,Chromosome sampleChromosome,List<EnergyPoint> energyList) throws InvalidConfigurationException{
		Population pop = new Population(conf, conf.getPopulationSize());
		
		for(EnergyPoint point : energyList){
			Gene[] newGenes = new Gene[ 3 ];
			
			newGenes = sampleChromosome.getGenes();
			
			newGenes[0].setAllele((int) point.getX() * 100);
			newGenes[1].setAllele((int) point.getY() * 100);
			newGenes[2].setAllele((int) point.getZ() * 100);
			
			IChromosome chrom = Chromosome.randomInitialChromosome(conf);
			chrom.setGenes(newGenes);
			
			pop.addChromosome(chrom);
		}
		
		return pop;
	}
	
	private static EnergyPoint findMinimum(List<EnergyPoint> points){
		EnergyPoint minimum = null;
		double energy = Double.MAX_VALUE;
		
		for(EnergyPoint point : points){
			if(point.getEnergy() < energy){
				energy = point.getEnergy();
				minimum = point;
			}
		}
		
		
		return minimum;
	}
	
}
