import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import javax.media.j3d.*;
import javax.vecmath.*;
import javax.swing.*;
import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.universe.SimpleUniverse;
import com.sun.j3d.utils.geometry.Primitive;
import com.sun.j3d.utils.geometry.Sphere;
import com.sun.j3d.utils.picking.PickCanvas;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickTool;
import com.sun.j3d.loaders.objectfile.ObjectFile;
import com.sun.j3d.utils.behaviors.vp.OrbitBehavior;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import com.sun.j3d.loaders.Scene;

/*
 * INPUTS
 * - xmlFilePath: a String defining the path of the XML file that describes the robot structure;
 * - p: a Path object defining the interpolated trajectory we want the robot to perform;
 * - we: a WorkingEnvelope object defining the volume in which we want the algorithm to consider the initial position of the trajectory.
 * 
 * OUTPUT
 * 3D visualization to help the user deciding where it's best(from an energy consumption point of view)
 * to set the initial position of the predefined trajectory
 */

public class Map3D extends MouseAdapter{


	private SimpleUniverse u = null;
	private Canvas3D canvas = null;
	private PickCanvas pickCanvas;
	private List<EnergyPoint> energyCloud;

	public Map3D(String xmlFilePath, Path p, WorkingEnvelope we) {

		Robot r = new Robot(xmlFilePath, true);

		/*
		 * Here I create the frame and its basic elements
		 */
		Frame frame = new Frame("Map3D");
		GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
		canvas = new Canvas3D(config);
		canvas.setSize(400, 400);

		/*
		 * Here I create the universe to fill in with the objects that I'm going to add calling the method createSceneGraph 
		 */
		u = new SimpleUniverse(canvas);
		BranchGroup scene = createSceneGraph(r, p, we);
		u.addBranchGraph(scene);

		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView());

		frame.addWindowListener(new WindowAdapter() {

			public void windowClosing(WindowEvent winEvent) {

				System.exit(0);

			}

		});

		frame.add(canvas);

		pickCanvas = new PickCanvas(canvas, scene);

		pickCanvas.setMode(PickCanvas.GEOMETRY);

		canvas.addMouseListener(this);

		frame.pack();

		frame.show();

	}

	public void mouseClicked(MouseEvent e){

		pickCanvas.setShapeLocation(e);

		PickResult result = pickCanvas.pickClosest();

		if (result == null) {

			System.out.println("Nothing picked");

		} else {

			Primitive sphere = (Primitive)result.getNode(PickResult.PRIMITIVE);

			Shape3D robot = (Shape3D)result.getNode(PickResult.SHAPE3D);

			if (sphere != null) {
				
				energyCloud.sort(null);
				double energyMin = energyCloud.get(0).getEnergy();
				double energy = 0d;
				double energyWaste = 0d;
				double[] h = new double[16];
				Transform3D t3d = new Transform3D();
				sphere.getLocalToVworld(t3d);
				t3d.get(h);

				for(EnergyPoint ep : energyCloud){

					if( ep.getPosition()[0] == h[3] && ep.getPosition()[1] == h[7] && ep.getPosition()[2] == h[11]){

						energy = ep.getEnergy();
						energyWaste = Math.round(energy/energyMin);

					}
				}

				System.out.println("");
				System.out.println("CURRENT INITIAL POSITION");
				System.out.println("x = " + h[3] + " y = " + h[7] + " z = " + h[11]);
				System.out.println("CURRENT ENERGY CONSUMPTION: " + energy + " J");
				System.out.println("That is approximately " + energyWaste + " times more consuming of the minimum value of energy consumption.");
				System.out.println("INITIAL POSITION OF MIN ENERGY VALUE");
				System.out.println("x = " + energyCloud.get(0).getPosition()[0] + " y = " + energyCloud.get(0).getPosition()[1] + " z = " + energyCloud.get(0).getPosition()[2]);
				System.out.println("MIN VALUE: " + energyMin + " J");
				System.out.println("INITIAL POSITION OF MIN2 ENERGY VALUE");
				System.out.println("x = " + energyCloud.get(1).getPosition()[0] + " y = " + energyCloud.get(1).getPosition()[1] + " z = " + energyCloud.get(1).getPosition()[2]);
				System.out.println("MIN VALUE: " + energyCloud.get(1).getEnergy() + " J");
				System.out.println("INITIAL POSITION OF MIN3 ENERGY VALUE");
				System.out.println("x = " + energyCloud.get(2).getPosition()[0] + " y = " + energyCloud.get(2).getPosition()[1] + " z = " + energyCloud.get(2).getPosition()[2]);
				System.out.println("MIN VALUE: " + energyCloud.get(2).getEnergy() + " J");

			} else if (robot != null) {

				System.out.println(robot.getClass().getName());

			} else{

				System.out.println("null");

			}

		}

	}

	public Map3D(String xmlFilePath, List<EnergyPoint> energyList){
		double distanceView = 110;
		double[] centerView = { -0.1, 0.45, 3.5 };

		Robot r = new Robot(xmlFilePath, true);

		Frame frame = new Frame("Map3D");

		GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
		canvas = new Canvas3D(config);
		canvas.setSize(400, 400);

		u = new SimpleUniverse(canvas);




		BranchGroup scene = createSceneGraph(energyList);
		System.out.println("End creating SceneGraph");

		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView());

		u.addBranchGraph(scene);

		frame.addWindowListener(new WindowAdapter() {

			public void windowClosing(WindowEvent winEvent) {

				System.exit(0);

			}

		});

		frame.add(canvas);

		pickCanvas = new PickCanvas(canvas, scene);

		pickCanvas.setMode(PickCanvas.GEOMETRY);

		canvas.addMouseListener(this);

		frame.pack();

		frame.show();

		System.out.println("End of Map3D");
	}

//	public static void main(String[] args) {
//
//		//		JProgressBar progressBar = new JProgressBar(0, 100);
//		//		progressBar.setValue(0);
//		//        progressBar.setStringPainted(true);
//		//		JApplet applet = new JApplet();
//		//		applet.add(progressBar);
//
//		/*
//		 * TEST
//		 */
//
//		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";
//
//		//		double[][] h1 = { { 0.0, 0.0, 1.0, 0.59627 }, 
//		//				{ 1.0, 0.0, 0.0, 0.0 }, 
//		//				{ 0.0, 1.0, 0.0, 0.66282 },
//		//				{ 0.0, 0.0, 0.0, 1.0 } };
//		//		double[][] h2 = { { 0.0, 0.0, 1.0, 0.59627 },
//		//				{ 1.0, 0.0, 0.0, -0.21132 }, 
//		//				{ 0.0, 1.0, 0.0, 0.58041 },
//		//				{ 0.0, 0.0, 0.0, 1.0 } };
//		double[][] h1 = { 	{ 1.0, 0.0, 0.0, 0.0 }, 
//				{ 0.0, 1.0, 0.0, 0.2 }, 
//				{ 0.0, 0.0, -1.0, 0.0 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
//		double[][] h2 = 	{ 	{ 1.0, 0.0, 0.0, 0.0 },
//				{ 0.0, 1.0, 0.0, -0.2 }, 
//				{ 0.0, 0.0, -1.0, 0.0 },
//				{ 0.0, 0.0, 0.0, 1.0 } };
//
//		//		Quat4d q = new Quat4d(0.28071, -0.46957, -0.82150, -0.16078);
//		//		double[][] rot = Path.quatToRotm(q);
//		//		double[] pos0 = {0d, 0d, 0d};
//		//		double[] pos1 = {0d, 1.9, 340.7};
//		//		double[] pos2 = {-5.5, -633.1, 337.3};
//		//		double[] pos3 = {-5.5, -633.1, 0d};
//		//		Target testTarget1 = new Target(pos0, rot);
//		//		Target testTarget2 = new Target(pos1, rot);
//		//		Target testTarget3 = new Target(pos2, rot);
//		//		Target testTarget4 = new Target(pos3, rot);
//
//		Target testTarget1 = new Target(h1);
//		Target testTarget2 = new Target(h2);
//		//		double[] backToCow = { - testTarget1.getPosition()[0],  - testTarget1.getPosition()[1],  - testTarget1.getPosition()[2]};
//		//		testTarget1.translateTarget(backToCow);
//		//		testTarget2.translateTarget(backToCow);
//		double x_sample = 1E-3;
//		double t_sample = 0.01;
//		double acceleration = 1.0;
//		double velocity = 0.5;
//		Path p = new Path(testTarget1, testTarget2, t_sample, x_sample);
//		//		Path p = new SquarePath(testTarget1, testTarget2, testTarget3, testTarget4, t_sample, x_sample);
//
//		p.setMaxAcc(acceleration);
//		p.setMaxVel(velocity);
//
//		//		double[] min = {-0.5, 0.0, 0.0};
//		//		double[] max = {1.0, 1.0, 1.0};
//		//		double res = 0.05;
//		double[] min = { -1.0, -1.0, 0.0 };
//		double[] max = { 1.0, 1.0, 1.0 };
//		double res = 0.05;
//
//		long time = System.currentTimeMillis();
//
//		WorkingEnvelope we = new WorkingEnvelope(min, max, res);
//
//		//		new Map3d (xmlFilePath, p, we);
//
//		new Map3D(xmlFilePath, p, we);
//
//		//		Matrix.displayVector(backToCow);
//		//		Matrix.displayMatrix(testTarget1.getHomMatrix());
//
//		time = System.currentTimeMillis() - time;
//
//		System.out.println("Overall time: " + time/1000 + " sec");
//
//
//	}

	private BranchGroup createSceneGraph(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It adds all the objects to the world: 
		 * - energy point cloud; 
		 * - robot model;
		 * - reference frame.
		 */

		energyCloud = createEnergyCloudThreadedNEWNEW(r, p, we);
		//		energyCloud = createEnergyCloudThreaded(r, p, we);
//		energyCloud = createEnergyCloudThreadedLoopNEW(r, p, we);

		int numFeasibleSol = energyCloud.size();
		long[] weSize = we.getSize();
		long weSizeTot = weSize[0] * weSize[1] * weSize[2];

		BranchGroup objRoot = new BranchGroup();
		TransformGroup cow = new TransformGroup();
		cow.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		objRoot.addChild(cow);

		/*
		 * Mouse interaction
		 */
		OrbitBehavior ob = new OrbitBehavior(u.getCanvas());
		ob.setReverseRotate(true);
		ob.setReverseTranslate(true);
		ob.setSchedulingBounds(new BoundingSphere(new Point3d(0.0,0.0,0.0),Double.MAX_VALUE));
		u.getViewingPlatform().setViewPlatformBehavior(ob);

		int numClusters = 5;
		computeEnergyColors(energyCloud, numClusters);

		for(EnergyPoint ep : energyCloud)
			cow.addChild(energyPoint(ep));

		/*
		 * Here I call the method robotModel to add the robot to the branch
		 */

		cow.addChild(robotModel());

		cow.addChild(axis());

		/*
		 * Finally, I set the Lights so that the robot can be seen from any angle
		 */

		objRoot.addChild(setDirectionalLight(new Vector3f(1.0f, 1.0f, 1.0f)));
		objRoot.addChild(setDirectionalLight(new Vector3f(-1.0f, 1.0f, -1.0f)));
		objRoot.addChild(setDirectionalLight(new Vector3f(-1.0f, -1.0f, 1.0f)));

		System.out.println("Number of feasible solutions: " + numFeasibleSol + "/" + weSizeTot);

		return objRoot;

	}

	private BranchGroup createSceneGraph( List<EnergyPoint> energyCloud) {
		int numFeasibleSol = energyCloud.size();
		int numPoint = 0;

		float[] energyColor = new float[3];

		double[] CenterOfWorld = {0.0, 0.0, 0.0};
		double[] currPos = new double[3];

		BranchGroup objRoot = new BranchGroup();
		TransformGroup cow = new TransformGroup();
		TransformGroup rotObj = new TransformGroup();
		Transform3D rot = new Transform3D();
		rot.rotX(-Math.PI/2);
		rotObj.setTransform(rot);
		cow.addChild(rotObj);
		cow.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		objRoot.addChild(cow);


		OrbitBehavior ob = new OrbitBehavior(u.getCanvas());
		ob.setReverseRotate(true);
		ob.setReverseTranslate(true);
		ob.setSchedulingBounds(new BoundingSphere(new Point3d(0.0,0.0,0.0),Double.MAX_VALUE));
		u.getViewingPlatform().setViewPlatformBehavior(ob);

		int numClusters = 5;
		computeEnergyColors(energyCloud, numClusters);

		for(EnergyPoint ep : energyCloud)
			rotObj.addChild(energyPoint(ep));

		/*
		 * Here I call the method robotModel to add the robot to the branch
		 */

		rotObj.addChild(robotModel());

		rotObj.addChild(axis());

		/*
		 * Finally, I set the Lights
		 */

		objRoot.addChild(setDirectionalLight(new Vector3f(1.0f, 1.0f, 1.0f)));
		objRoot.addChild(setDirectionalLight(new Vector3f(-1.0f, 1.0f, -1.0f)));
		objRoot.addChild(setDirectionalLight(new Vector3f(-1.0f, -1.0f, 1.0f)));

		System.out.println("Number of feasible solutions: " + numFeasibleSol );

		return objRoot;
	}

	private class Result{

		int numClusters;
		int[] cluster;
		int[] size;
		double[] centers;
		double[] withinss;

		private Result(int nc, int dimCluster, int dimSize, int dimCenters, int dimWithinss){

			numClusters = nc;
			cluster = new int[dimCluster];
			size = new int[dimSize];
			centers = new double[dimCenters];
			withinss = new double[dimWithinss];

		}

	}

	private void computeEnergyColors(List<EnergyPoint> energyCloud, int numClusters) {

		long time = System.currentTimeMillis();

		int[][] B = new int[numClusters + 1][energyCloud.size() + 1];

		double[][] D = new double[numClusters + 1][energyCloud.size() + 1];

		Result r;

		EnergyPoint[] energyCloudSorted = new EnergyPoint[energyCloud.size() + 1];
		EnergyPoint[] energyCloudArray = new EnergyPoint[energyCloud.size()];
		for(int i = 0; i < energyCloud.size(); i++)
			energyCloudArray[i] = energyCloud.get(i);

		System.arraycopy(energyCloudArray, 0, energyCloudSorted, 1, energyCloudArray.length);
		Arrays.sort(energyCloudSorted, 1, energyCloudSorted.length);

		fillDpMatrix(energyCloudSorted, D, B);

		r = backtrack(energyCloudSorted, B);

		// Perform clustering on the original data
		for(int i = 1; i < energyCloud.size(); ++i) {

			int indexLeft = 1;
			int indexRight;

			for (int k = 1; k < r.size.length; ++k) {

				indexRight = indexLeft + r.size[k] - 1;

				if ( energyCloud.get(i).getEnergy() <= energyCloudSorted[indexRight].getEnergy() ) {

					r.cluster[i] = k;
					break;

				}

				indexLeft = indexRight + 1;

			}

		}

		System.out.println("Number of clusters: " + r.numClusters);
		double[] dimClusters = new double[r.size.length - 1];
		for(int i = 0; i < dimClusters.length; i++)
			dimClusters[i] = (double) r.size[i + 1];
		System.out.println("Dimension of clusters:");
		Matrix.displayVector(dimClusters);

		/*
		 * Choosing color for each cluster
		 */

		for(int i = 0; i < energyCloud.size(); i++){
			//        	
			float[] rgb = {0f, 1f, 0f};

			if(r.cluster[i + 1] > Math.ceil(numClusters/2)){
				rgb[0] = 1f;
				rgb[1] -= (r.cluster[i + 1] - Math.ceil(numClusters/2))*0.5f;
			}else
				rgb[0] += r.cluster[i + 1]*0.5f;


			Color3f color = new Color3f(rgb);

			energyCloud.get(i).setColor(color);

		}

		time = System.currentTimeMillis() - time;

		System.out.println("Time for computing clustering: " + time/1000 + " sec");


	}

	private Result backtrack(EnergyPoint[] x, int[][] B) {

		int K = B.length - 1;
		int N = B[0].length - 1;
		int cluster_right = N;
		int cluster_left;

		Result r = new Result(K, N + 1, K + 1, N + 1, K + 1);

		// Backtrack the clusters from the dynamic programming matrix
		for(int k = K; k >= 1; --k) {
			cluster_left = B[k][cluster_right];

			for(int i = cluster_left; i <= cluster_right; ++i)
				r.cluster[i] = k;

			double sum = 0d;

			for(int i = cluster_left; i <= cluster_right; ++i)
				sum += x[i].getEnergy();

			r.centers[k] = sum/(cluster_right - cluster_left+1);

			for(int i = cluster_left; i <= cluster_right; ++i)
				r.withinss[k] += (x[i].getEnergy() - r.centers[k]) *
				(x[i].getEnergy() - r.centers[k]);

			r.size[k] = cluster_right - cluster_left + 1;

			if(k > 1) {
				cluster_right = cluster_left - 1;
			}
		}

		return r;

	}	


	private void fillDpMatrix(EnergyPoint[] x, double[][] D, int[][] B) {
		int K = D.length - 1;
		int N = D[0].length -1;
		double d;

		for (int i = 1; i <= K; ++i) {
			D[i][1] = 0d;
			B[i][1] = 1;
		}

		double meanX1;
		double meanXj;

		for (int k = 1; k <= K; ++k) {
			meanX1 = x[1].getEnergy();

			for (int i = Math.max(2, k); i <= N; ++i) {
				if (k == 1) {
					D[1][i] = D[1][i-1] + (i-1) / (double) i *
							(x[i].getEnergy() - meanX1) * (x[i].getEnergy() - meanX1);
					meanX1 = ((i - 1) * meanX1 + x[i].getEnergy()) / (double)i;

					B[1][i] = 1;
				} else {
					D[k][i] = -1d;
					d = 0d;
					meanXj = 0d;

					for (int j = i; j >= k; --j) {
						d += (i - j) / (double) (i - j + 1) * 
								(x[j].getEnergy() - meanXj) * (x[j].getEnergy() - meanXj);
						meanXj = (x[j].getEnergy() + (i - j) * meanXj) / (double)(i - j + 1);

						if (D[k][i] == -1d) {
							if(j == 1) {
								D[k][i] = d;
								B[k][i] = j;
							} else {
								D[k][i] = d + D[k - 1][j - 1];
								B[k][i] = j;
							}
						} else {
							if(j == 1) {
								if(d <= D[k][i]) {
									D[k][i] = d;
									B[k][i] = j;
								}
							} else {
								if(d + D[k - 1][j - 1] < D[k][i]) {
									D[k][i] = d + D[k - 1][j - 1];
									B[k][i] = j;
								}
							}
						}
					}
				}

			}

		}

	}

	private TransformGroup energyPoint(EnergyPoint ep) {

		TransformGroup tg = setPosition(ep.getPosition());


		//		Sphere s = new Sphere(0.04f);
		Sphere s = new Sphere(0.14f);
		/*
		 * Setting the appearance of the sphere (color and transparency)
		 */

		Appearance app = new Appearance();
		Color3f color = ep.getColor();
		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.FASTEST);
		app.setColoringAttributes(ca);

		s.setCapability(Node.ENABLE_PICK_REPORTING);

		//		PickTool.setCapabilities(s, PickTool.INTERSECT_FULL);

		TransparencyAttributes ta = new TransparencyAttributes();
		ta.setTransparencyMode(TransparencyAttributes.BLENDED);
		ta.setTransparency(0.15f);
		app.setTransparencyAttributes(ta);

		s.setAppearance(app);

		tg.addChild(s);

		return tg;

	}

	private List<EnergyPoint> createEnergyCloud(Robot r, Path p, WorkingEnvelope we) {

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

		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */

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

//					targets[0].setHomMatrix(inPosH);
//					targets[1].setHomMatrix(finPosH);

					Point3D point = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());
					
					currPos = Matrix.copyVector(point.getPosition());
					curr[0].translateTarget(currPos);
					curr[1].translateTarget(currPos);

					/*
					 * Here I create the thread to calculate the energy value
					 * with current initial and final position
					 */
					PlannerList thread = new PlannerList(curr[0], curr[1], r, point);

					Future<EnergyPoint> result = execs.submit(thread);

					energyCloud.add(result);

					/*
					 * Every time the Online Planner calculates an energy value
					 * it stores it in energyList. In order to get the current
					 * energy value I have to "ask" for the last element added
					 * to the list
					 */

					//					System.out.println(
					//							"energyMatrix[" + x + "]" + "[" + y + "]" + "[" + z + "]" + " = " + energyMatrix[x][y][z]);

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

	private List<EnergyPoint> createEnergyCloudThreaded(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		long time = System.currentTimeMillis();

		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];

		List<Point3D> weList = we.getWeList();
		List<Future<EnergyPoint>> energyCloud = new LinkedList<Future<EnergyPoint>>();
		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());


		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */
		ExecutorService execs = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

		for(Point3D point : weList){

			Target[] curr = new Target[2];
			curr[0] = new Target(inPosH);
			curr[1] = new Target(finPosH);

			targets[0].setHomMatrix(inPosH);
			targets[1].setHomMatrix(finPosH);

			currPos = Matrix.copyVector(point.getPosition());

			//					System.out.println("currPos = ");
			//					Matrix.displayVector(currPos);

			curr[0].translateTarget(currPos);
			curr[1].translateTarget(currPos);

			//					System.out.println("Current initial position: ");
			//					Matrix.displayMatrix(curr[0].getHomMatrix());

			//					System.out.println("Current final position: ");
			//					Matrix.displayMatrix(curr[1].getHomMatrix());

			/*
			 * Here I create the thread to calculate the energy value
			 * with current initial and final position
			 */

			//					dynSolution.run(curr[0], curr[1], point);
			//			OnlinePlannerCallable thread = new OnlinePlannerCallable(curr[0], curr[1], r, point);
			PlannerList thread = new PlannerList(curr[0], curr[1], r, point);

			Future<EnergyPoint> result = execs.submit(thread);

			energyCloud.add(result);

			/*
			 * Every time the Online Planner calculates an energy value
			 * it stores it in energyList. In order to get the current
			 * energy value I have to "ask" for the last element added
			 * to the list
			 */
			//					double energy = dynSolution.getEnergyList().get(dynSolution.getEnergyList().size() - 1).getEnergy();
			//					
			//					if(energy != 0.0){
			//						
			//					EnergyPoint ep = new EnergyPoint();
			//					
			//					ep.setPosition(point.getPosition());
			//					ep.setEnergy(energy);
			//					
			//					energyCloudFinal.add(ep);
			//					
			//					}

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

	private List<EnergyPoint> createEnergyCloudThreadedLoop(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		/*
		 * PROGRESS BAR
		 */
		JFrame f = new JFrame("Processing...");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container content = f.getContentPane();
		JProgressBar progressBar = new JProgressBar(0, 100);
		progressBar.setValue(0);
		progressBar.setStringPainted(true);
		content.add(progressBar, BorderLayout.CENTER);
		f.setSize(400, 100);
		f.setVisible(true);

		long time = System.currentTimeMillis();

		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];

		List<Point3D> weList = we.getWeList();

		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */
		int numberOfThreads = Runtime.getRuntime().availableProcessors();
		int portion = weList.size()/numberOfThreads;
		int progress = 0;
		int dProgress = 100/numberOfThreads;
		List<Thread> threads = new LinkedList<Thread>();
		List<OnlinePlanner> planners = new LinkedList<OnlinePlanner>();
		//		List<PlannerList> planners = new LinkedList<PlannerList>();

		for(int i = 0; i < numberOfThreads; i++){
			int from = i*portion;
			int to = i < numberOfThreads-1 ? (i+1)*portion : weList.size();

			OnlinePlanner planner = new OnlinePlanner(p, r, weList.subList(from, to-1));
			//			PlannerList planner = new PlannerList(p, r, weList.subList(from, to-1));
			Thread thread = new Thread( planner );
			planners.add(planner);
			threads.add(thread);
		}

		for(Thread currThread : threads){
			currThread.start();
		}

		for(Thread currThread : threads){
			try {
				currThread.join();
				if(currThread.getState() == Thread.State.TERMINATED){
					System.out.println("Thread is terminated.");
					progress += dProgress;			
					progressBar.setValue(progress);
					progressBar.setStringPainted(true);
					progressBar.paint(progressBar.getGraphics());
					if(progress == 100)
						f.setVisible(false);
				}
			}
			catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		for(OnlinePlanner elem : planners){
			for(EnergyPoint point : elem.getEnergyList()){
				if(point.getEnergy() != 0.0){
					energyCloudFinal.add(point);
				}
			}
		}

		//		for(PlannerList elem : planners){
		//			for(EnergyPoint point : elem.getEnergyList()){
		//				if(point.getEnergy() != 0.0){
		//					energyCloudFinal.add(point);
		//				}
		//			}
		//		}
		energyCloudFinal.sort(null);

		time = System.currentTimeMillis() - time;
		System.out.println("Time for calculating energy values: " + time/1000 + " sec");

		return energyCloudFinal;


	}

	private List<EnergyPoint> createEnergyCloudThreadedLoopNEW(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		/*
		 * PROGRESS BAR
		 */
		JFrame f = new JFrame("Processing...");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container content = f.getContentPane();
		JProgressBar progressBar = new JProgressBar(0, 100);
		progressBar.setValue(0);
		progressBar.setStringPainted(true);
		content.add(progressBar, BorderLayout.CENTER);
		f.setSize(400, 100);
		f.setVisible(true);

		long time = System.currentTimeMillis();

		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];

		//		List<Point3D> weList = we.getWeList();

		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loop
		 */
		int numberOfThreads = Runtime.getRuntime().availableProcessors();
		int counter = 0;
		long[] weSize = we.getSize();
		double weResolution = we.getResolution();
		double[] weMinValues = we.getMinValues();
		long weSizeTot = weSize[0] * weSize[1] * weSize[2];
		int portion = (int) weSizeTot/numberOfThreads;

		int progress = 0;
		int dProgress = 100/numberOfThreads;
		List<Thread> threads = new LinkedList<Thread>();
		List<OnlinePlanner> planners = new LinkedList<OnlinePlanner>();

		List<Point3D> currPoints = new ArrayList<Point3D>(portion);

		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		int[] index = new int[3];

		for (index[0] = 0; index[0] < weSize[0]; index[0]++)
			for (index[1] = 0; index[1] < weSize[1]; index[1]++)
				for (index[2] = 0; index[2] < weSize[2]; index[2]++){

					for(int j = 0; j < 3; j++)					
						buffer[j] = new BigDecimal(weMinValues[j] + index[j] * weResolution, mt);

					Point3D p3d = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());

					currPoints.add(counter, p3d);

					counter++;

					if(counter < portion)
						continue;

					counter = 0;

					OnlinePlanner planner = new OnlinePlanner(p, r, currPoints);

					Thread thread = new Thread( planner );

					planners.add(planner);
					threads.add(thread);

				}

		for(Thread currThread : threads)
			currThread.start();

		for(Thread currThread : threads){

			try {

				currThread.join();

				if(currThread.getState() == Thread.State.TERMINATED){

					System.out.println("Thread is terminated.");
					progress += dProgress;			
					progressBar.setValue(progress);
					progressBar.setStringPainted(true);
					progressBar.paint(progressBar.getGraphics());
					if(progress == 100)
						f.setVisible(false);

				}

			}

			catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

		for(OnlinePlanner elem : planners){

			for(EnergyPoint point : elem.getEnergyList())
				if(point.getEnergy() != 0.0)
					energyCloudFinal.add(point);

		}

		energyCloudFinal.sort(null);

		time = System.currentTimeMillis() - time;
		System.out.println("Time for calculating energy values: " + time/1000 + " sec");

		return energyCloudFinal;


	}

	private List<EnergyPoint> createEnergyCloudThreadedNEWNEW(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		/*
		 * PROGRESS BAR
		 */
		JFrame f = new JFrame("Processing...");
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container content = f.getContentPane();
		JProgressBar progressBar = new JProgressBar(0, 100);
		progressBar.setValue(0);
		progressBar.setStringPainted(true);
		content.add(progressBar, BorderLayout.CENTER);
		f.setSize(400, 100);
		f.setVisible(true);

		long time = System.currentTimeMillis();

		int[] index = new int[3];
		long[] weSize = we.getSize();
		long weSizeTot = weSize[0]*weSize[1]*weSize[2];
		System.out.println("weSize: " + weSize[0] + " " + weSize[1] + " " +  weSize[2]);
		System.out.println("weSizeTot: " + weSizeTot);
		long dProgress = (weSize[0]*weSize[1]*weSize[2])/100;
		System.out.println("dProgress = " + dProgress);
		int progress = 0;
		int counter = 0;
		int portion = (int) weSizeTot/100;
		int i;
		int iteration = 1;

		double weResolution = we.getResolution();
		double[] weMinValues = we.getMinValues();
		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];
		
		double[][] pos3H = new double[4][4];
		double[][] pos4H = new double[4][4];
		Target t3 = new Target(pos3H);
		Target t4 = new Target(pos4H);

		BigDecimal[] buffer = new BigDecimal[3];
		MathContext mt = new MathContext(2, RoundingMode.HALF_DOWN);

		//		List<Future<EnergyPoint>> energyCloud = new LinkedList<Future<EnergyPoint>>();
		List<Future<EnergyPoint>> energyCloud = new ArrayList<Future<EnergyPoint>>(portion);
		List<EnergyPoint> energyCloudFinal = new LinkedList<EnergyPoint>();

		Target[] targets = p.getPathPositions();
		
		System.out.println("targets.length = " + targets.length);

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());
		
		if(targets.length > 2){
			
			pos3H = Matrix.copyMatrix(targets[2].getHomMatrix());
			pos4H = Matrix.copyMatrix(targets[3].getHomMatrix());
			
		}


		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */
		ExecutorService execs = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

		for (index[0] = 0; index[0] < weSize[0]; index[0]++)
			for (index[1] = 0; index[1] < weSize[1]; index[1]++)
				for (index[2] = 0; index[2] < weSize[2]; index[2]++){

					for(int j = 0; j < 3; j++)					
						buffer[j] = new BigDecimal(weMinValues[j] + index[j] * weResolution, mt);
					
					Point3D point = new Point3D(buffer[0].doubleValue(), buffer[1].doubleValue(), buffer[2].doubleValue());

					Target[] curr = new Target[2];
					curr[0] = new Target(inPosH);
					curr[1] = new Target(finPosH);

					targets[0].setHomMatrix(inPosH);
					targets[1].setHomMatrix(finPosH);

					currPos = Matrix.copyVector(point.getPosition());

					curr[0].translateTarget(currPos);
					curr[1].translateTarget(currPos);
					
					if(targets.length > 2){
						
						targets[2].setHomMatrix(pos3H);
						targets[3].setHomMatrix(pos4H);

						t3.translateTarget(currPos);
						t4.translateTarget(currPos);
						
					}

//					OnlinePlannerCallable thread = new OnlinePlannerCallable(curr[0], curr[1], r, point);
					
					PlannerList thread = null;
					
					if(targets.length == 2){
					
					thread = new PlannerList(curr[0], curr[1], r, point);
					
					}else{
					
						Path sq = new SquarePath(curr[0], curr[1], t3, t4);
						sq.setMaxAcc(p.getMaxAcc());
						sq.setMaxVel(p.maxVel);
						thread = new PlannerList(sq, r, point);
						
					}
					
					Future<EnergyPoint> result = execs.submit(thread);

					energyCloud.add(counter, result);

					if(counter < portion - 1){

						counter++;
						continue;

					}
					
//					System.out.println("counter at the end of a list: " + counter);
//					System.out.println("Size List Future EP: " + energyCloud.size());

					execs.shutdown();
					
					counter = 0;

					i = 0;
					for(Future<EnergyPoint> ep : energyCloud){

						try {
							if(ep.get().getEnergy() != 0d){
								energyCloudFinal.add(energyCloud.get(i).get());
							}
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						} catch (ExecutionException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}

						i++;

					}

					progress++;			
					progressBar.setValue(progress);
					progressBar.setStringPainted(true);
					progressBar.paint(progressBar.getGraphics());
					if(progress == 100)
						f.setVisible(false);
					
					energyCloud = new ArrayList<Future<EnergyPoint>>(portion);
					
					if(iteration < weSizeTot)
					execs = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

				}

		time = System.currentTimeMillis() - time;
		
		System.out.println("Size of energyCloudFinal = " + energyCloudFinal.size());

		System.out.println("Time for calculating energy values: " + time/1000 + " sec");

		return energyCloudFinal;


	}
	
	private TransformGroup setPosition(double[] position) {

		/*
		 * just an easy method to avoid writing the following lines every time I
		 * need to set a certain position for an object
		 */

		TransformGroup tg = new TransformGroup();
		//		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		//		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);

		Transform3D t3d = new Transform3D();

		t3d.setTranslation(new Vector3d(position[0], position[1], position[2]));

		/*
		 * TO CHECK: does it improve the performances when the working envelope
		 * has > 20K points?
		 */
		t3d.setScale(0.05);

		tg.setTransform(t3d);

		return tg;

	}

	private Transform3D setView() {

		/*
		 * just an easy method to improve readability of the code It sets the
		 * center of the view
		 */
		double[] center = { 0.0, 0.45, 5.0 };

		Transform3D t3d = new Transform3D();
		Transform3D temp = new Transform3D();
		t3d.rotX(Math.PI/2);
		temp.setTranslation(new Vector3d(center[0], center[1], center[2]));
		t3d.mul(temp);

		return t3d;

	}

	private Light setDirectionalLight(Vector3f direction) {

		/*
		 * If you want the robot to be orange
		 */
		//		Color3f color = new Color3f(1.0f, 0.5f, 0.0f);

		/*
		 * If you want the robot to be white
		 */
		Color3f color = new Color3f(1.0f, 1.0f, 1.0f);

		/*
		 * It defines the sphere inside which the lights have influence
		 */
		BoundingSphere bounds = new BoundingSphere(new Point3d(0, 0, 0), 100);

		DirectionalLight light = new DirectionalLight(color, direction);

		light.setInfluencingBounds(bounds);

		return light;

	}

	private TransformGroup robotModel() {

		String path = "src/Robot 3D models/IRB 1600.obj";

		ObjectFile loader = new ObjectFile();
		Scene s = null;

		TransformGroup tg = new TransformGroup();
		Transform3D t3d = new Transform3D();

		//		t3d.rotZ(Math.PI/2);

		if(path.contains("1600")){

			t3d.setTranslation(new Vector3d(0.0, 0.0, 0.0));
			t3d.setScale(0.0011);

		}else{

			t3d.setTranslation(new Vector3d(-0.15, -0.15, 0.0));
			t3d.setScale(0.0007);

		}

		tg.setTransform(t3d);

		try {

			s = loader.load(path);

		} catch (Exception e) {

			System.err.println(e);
			System.exit(1);

		}

		tg.addChild(s.getSceneGroup());

		return tg;

	}

	private BranchGroup axis() {

		BranchGroup axisBG = new BranchGroup();

		// create line for X axis
		LineArray axisXLines = new LineArray(2, LineArray.COORDINATES | LineArray.COLOR_3);
		axisBG.addChild(new Shape3D(axisXLines));

		axisXLines.setCoordinate(0, new Point3f(-1.0f, 0.0f, 0.0f));
		axisXLines.setCoordinate(1, new Point3f(1.0f, 0.0f, 0.0f));

		// create line for Y axis
		LineArray axisYLines = new LineArray(2, LineArray.COORDINATES | LineArray.COLOR_3);
		axisBG.addChild(new Shape3D(axisYLines));

		axisYLines.setCoordinate(0, new Point3f(0.0f, -1.0f, 0.0f));
		axisYLines.setCoordinate(1, new Point3f(0.0f, 1.0f, 0.0f));

		// create line for Z axis
		Point3f z1 = new Point3f(0.0f, 0.0f, 0.0f);
		Point3f z2 = new Point3f(0.0f, 0.0f, 1.1f);

		LineArray axisZLines = new LineArray(10, LineArray.COORDINATES | LineArray.COLOR_3);
		axisBG.addChild(new Shape3D(axisZLines));

		axisZLines.setCoordinate(0, z1);
		axisZLines.setCoordinate(1, z2);
		axisZLines.setCoordinate(2, z2);
		axisZLines.setCoordinate(3, new Point3f(0.07f, 0.07f, 0.97f));
		axisZLines.setCoordinate(4, z2);
		axisZLines.setCoordinate(5, new Point3f(-0.07f, 0.07f, 0.97f));
		axisZLines.setCoordinate(6, z2);
		axisZLines.setCoordinate(7, new Point3f(0.07f, -0.07f, 0.97f));
		axisZLines.setCoordinate(8, z2);
		axisZLines.setCoordinate(9, new Point3f(-0.07f, -0.07f, 0.97f));

		Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
		Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		Color3f colors[] = new Color3f[9];

		for (int v = 0; v < 9; v++)
			colors[v] = white;

		axisXLines.setColor(0, white);
		axisXLines.setColor(1, red);
		axisYLines.setColor(0, white);
		axisYLines.setColor(1, white);

		axisZLines.setColors(1, colors);

		return axisBG;

	}

}




