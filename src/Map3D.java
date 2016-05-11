
import java.applet.*;
import java.awt.*;

import javax.media.j3d.*;
import javax.vecmath.*;

import javax.swing.*;

import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.universe.SimpleUniverse;
import com.sun.j3d.utils.universe.PlatformGeometry;
import com.sun.j3d.utils.behaviors.keyboard.*;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.sun.glass.events.MouseEvent;
import com.sun.j3d.loaders.Scene;

import java.awt.event.KeyListener;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.awt.event.KeyEvent;
import com.sun.j3d.utils.behaviors.mouse.MouseRotate;
import com.sun.j3d.utils.behaviors.mouse.MouseZoom;
import com.sun.j3d.utils.behaviors.mouse.MouseTranslate;

import java.util.*;

import com.sun.j3d.utils.geometry.Primitive;
import com.sun.j3d.utils.geometry.Sphere;
import com.sun.j3d.utils.picking.PickResult;

import org.jdesktop.j3d.loaders.vrml97.*;
import com.sun.j3d.loaders.objectfile.ObjectFile;

import com.sun.j3d.utils.behaviors.vp.OrbitBehavior;


//import wise.wiseshopfloor.core.J3DViewer;

//import wise.wiseshopfloor.modelbase.*;
//import wise.wiseshopfloor.universe.CameraState;


/*
 * TO DO:
 * - to add robot model;
 * - to create TransformGroup to set origin of the world (and hierarchical transformations);
 * - to add frame to show where the origin is;
 * - to add description of values of energy according to color of the sphere.
 */

public class Map3D extends JApplet{
	
	private static final Comparator Comparator = null;
	/*
	 * 
	 */
	private SimpleUniverse u = null;
	private Canvas3D canvas = null;

	public Map3D(String xmlFilePath, Path p, WorkingEnvelope we) {

		double distanceView = 110;
		double[] centerView = { -0.1, 0.45, 3.5 };

		Robot r = new Robot(xmlFilePath, true);

		/*
		 * here I create a simple canvas
		 */
		setLayout(new BorderLayout());
		GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
		canvas = new Canvas3D(config);
		add("North", new Label(""));
		add("Center", canvas);
		add("South", new Label(""));

		u = new SimpleUniverse(canvas);

		BranchGroup scene = createSceneGraph(r, p, we);
		System.out.println("End creating SceneGraph");

		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView(centerView, distanceView));

		u.addBranchGraph(scene);
		System.out.println("End of Map3D");

	}
	
	public Map3D(String xmlFilePath, List<EnergyPoint> energyList){
		double distanceView = 110;
		double[] centerView = { -0.1, 0.45, 3.5 };

		Robot r = new Robot(xmlFilePath, true);
		
		setLayout(new BorderLayout());
		GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
		canvas = new Canvas3D(config);
		add("North", new Label(""));
		add("Center", canvas);
		add("South", new Label(""));

		u = new SimpleUniverse(canvas);

		BranchGroup scene = createSceneGraph(energyList);
		System.out.println("End creating SceneGraph");

		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView(centerView, distanceView));

		u.addBranchGraph(scene);
		System.out.println("End of Map3D");
	}

	public static void main(String[] args) {
		
//		JProgressBar progressBar = new JProgressBar(0, 100);
//		progressBar.setValue(0);
//        progressBar.setStringPainted(true);
//		JApplet applet = new JApplet();
//		applet.add(progressBar);

		/*
		 * TEST
		 */

//		String xmlFilePath = "C:\\Users\\Leonardo Urbano\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML ABB IRB 140.xml";
//		String xmlFilePath = "C:\\Users\\Enrico\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML ABB IRB 140.xml";
		String xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";

		double[][] h1 = { { 0.0, 0.0, 1.0, 0.59627 }, 
				{ 1.0, 0.0, 0.0, 0.0 }, 
				{ 0.0, 1.0, 0.0, 0.66282 },
				{ 0.0, 0.0, 0.0, 1.0 } };
		double[][] h2 = { { 0.0, 0.0, 1.0, 0.59627 },
				{ 1.0, 0.0, 0.0, -0.21132 }, 
				{ 0.0, 1.0, 0.0, 0.58041 },
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

//		double[] min = {-0.5, 0.0, 0.0};
//		double[] max = {1.0, 1.0, 1.0};
//		double res = 0.05;
				double[] min = { -1.0, -1.0, 0.0 };
				double[] max = { 1.0, 1.0, 1.0 };
				double res = 0.05;

		long time = System.currentTimeMillis();
		
		WorkingEnvelope we = new WorkingEnvelope(min, max, res);

		//		new Map3d (xmlFilePath, p, we);
		
		Map3D demo = new Map3D(xmlFilePath, p, we);
		new MainFrame(demo, 400, 400);

//		Matrix.displayVector(backToCow);
//		Matrix.displayMatrix(testTarget1.getHomMatrix());
		
		time = System.currentTimeMillis() - time;
		
		System.out.println("Overall time: " + time/1000 + " sec");


	}

	private BranchGroup createSceneGraph(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It adds all the objects to the world: 
		 * - energy point cloud; 
		 * - robot model;
		 * - reference frame.
		 */
		
//		List<EnergyPoint> energyCloud = createEnergyCloud(r, p, we);
		List<EnergyPoint> energyCloud = createEnergyCloudThreaded(r, p, we);
//		List<EnergyPoint> energyCloud = createEnergyCloudThreadedLoop(r, p, we);
		
		int numFeasibleSol = energyCloud.size();
		int numPoint = 0;
		int[] size = we.getSize();

		float[] energyColor = new float[3];
		//		float[] energyColor = {0, 255, 0};

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

		/*
		 * Here I set Mouse interaction (rotation, translation and zoom)
		 */
//		BoundingSphere bounds = new BoundingSphere(new Point3d(), 10000);
//
//		MouseRotate mr = new MouseRotate();
//		mr.setTransformGroup(cow);
//		cow.addChild(mr);
//		mr.setSchedulingBounds(bounds);
//
//		MouseTranslate mt = new MouseTranslate();
//		mt.setTransformGroup(cow);
//		cow.addChild(mt);
//		mt.setSchedulingBounds(bounds);
//
//		MouseZoom mz = new MouseZoom();
//		mz.setTransformGroup(cow);
//		cow.addChild(mz);
//		mz.setSchedulingBounds(bounds);
		
	    OrbitBehavior ob = new OrbitBehavior(u.getCanvas());
	    ob.setReverseRotate(true);
	    ob.setReverseTranslate(true);
	    ob.setSchedulingBounds(new BoundingSphere(new Point3d(0.0,0.0,0.0),Double.MAX_VALUE));
	    u.getViewingPlatform().setViewPlatformBehavior(ob);
		
		int numClusters = 11;
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

		objRoot.addChild(setLight());

		System.out.println("Number of feasible solutions: " + numFeasibleSol + "/" + size[0] * size[1] * size[2]);

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
		
		int numClusters = 11;
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

		objRoot.addChild(setLight());

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
        	rgb[1] -= (r.cluster[i + 1] - Math.ceil(numClusters/2))*0.2f;
        }else
        	rgb[0] += r.cluster[i + 1]*0.2f;
        
        	        	
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

		Sphere s = new Sphere(0.14f);

		/*
		 * Setting the appearance of the sphere (color and transparency)
		 */

		Appearance app = new Appearance();
		Color3f color = ep.getColor();
		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.FASTEST);
		app.setColoringAttributes(ca);

		TransparencyAttributes ta = new TransparencyAttributes();
		ta.setTransparencyMode(TransparencyAttributes.BLENDED);
		ta.setTransparency(0.1f);
		app.setTransparencyAttributes(ta);

		s.setAppearance(app);

		tg.addChild(s);

		return tg;

	}

	private float[] getEnergyColor(double[][][] energyMatrix, int[] pos, double[] valuePropagation) {

		/*
		 * Method that maps energy values into colors. It first gets min and max
		 * values of energyMatrix to create the array energyIntervals. It then
		 * scans energyIntervals to get the right index to use in order to
		 * assign the right rgb value to the current energy value
		 */

		int numIntervals = 4;
		int index;
		int n;
		int[] indMax = Matrix.getIndMaxValue3D(energyMatrix);
		int[] indMin = Matrix.getIndMinValue3D(energyMatrix);

		float[] rgb = { 255, 0, 0 };

		System.out.println("current index = " + pos[0] + " " + pos[1] + " " + pos[2]);

		/*
		 * Here I calculate energy and rgb resolution according to the number of
		 * intervals I decided to have
		 */
		double dEnergy = (energyMatrix[indMax[0]][indMax[1]][indMax[2]] - energyMatrix[indMin[0]][indMin[1]][indMin[2]])
				/ numIntervals;
		//		
		//		double dRgb = 255;
		//		double[] energyIntervals = new double[numIntervals + 1];
		double[] energyIntervals = new double[numIntervals];

		for (int i = 0; i < energyIntervals.length; i++)

			energyIntervals[i] = energyMatrix[indMin[0]][indMin[1]][indMin[2]] + i * dEnergy;

		System.out.println("energyIntervals: ");
		Matrix.displayVector(energyIntervals);

		/*
		 * TO DO: see if it's possible to implement something like binary search to speed up 
		 * the search of the index
		 */
		index = energyIntervals.length - 1;
		while (index >= 0) {

			if (energyMatrix[pos[0]][pos[1]][pos[2]] >= energyIntervals[index])

				break;

			index--;

		}

		valuePropagation[index]++;

		n = energyIntervals.length - 1 - index;

		//		System.out.println(
		//				"(INSIDE METHOD, AFTER INTERVALS) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
		System.out.println("interval number: " + n);

		/*
		 * It starts from rgb(255, 0, 0) = red and it does the following: 
		 * - if n <= 3 it increases the value of "g" until you get (255, 255, 0) = yellow;
		 * - if n > 3 it start decreasing the value of "r" until you get (0,255,0) = green.
		 * Doing so you get red -> orange -> yellow -> green
		 */
		//		if (n > 3) {
		//
		//			rgb[0] -= (n - 3) * dRgb;
		//			rgb[1] += 3 * dRgb;
		//
		//		} else
		//
		//			rgb[1] += n * dRgb;

		if(index == 0){
			rgb[0] -= 255;
			rgb[1] += 255;
		}
		if(index == 1)
			rgb[1] += 125;
		if(index == 2)
			rgb[1] += 255;

		System.out.println("The rgb is: " + rgb[0] + ", " + rgb[1] + ", " + rgb[2]);

		return rgb;

	}

//	private float[] getEnergyColorKmeans(double[][][] energyMatrix, int[] pos, double[] valuePropagation) {
//
//		/*
//		 * Method that maps energy values into colors. It first gets min and max
//		 * values of energyMatrix to create the array energyIntervals. It then
//		 * scans energyIntervals to get the right index to use in order to
//		 * assign the right rgb value to the current energy value
//		 */
//
//		int numClusters = 10;
//
//		double[] centroids = new double[numClusters];
//
//		/*
//		 * Random initialization of centroids with values inside energyMatrix
//		 */
//		for(int i = 0; i < centroids.length; i++){
//
//			centroids[i] = energyMatrix[(int)(Math.random()*energyMatrix.length)][(int)(Math.random()*energyMatrix[0].length)][(int)(Math.random()*energyMatrix[0][0].length)];
//
//			/*
//			 * Need to ensure all centroids get initialized in different points!!
//			 */
//
//		}
//
//		for(int i = 0; i < energyMatrix.length; i++)
//			for(int j = 0; j < energyMatrix[0].length; j++)
//				for(int k = 0; k < energyMatrix[0][0].length; k++){
//
//
//
//				}


		//		System.out.println(
		//				"(INSIDE METHOD, AFTER INTERVALS) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
//		System.out.println("interval number: " + n);

		/*
		 * It starts from rgb(255, 0, 0) = red and it does the following: 
		 * - if n <= 3 it increases the value of "g" until you get (255, 255, 0) = yellow;
		 * - if n > 3 it start decreasing the value of "r" until you get (0,255,0) = green.
		 * Doing so you get red -> orange -> yellow -> green
		 */
		//		if (n > 3) {
		//
		//			rgb[0] -= (n - 3) * dRgb;
		//			rgb[1] += 3 * dRgb;
		//
		//		} else
		//
		//			rgb[1] += n * dRgb;
		//		
		//		if(index == 0){
		//			rgb[0] -= 255;
		//			rgb[1] += 255;
		//		}
		//		if(index == 1)
		//			rgb[1] += 125;
		//		if(index == 2)
		//			rgb[1] += 255;
		//
		//		System.out.println("The rgb is: " + rgb[0] + ", " + rgb[1] + ", " + rgb[2]);
		//
		//		return rgb;
//
//	}

	private List<EnergyPoint> createEnergyCloud(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];
		
		List<Point3D> weList = we.getWeList();
		List<EnergyPoint> energyCloud = new LinkedList<EnergyPoint>();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */
		OnlinePlanner dynSolution = new OnlinePlanner(targets, r);

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
					dynSolution.run(curr[0], curr[1], point);

					/*
					 * Every time the Online Planner calculates an energy value
					 * it stores it in energyList. In order to get the current
					 * energy value I have to "ask" for the last element added
					 * to the list
					 */
					
					double energy = dynSolution.getEnergyList().get(dynSolution.getEnergyList().size() - 1).getEnergy();
					
					if(energy != 0.0){
						
					EnergyPoint ep = new EnergyPoint();
					
					ep.setPosition(point.getPosition());
					ep.setEnergy(energy);
					
					energyCloud.add(ep);
					
					}

//					System.out.println(
//							"energyMatrix[" + x + "]" + "[" + y + "]" + "[" + z + "]" + " = " + energyMatrix[x][y][z]);

		}

		return energyCloud;

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
					OnlinePlannerCallable thread = new OnlinePlannerCallable(curr[0], curr[1], r, point);
					
					Future<EnergyPoint> result = execs.submit(thread);
					
					energyCloud.add(result);
				
	
		}
		
		execs.shutdown();
		
		int i = 0;
		for(Point3D point : weList){
			try {
				if(energyCloud.get(i).get().getEnergy() != 0.0){
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
//		ExecutorService execs = Executors.newFixedThreadPool(numberOfThreads);
		
		for(int i = 0; i < numberOfThreads; i++){
			int from = i*portion;
			int to = i < numberOfThreads-1 ? (i+1)*portion : weList.size();
			
			OnlinePlanner planner = new OnlinePlanner(p, r, weList.subList(from, to-1));
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
		
		time = System.currentTimeMillis() - time;
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

	private Transform3D setView(double[] center, double distance) {

		/*
		 * just an easy method to improve readability of the code It sets the
		 * center of the view and then the distanceof the viewer from it
		 */

		Transform3D t3d = new Transform3D();
		t3d.setTranslation(new Vector3d(center[0], center[1], center[2]));
		u.getViewer().getView().setBackClipDistance(distance);

		return t3d;

	}

	private Light setLight() {

		Color3f color = new Color3f(1.5f, 1.5f, 1.5f);

		/*
		 * It defines the sphere inside which the lights have influence
		 */
		BoundingSphere bounds = new BoundingSphere(new Point3d(0, 0, 0), 100);

		Vector3f direction = new Vector3f(4.0f, -7.0f, -12.0f);

		DirectionalLight light = new DirectionalLight(color, direction);
		//		AmbientLight light = new AmbientLight(color);

		light.setInfluencingBounds(bounds);

		return light;

	}

	private BranchGroup robotModel() {

		BranchGroup objRoot = new BranchGroup();

		TransformGroup tg = new TransformGroup();

		Transform3D t3d = new Transform3D();
		Transform3D temp = new Transform3D();

		t3d.setTranslation(new Vector3d(0.0, 0.0, 0.0));
		
		t3d.rotZ(Math.PI/2);
		temp.rotZ(-Math.PI/2);

		t3d.mul(temp);

		t3d.setScale(0.0007);

		tg.setTransform(t3d);

		ObjectFile loader = new ObjectFile();

		Scene s = null;

		try {

			s = loader.load("src/Robot 3D models/IRB 1600.obj");

		} catch (Exception e) {

			System.err.println(e);
			System.exit(1);

		}
		
//		Appearance app = new Appearance();
//		TransparencyAttributes ta = new TransparencyAttributes();
//		ta.setTransparencyMode(TransparencyAttributes.BLENDED);
//		ta.setTransparency(0.3f);
//		app.setTransparencyAttributes(ta);

		tg.addChild(s.getSceneGroup());

		objRoot.addChild(tg);

		objRoot.compile();

		return objRoot;



		//		BranchGroup objRoot = new BranchGroup();
		//
		//		TransformGroup tg = new TransformGroup();
		//
		//		Transform3D t3d = new Transform3D();
		//
		//		t3d.setTranslation(new Vector3d(0.0, 0.0, 0.0));
		//		
		//		VrmlLoader loader = new VrmlLoader();
		//
		//		Scene s = null;
		//
		//		try {
		//
		//			s = loader.load("irb1600/IRB1600.3ds");
		//
		//		} catch (Exception e) {
		//
		//			System.err.println(e);
		//			System.exit(1);
		//
		//		}
		//
		//		tg.addChild(s.getSceneGroup());
		//
		//		objRoot.addChild(tg);
		//
		//		objRoot.compile();
		//
		//		return objRoot;

	}

	private BranchGroup axis() {

		BranchGroup axisBG = new BranchGroup();

		// create line for X axis
		LineArray axisXLines = new LineArray(2, LineArray.COORDINATES);
		axisBG.addChild(new Shape3D(axisXLines));

		axisXLines.setCoordinate(0, new Point3f(-1.0f, 0.0f, 0.0f));
		axisXLines.setCoordinate(1, new Point3f(1.0f, 0.0f, 0.0f));

		//        Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		//        Color3f green = new Color3f(0.0f, 1.0f, 0.0f);
		//        Color3f blue = new Color3f(0.0f, 0.0f, 1.0f);
		Color3f white = new Color3f(1.0f, 1.0f, 1.0f);

		// create line for Y axis
		LineArray axisYLines = new LineArray(2, LineArray.COORDINATES
				| LineArray.COLOR_3);
		axisBG.addChild(new Shape3D(axisYLines));

		axisYLines.setCoordinate(0, new Point3f(0.0f, -1.0f, 0.0f));
		axisYLines.setCoordinate(1, new Point3f(0.0f, 1.0f, 0.0f));

		axisYLines.setColor(0, white);
		axisYLines.setColor(1, white);

		// create line for Z axis
		Point3f z1 = new Point3f(0.0f, 0.0f, -1.0f);
		Point3f z2 = new Point3f(0.0f, 0.0f, 1.0f);

		LineArray axisZLines = new LineArray(10, LineArray.COORDINATES
				| LineArray.COLOR_3);
		axisBG.addChild(new Shape3D(axisZLines));

		axisZLines.setCoordinate(0, z1);
		axisZLines.setCoordinate(1, z2);
		axisZLines.setCoordinate(2, z2);
		axisZLines.setCoordinate(3, new Point3f(0.1f, 0.1f, 0.9f));
		axisZLines.setCoordinate(4, z2);
		axisZLines.setCoordinate(5, new Point3f(-0.1f, 0.1f, 0.9f));
		axisZLines.setCoordinate(6, z2);
		axisZLines.setCoordinate(7, new Point3f(0.1f, -0.1f, 0.9f));
		axisZLines.setCoordinate(8, z2);
		axisZLines.setCoordinate(9, new Point3f(-0.1f, -0.1f, 0.9f));

		Color3f colors[] = new Color3f[9];

		colors[0] = new Color3f(0.0f, 1.0f, 1.0f);
		for (int v = 0; v < 9; v++) {
			colors[v] = white;
		}

		axisZLines.setColors(1, colors);

		return axisBG;

	}
	
}



