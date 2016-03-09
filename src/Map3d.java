import java.applet.*;
import java.awt.*;

import javax.media.j3d.*;
import javax.vecmath.*;

import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.universe.SimpleUniverse;
import com.sun.j3d.utils.universe.PlatformGeometry;
import com.sun.j3d.utils.behaviors.keyboard.*;

import com.sun.j3d.loaders.Scene;

import java.awt.event.KeyListener;
import java.awt.event.KeyEvent;

import java.util.*;
import com.sun.j3d.utils.geometry.Sphere;

import org.jdesktop.j3d.loaders.vrml97.*;

import com.sun.j3d.utils.behaviors.vp.OrbitBehavior;

/*
 * TO DO:
 * - to add robot model;
 * - to create TransformGroup to set origin of the world (and hierarchical transformations);
 * - to add frame to show where the origin is;
 * - to add description of values of energy according to color of the sphere.
 */

public class Map3d {

	/*
	 * 
	 */

	private SimpleUniverse u = null;
	private Canvas3D canvas = null;
	// private TransformGroup viewtrans = null;

	public Map3d(String robotDirectory, Path p, WorkingEnvelope we) {

		double distanceView = 100;
		double[] centerView = { -0.1, 0.1, 2.0 };

		Robot r = new Robot(robotDirectory, true);

		u = new SimpleUniverse(canvas);

		BranchGroup scene = createSceneGraph(r, p, we);
		System.out.println("End creating SceneGraph");

		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView(centerView, distanceView));
		// u.getViewingPlatform().setNominalViewingTransform();
		// u.getViewer().getView().setBackClipDistance(100);

		setViewInteraction(canvas);

		u.addBranchGraph(scene);
		System.out.println("End of map3d");

	}

	public static void main(String[] args) {

		/*
		 * TEST
		 */

		String directory = "C:\\Users\\Leonardo Urbano\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML Costruttore Robot.xml";

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
		targets[0] = testTarget1;
		targets[1] = testTarget2;
		double x_sample = 1E-3;
		double t_sample = 0.01;
		double acceleration = 1.0;
		double velocity = 0.2;
		Path p = new Path(targets[0], targets[1], t_sample, x_sample);
		p.setMaxAcc(acceleration);
		p.setMaxVel(velocity);

		// double[] min = {-0.5, 0.0, 0.0};
		// double[] max = {1.0, 1.0, 1.0};
		// double res = 0.05;
		double[] min = { -0.5, 0.0, 0.0 };
		double[] max = { 0.5, 0.5, 0.5 };
		double res = 0.05;

		WorkingEnvelope we = new WorkingEnvelope(min, max, res);

		new Map3d (directory, p, we);

	}

	private BranchGroup createSceneGraph(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It adds all the objects to the world: 
		 * - energy point cloud; 
		 * - robot model;
		 */
		int[] indEnergyMatrix = new int[3];
		int[] size = we.getSize();
		int numFeasibleSol = size[0] * size[1] * size[2];
		int numPoint = 0;

		float[] energyColor = new float[3];
		// float[] CenterOfWorld = {0, 0, 0};

		double[] currPos = new double[3];
		double[][][] energyMatrix = createEnergyMatrix(r, p, we);
		double[][][][] weMatrix = we.getWorkingEnvelopeMatrix();

		BranchGroup objRoot = new BranchGroup();
		// TransformGroup cow = setPosition(CenterOfWorld);

		/*
		 * Here I call the method to add mouse-view interaction
		 */
		// setViewInteraction(canvas);

		for (int x = 0; x < weMatrix.length; x++) {
			for (int y = 0; y < weMatrix[0].length; y++) {
				for (int z = 0; z < weMatrix[0][0].length; z++) {

					numPoint++;

					indEnergyMatrix[0] = x;
					indEnergyMatrix[1] = y;
					indEnergyMatrix[2] = z;
					
					currPos = Matrix.copyVector(weMatrix[x][y][z]);

					if (energyMatrix[x][y][z] == 0.0) {

						numFeasibleSol--;
						continue;

					}

					System.out.println("(OUTSIDE METHOD, INSIDE LOOP) current energy value: "
							+ energyMatrix[indEnergyMatrix[0]][indEnergyMatrix[1]][indEnergyMatrix[2]]);

					energyColor = getEnergyColor(energyMatrix, indEnergyMatrix);

					objRoot.addChild(energyPoint(currPos, energyColor));

					System.out.println("End point number " + numPoint);

				}

			}

		}

		/*
		 * here I call the method robotModel to add the robot to the branch
		 */
		// objRoot.addChild(robotModel());

		System.out.println("Number of feasible solutions: " + numFeasibleSol + "/" + size[0] * size[1] * size[2]);

		return objRoot;

	}

	private BranchGroup energyPoint(double[] position, float[] rgb) {

		System.out.println("(INSIDE ENERGYPOINT) rgb: " + rgb[0] + " " + rgb[1] + " " + rgb[2]);

		BranchGroup objRoot = new BranchGroup();

		TransformGroup tg = setPosition(position);

		Sphere s = new Sphere(0.15f);

		/*
		 * setting the appearance of the sphere (the color)
		 */
		Appearance app = new Appearance();
		Color3f color = new Color3f(rgb[0], rgb[1], rgb[2]);
		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.FASTEST);
		app.setColoringAttributes(ca);
		s.setAppearance(app);

		tg.addChild(s);

		objRoot.addChild(tg);

		/*
		 * Here I call the method setLight to add lights to the branch
		 */
		objRoot.addChild(setLight());

		objRoot.compile();

		return objRoot;

	}

	private float[] getEnergyColor(double[][][] energyMatrix, int[] pos) {

		/*
		 * Method that maps energy values into colors. It first gets min and max
		 * values of energyMatrix to create the array energyIntervals. It then
		 * scans energyIntervals to get the right index to use in order to
		 * assign the right rgb value to the current energy value
		 */

		int numIntervals = 6;
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
		double dRgb = 510 / numIntervals;
		double[] energyIntervals = new double[numIntervals + 1];

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

		n = energyIntervals.length - 1 - index;

		System.out.println(
				"(INSIDE METHOD, AFTER INTERVALS) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
		System.out.println("interval number: " + n);

		/*
		 * It starts from rgb(255, 0, 0) = red and it does the following: 
		 * - if n <= 3 it increases the value of "g" until you get (255, 255, 0) = yellow;
		 * - if n > 3 it start decreasing the value of "r" until you get (0,255,0) = green.
		 * Doing so you get red -> orange -> yellow -> green
		 */
		if (n > 3) {

			rgb[0] -= (n - 3) * dRgb;
			rgb[1] += 3 * dRgb;

		} else

			rgb[1] += n * dRgb;

		System.out.println("The rgb is: " + rgb[0] + ", " + rgb[1] + ", " + rgb[2]);

		return rgb;

	}

	private double[][][] createEnergyMatrix(Robot r, Path p, WorkingEnvelope we) {

		/*
		 * It creates a 3D matrix with energy values for every single starting
		 * position inside the working envelope
		 */

		int[] size = we.getSize();

		double[] currPos = new double[3];
		double[][] inPosH = new double[4][4];
		double[][] finPosH = new double[4][4];
		double[][][] energyMatrix = new double[size[0]][size[1]][size[2]];
		double[][][][] weMatrix = we.getWorkingEnvelopeMatrix();

		Target[] targets = p.getPathPositions();

		inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
		finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());

		/*
		 * Here I initialize the OnlinePlanner to then create a thread each
		 * iteration of the following loops
		 */
		OnlinePlanner dynSolution = new OnlinePlanner(targets, r);

		for (int x = 0; x < weMatrix.length; x++) {
			for (int y = 0; y < weMatrix[0].length; y++) {
				for (int z = 0; z < weMatrix[0][0].length; z++) {

					Target[] curr = new Target[2];
					curr[0] = new Target(inPosH);
					curr[1] = new Target(finPosH);

					targets[0].setHomMatrix(inPosH);
					targets[1].setHomMatrix(finPosH);

					currPos = Matrix.copyVector(weMatrix[x][y][z]);

					System.out.println("currPos = ");
					Matrix.displayVector(currPos);

					curr[0].translateTarget(currPos);
					curr[1].translateTarget(currPos);

					System.out.println("Current initial position: ");
					Matrix.displayMatrix(curr[0].getHomMatrix());

					System.out.println("Current final position: ");
					Matrix.displayMatrix(curr[1].getHomMatrix());

					/*
					 * Here I create the thread to calculate the energy value
					 * with current initial and final position
					 */
					dynSolution.run(curr[0], curr[1]);

					/*
					 * Every time the Online Planner calculates an energy value
					 * it stores it in energyList. In order to get the current
					 * energy value I have to "ask" for the last element added
					 * to the list
					 */
					energyMatrix[x][y][z] = dynSolution.getEnergyList().get(dynSolution.getEnergyList().size() - 1);

					System.out.println(
							"energyMatrix[" + x + "]" + "[" + y + "]" + "[" + z + "]" + " = " + energyMatrix[x][y][z]);

				}
			}
		}

		return energyMatrix;

	}

	private TransformGroup setPosition(double[] position) {

		/*
		 * just an easy method to avoid writing the following lines every time I
		 * need to set a certain position for an object
		 */

		TransformGroup tg = new TransformGroup();

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

	private void setViewInteraction(Canvas3D c) {

		/*
		 * TO DO: to add mouse interactions instead of moving the view with the
		 * keyboard
		 */

		/*
		 * It adds mouse-view interaction
		 */

		// BoundingSphere bounds = new BoundingSphere(new Point3d(), 10000);
		//
		// OrbitBehavior orbit = new OrbitBehavior(c,
		// OrbitBehavior.REVERSE_ROTATE);
		// orbit.setSchedulingBounds(bounds);
		//
		// u.getViewingPlatform().setViewPlatformBehavior(orbit);

		BoundingSphere bounds = new BoundingSphere(new Point3d(), 10000);

		// viewtrans = u.getViewingPlatform().getViewPlatformTransform();

		KeyNavigatorBehavior key = new KeyNavigatorBehavior(u.getViewingPlatform().getViewPlatformTransform());
		key.setSchedulingBounds(bounds);

		PlatformGeometry platformGeom = new PlatformGeometry();
		platformGeom.addChild(key);

		u.getViewingPlatform().setPlatformGeometry(platformGeom);

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

		Color3f color = new Color3f(1.8f, 0.1f, 0.1f);

		/*
		 * It defines the sphere inside which the lights have influence
		 */
		BoundingSphere bounds = new BoundingSphere(new Point3d(0, 0, 0), 100);

		Vector3f direction = new Vector3f(4.0f, -7.0f, -12.0f);

		DirectionalLight light = new DirectionalLight(color, direction);

		light.setInfluencingBounds(bounds);

		return light;

	}

	private BranchGroup robotModel() {

		/*
		 * TO DO: to add the other parts!!!
		 */

		/*
		 * It creates the robot objects by adding each robot part with the right
		 * positioning (one part wrt the adiacent one)
		 */

		BranchGroup objRoot = new BranchGroup();

		double[] pos = { 0, 0, 0 };

		TransformGroup tg = setPosition(pos);

		VrmlLoader loader = new VrmlLoader();

		Scene s = null;

		try {

			s = loader.load("irb1600/Fundament.wrl");

		} catch (Exception e) {

			System.err.println(e);
			System.exit(1);

		}

		tg.addChild(s.getSceneGroup());

		objRoot.addChild(tg);

		objRoot.compile();

		return objRoot;

	}

}
