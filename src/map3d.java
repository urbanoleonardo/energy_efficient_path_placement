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
 * - funzione che definisce colori in base a valore energia
 * - aggiungere tutto il robot, non solo un pezzo
 * - posizionare robot dove si vuole
 * - creare trasformate gerarchiche, in particolare settare 0 sullo 0 della frame zero
 * - legenda per descrivere relazione colori-energia
 */

public class map3d{

	private SimpleUniverse u = null;
	private Canvas3D canvas = null;
	//	private TransformGroup viewtrans = null;

	public map3d(String robotDirectory, Path p, Working_Envelope we){


		Robot r = new Robot(robotDirectory, true);

		u = new SimpleUniverse(canvas);

		BranchGroup scene = createSceneGraph(r, p, we);
		System.out.println("End creating Scene Graph");
		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView());
		//u.getViewingPlatform().setNominalViewingTransform();
		u.getViewer().getView().setBackClipDistance(100);
		setViewInteraction(canvas);
		u.addBranchGraph(scene);
		System.out.println("End of everything");

	}

	private BranchGroup createSceneGraph(Robot r, Path p, Working_Envelope we){

		int numPoint = 0;
		int[] size = we.getSize();
		int numFeasibleSol = size[0]*size[1]*size[2];

//		int i = 0;
//		int j = 0;
//		int k = 0;
		int[] indEnergyMatrix = new int[3];
		//		float[] CenterOfWorld = {0, 0, 0};
		float[] energyColor = new float[3];
		double[] currPos = new double[3];
		double[][][] energyMatrix = createEnergyMatrix(r, p, we);
		double[][][][] weMatrix = we.getWorkingEnvelopeMatrix();
		BranchGroup objRoot = new BranchGroup();
		//		TransformGroup cow = setPosition(CenterOfWorld);

		//		setViewInteraction(canvas);
		
//		for(int x = 0; x < weMatrix.length; x++){
//			for(int y = 0; y < weMatrix[0].length; y++){
//				for(int z = 0; z < weMatrix[0][0].length; z++){
		for(int x = 15; x < weMatrix.length; x++){
			for(int y = 3; y < weMatrix[0].length; y++){
				for(int z = 2; z < weMatrix[0][0].length; z++){
//					int x = 1;
//					int y = 3; 
//					int z = 1;
					
					numPoint++;
					
					indEnergyMatrix[0] = x;
					indEnergyMatrix[1] = y;
					indEnergyMatrix[2] = z;
		
					currPos = Matrix.copyVector(weMatrix[x][y][z]);
					
//					System.out.println("currPos: " + currPos[0] + " " + currPos[1] + " " + currPos[2]);
					
					if(energyMatrix[x][y][z] == 0.0){
						numFeasibleSol--;
						continue;
					}

					System.out.println("(OUTSIDE METHOD, INSIDE LOOP) current indEnergyMatrix: " + indEnergyMatrix[0] + " " + indEnergyMatrix[1] + " " + indEnergyMatrix[2]);
					
					System.out.println("(OUTSIDE METHOD, INSIDE LOOP) current energy value: " + energyMatrix[indEnergyMatrix[0]][indEnergyMatrix[1]][indEnergyMatrix[2]]);
					System.out.println("(OUTSIDE METHOD, INSIDE LOOP) current energy value in [1][3][1]: " + energyMatrix[1][3][1]);

					
					energyColor = getEnergyColor(energyMatrix, indEnergyMatrix);
					//objRoot.addChild(energyPoint(currPos, energyColor));
					
					System.out.println("Start drawing point number " + numPoint);
					
					objRoot.addChild(energyPoint(currPos, energyColor));
					
					System.out.println("End point number " + numPoint);

				}
			}
		}


		//		objRoot.addChild(robotModel());

		System.out.println("Number of feasible solutions: " + numFeasibleSol + "/" + size[0]*size[1]*size[2]);
		return objRoot;

	}

	private BranchGroup energyPoint(double[] position, float[] rgb){
	
		System.out.println("(INSIDE ENERGYPOINT) rgb: " + rgb[0] + " " + rgb[1] + " " + rgb[2]);
		
		BranchGroup objRoot = new BranchGroup();
		TransformGroup tg = setPosition(position);
	
//		Sphere s = new Sphere(0.01f);
				Sphere s = new Sphere(0.15f);
	
		Appearance app = new Appearance();
		Color3f color = new Color3f(rgb[0], rgb[1], rgb[2]);
		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.FASTEST);
		app.setColoringAttributes(ca);
		s.setAppearance(app);
	
		tg.addChild(s);
	
		objRoot.addChild(tg);
	
		objRoot.addChild(setLight());
	
		objRoot.compile();
	
		return objRoot;
	
	}

	private float[] getEnergyColor(double[][][] energyMatrix, int[] pos){
		
			System.out.println("(INSIDE METHOD) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
			System.out.println("(INSIDE METHOD) current energy value in [1][3][1]: " + energyMatrix[1][3][1]);
		
			float[] rgb = {255, 0, 0};
			int[] indMax = Matrix.getIndMaxValue3D(energyMatrix);
			int[] indMin = Matrix.getIndMinValue3D(energyMatrix);
	
			System.out.println("(INSIDE METHOD, AFTER MIN, MAX) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
			System.out.println("(INSIDE METHOD, AFTER MIN, MAX) current energy value in [1][3][1]: " + energyMatrix[1][3][1]);
			
			System.out.println("indMax = " + indMax[0] + " " + indMax[1]  + " " + indMax[2]);
			System.out.println("max = " + energyMatrix[indMax[0]][indMax[1]][indMax[2]]);
			System.out.println("indMin = " + indMin[0] + " " + indMin[1]  + " " + indMin[2]);
			System.out.println("min = " + energyMatrix[indMin[0]][indMin[1]][indMin[2]]);
			System.out.println("current index = " + pos[0] + " " + pos[1] + " " + pos[2]);
	
			int n = 10;
			int k;
			double dEnergy = (energyMatrix[indMax[0]][indMax[1]][indMax[2]] - energyMatrix[indMin[0]][indMin[1]][indMin[2]])/n;
			double dRgb = 510/n;
			
			double[] energyIntervals = new double[n + 1];
			for(int i = 0; i < energyIntervals.length; i++)
				energyIntervals[i] = energyMatrix[indMin[0]][indMin[1]][indMin[2]] + i*dEnergy;
			
			System.out.println("energyIntervals: ");
			Matrix.displayVector(energyIntervals);
	
					for(k = energyIntervals.length - 1; k >= 0; k--)
						if(energyMatrix[pos[0]][pos[1]][pos[2]] >= energyIntervals[k])
								break;
			System.out.println("(INSIDE METHOD, AFTER INTERVALS) current energy value: " + energyMatrix[pos[0]][pos[1]][pos[2]]);
			System.out.println("interval number: " + (energyIntervals.length-1-k));
	
					if((energyIntervals.length - 1 - k) > 5){
						rgb[0] -= 5*dRgb;
						rgb[1] += (energyIntervals.length - 1 - k - 5)*dRgb;
					}else
						rgb[0] -= (energyIntervals.length - 1 - k)*dRgb;
	
			System.out.println("Then the rgb is: " + rgb[0] + ", " + rgb[1] + ", " + rgb[2]);
	
			return rgb;
	
		}

	private double[][][] createEnergyMatrix(Robot r, Path p, Working_Envelope we){
	
			int[] size = we.getSize();
			//		int i = 0;
			//		int j = 0;
			//		int k = 0;
			//		double x = we.getMinValues()[0];
			//		double y = we.getMinValues()[1];
			//		double z = we.getMinValues()[2];
			//		double[] currPos = {0, 0, 0};
			double[] currPos = new double[3];
			double[][] inPosH = new double[4][4];
			double[][] finPosH = new double[4][4];
			double[][][] energyMatrix = new double[size[0]][size[1]][size[2]];
			double[][][][] weMatrix = we.getWorkingEnvelopeMatrix();
	
			Target[] targets = p.getPathPositions();
							
							inPosH = Matrix.copyMatrix(targets[0].getHomMatrix());
							finPosH = Matrix.copyMatrix(targets[1].getHomMatrix());
	
	//				System.out.println("Initial target: ");
	//				Matrix.displayMatrix(inPosH);
	//				System.out.println("Final target: ");
	//				Matrix.displayMatrix(finPosH);
	
			OnlinePlanner dynSolution = new OnlinePlanner(targets, r);
	
	//		for(int x = 6; x < 12; x++){
			for(int x = 0; x < weMatrix.length; x++){
				for(int y = 0; y < weMatrix[0].length; y++){
					for(int z = 0; z < weMatrix[0][0].length; z++){
	//					int x = 8;
	//					int y = 0;
	//					int z = 0;
				
						Target[] curr = new Target[2];
						curr[0] = new Target(inPosH);
						curr[1] = new Target(finPosH);
						
						targets[0].setHomMatrix(inPosH);
						targets[1].setHomMatrix(finPosH);
						
						currPos = Matrix.copyVector(weMatrix[x][y][z]);
						
						System.out.println("currPos = ");
						Matrix.displayVector(currPos);
						
	//					System.out.println("(BEFORE) current initial position: ");
	//					Matrix.displayMatrix(curr[0].getHomMatrix());
	//					
	//					System.out.println("(BEFORE) current final position: ");
	//					Matrix.displayMatrix(curr[1].getHomMatrix());
						
						curr[0].translateTarget(currPos);
						curr[1].translateTarget(currPos);
						
						System.out.println("Current initial position: ");
						Matrix.displayMatrix(curr[0].getHomMatrix());
						
						System.out.println("Current final position: ");
						Matrix.displayMatrix(curr[1].getHomMatrix());
						
//						OnlinePlanner dynSolution = new OnlinePlanner(curr, r);
	
						dynSolution.run(curr[0], curr[1]);
	
						energyMatrix[x][y][z] = dynSolution.getEnergyList().get(dynSolution.getEnergyList().size() - 1);
	
	
						System.out.println("energyMatrix[" + x + "]" + "[" + y + "]" + "[" + z + "]" + " = " + energyMatrix[x][y][z]);
	
					}
				}
			}	
	
			System.out.println("Size of energyMatrix: [" + size[0] + " " + size[1] + " " + size[2] + "]");
			
			return energyMatrix;
	
		}

	private TransformGroup setPosition(double[] position){
	
		TransformGroup tg = new TransformGroup();
	
		Transform3D t3d = new Transform3D();
	
		t3d.setTranslation(new Vector3d(position[0], position[1], position[2]));
		
		t3d.setScale(0.05);
	
		tg.setTransform(t3d);
	
		return tg;
	
	}

	private void setViewInteraction(Canvas3D c){

		//		BoundingSphere bounds = new BoundingSphere(new Point3d(), 10000);
		//
		//		OrbitBehavior orbit = new OrbitBehavior(c, OrbitBehavior.REVERSE_ROTATE);
		//		orbit.setSchedulingBounds(bounds);
		//
		//		u.getViewingPlatform().setViewPlatformBehavior(orbit);

		BoundingSphere bounds = new BoundingSphere(new Point3d(), 10000);
		//		viewtrans = u.getViewingPlatform().getViewPlatformTransform();
		KeyNavigatorBehavior key = new KeyNavigatorBehavior(u.getViewingPlatform().getViewPlatformTransform());
		key.setSchedulingBounds(bounds);
		PlatformGeometry platformGeom = new PlatformGeometry();
		platformGeom.addChild(key);
		u.getViewingPlatform().setPlatformGeometry(platformGeom);

	}

	private Transform3D setView(){

		Transform3D t3d = new Transform3D();
		t3d.setTranslation(new Vector3d(-0.1, 0.1, 2));

		return t3d;

	}	

	private Light setLight(){

		Color3f color = new Color3f(1.8f, 0.1f, 0.1f);

		BoundingSphere bounds = new BoundingSphere(new Point3d(0, 0, 0), 100);

		Vector3f direction = new Vector3f(4.0f, -7.0f, -12.0f);

		DirectionalLight light = new DirectionalLight(color, direction);

		light.setInfluencingBounds(bounds);

		return light;

	}

	private BranchGroup robotModel(){

		BranchGroup objRoot = new BranchGroup();

		VrmlLoader loader = new VrmlLoader();
		Scene s = null;

		try{
			s = loader.load("irb1600/Fundament.wrl");
		} catch (Exception e){
			System.err.println(e);
			System.exit(1);
		}

		double[] pos = {0, 0, 0};
		TransformGroup tg = setPosition(pos);

		tg.addChild(s.getSceneGroup());

		objRoot.addChild(tg);

		objRoot.compile();

		return objRoot;

	}

	//	private BranchGroup robotModel(TransformGroup tg){
	//		
	//		BranchGroup objRoot = new BranchGroup();
	//		
	//		VrmlLoader loader = new VrmlLoader();
	//		Scene s = null;
	//		
	//		try{
	//			s = loader.load("irb1600/IRB1600_Base.wrl");
	//		} catch (Exception e){
	//			System.err.println(e);
	//			System.exit(1);
	//		}
	//		
	//		//float[] pos = {0, 0, 0};
	//		//TransformGroup tg = setPosition(pos);
	//		
	//		tg.addChild(s.getSceneGroup());
	//		
	//		objRoot.addChild(tg);
	//		
	//		objRoot.compile();
	//
	//		return objRoot;
	//		
	//	}

	public static void main(String[] args) {

		String directory = "C:\\Users\\Leonardo Urbano\\Google Drive\\Tesi\\Prova Costruttore Robot\\XML Costruttore Robot.xml";

		double[][] h1 = {
				{0.0, 0.0, 1.0, 0.59627},
				{1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.66282},
				{0.0, 0.0, 0.0, 1.0}
		};
		double[][] h2 = {
				{0.0, 0.0, 1.0, 0.59627},
				{1.0, 0.0, 0.0, -0.21132},
				{0.0, 1.0, 0.0, 0.58041},
				{0.0, 0.0, 0.0, 1.0}};

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

//		double[] minValues = {-0.5, 0.0, 0.0};
//		double[] maxValues = {1.0, 1.0, 1.0};
//		double resolution = 0.05;
				double[] minValues = {-0.5, 0.0, 0.0};
				double[] maxValues = {0.3, 0.3, 0.3};
				double resolution = 0.05;
		
		Working_Envelope we = new Working_Envelope(minValues, maxValues, resolution);
		
		
		
//		double[][][][] m = we.getWorkingEnvelopeMatrix();
//		Robot r = new Robot(directory, true);
//		Target[] t = p.getPathPositions();
//		double[] pos = m[0][0][0];
//		System.out.println("pos = ");
//		Matrix.displayVector(pos);
//		t[0].translateTarget(pos);
//		t[1].translateTarget(pos);
//		OnlinePlanner dynSolution = new OnlinePlanner(t, r);
//		dynSolution.run(t[0], t[1]);
		
//
		new map3d(directory, p, we);

	}

}
