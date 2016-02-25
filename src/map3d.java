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
		u.getViewingPlatform().getViewPlatformTransform().setTransform(setView());
		//u.getViewingPlatform().setNominalViewingTransform();
		u.getViewer().getView().setBackClipDistance(100);
		setViewInteraction(canvas);
		u.addBranchGraph(scene);

	}

	private BranchGroup createSceneGraph(Robot r, Path p, Working_Envelope we){
		
		int i = 0;
		int j = 0;
		int k = 0;
		int[] indEnergyMatrix = new int[3];
		//		float[] CenterOfWorld = {0, 0, 0};
		float[] energyColor = new float[3];
		double[] currPos = new double[3];
		double[][][] energyMatrix = createEnergyMatrix(r, p, we);
		BranchGroup objRoot = new BranchGroup();
		//		TransformGroup cow = setPosition(CenterOfWorld);

		//		setViewInteraction(canvas);

		for(double x = we.getMinValues()[0]; x <= we.getMaxValues()[0]; x += we.getResolution()){
			for(double y = we.getMinValues()[1]; y <= we.getMaxValues()[1]; y += we.getResolution())
			{
				for(double z = we.getMinValues()[2]; z <= we.getMaxValues()[2]; z += we.getResolution())
				{
					
					indEnergyMatrix[0] = i;
					indEnergyMatrix[1] = j;
					indEnergyMatrix[2] = k;
					
					currPos[0] = x;
					currPos[1] = y;
					currPos[2] = z;
					
					energyColor = getEnergyColor(energyMatrix, indEnergyMatrix);
					//objRoot.addChild(energyPoint(currPos, energyColor));
					objRoot.addChild(energyPoint(currPos, energyColor));
					
					k++;
					
				}
				
				j++;
				
			}
			
			i++;
		}
			

//		objRoot.addChild(robotModel());

		return objRoot;

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
		t3d.setTranslation(new Vector3d(0, 0, 0));

		return t3d;

	}	

	private BranchGroup energyPoint(double[] position, float[] rgb){

		BranchGroup objRoot = new BranchGroup();
		TransformGroup tg = setPosition(position);
		
		Sphere s = new Sphere(0.01f);
		
		Appearance app = new Appearance();
		Color3f color = new Color3f(rgb[0], rgb[1], rgb[2]);
		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.NICEST);
		app.setColoringAttributes(ca);
		s.setAppearance(app);

		tg.addChild(s);

		objRoot.addChild(tg);

		objRoot.addChild(setLight());

		objRoot.compile();

		return objRoot;

	}
	
//	private BranchGroup energyPoint(double[] position, float[] rgb){
//
//		BranchGroup objRoot = new BranchGroup();
//		TransformGroup tg = setPosition(position);
//
//		Sphere s = new Sphere(0.01f);
//		
////		float[] rgb = energyColor(position);
//
//		Appearance app = new Appearance();
//		Color3f color = new Color3f(rgb[0], rgb[1], rgb[2]);
//		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.NICEST);
//		app.setColoringAttributes(ca);
//		s.setAppearance(app);
//
//		tg.addChild(s);
//
//		objRoot.addChild(tg);
//
//		objRoot.addChild(setLight());
//
//		objRoot.compile();
//
//		return objRoot;
//
//	}


	//	private BranchGroup energyPoint(TransformGroup tg1, float[] position, int[] rgb){
	//
	//		BranchGroup objRoot = new BranchGroup();
	//		TransformGroup tg2 = setPosition(position);
	//
	//		Sphere s = new Sphere(0.01f);
	//
	//		Appearance app = new Appearance();
	//		Color3f color = new Color3f(rgb[0], rgb[1], rgb[2]);
	//		ColoringAttributes ca = new ColoringAttributes(color, ColoringAttributes.NICEST);
	//		app.setColoringAttributes(ca);
	//		s.setAppearance(app);
	//
	//		tg2.addChild(s);
	//		
	//		tg1.addChild(tg2);
	//		
	//		objRoot.addChild(tg1);
	//		
	//		objRoot.addChild(setLight());
	//		
	//		objRoot.compile();
	//		
	//		return objRoot;
	//
	//	}
	
	private float[] getEnergyColor(double[][][] energyMatrix, int[] pos){
		
		float[] rgb = {0, 255, 0};
		int[] indMax = Matrix.getIndMaxValue3D(energyMatrix);
		int[] indMin = Matrix.getIndMinValue3D(energyMatrix);
		int n = 10;
		int k;
		double dEnergy = (energyMatrix[indMax[0]][indMax[1]][indMax[2]] - energyMatrix[indMin[0]][indMin[1]][indMin[2]])/n;
		double dRgb = 510/n;
		
		double[] energyIntervals = new double[n];
		for(int i = 0; i < energyIntervals.length; i++)
			energyIntervals[i] = energyMatrix[indMin[0]][indMin[1]][indMin[2]] + dEnergy*energyMatrix[indMin[0]][indMin[1]][indMin[2]];
		
		for(k = energyIntervals.length - 1; k >= 0; k--)
			if(energyIntervals[k - 1] >= energyMatrix[pos[0]][pos[1]][pos[2]])
				break;
		
		if(k > 5){
			rgb[0] += 5*dRgb;
			rgb[1] -= (k - 5)*dRgb;
		}else
			rgb[0] += k*dRgb;
		
		return rgb;
		
	}
	
	private double[][][] createEnergyMatrix(Robot r, Path p, Working_Envelope we){
		
		int[] size = we.getSize();
		int i = 0;
		int j = 0;
		int k = 0;
		double[] currPos = new double[3];

		double[][][] energyMatrix = new double[size[0]][size[1]][size[2]];
		
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
		
		/*
		 * First of all, I need to get the 3d matrix with energy values
		 */
		for(double x = we.getMinValues()[0]; x <= we.getMaxValues()[0]; x += we.getResolution())
		{
			for(double y = we.getMinValues()[1]; y <= we.getMaxValues()[1]; y += we.getResolution())
			{
				for(double z = we.getMinValues()[2]; z <= we.getMaxValues()[2]; z += we.getResolution())
				{
					
					currPos[0] = x;
					currPos[1] = y;
					currPos[2] = z;
					
					targets[0].translateTarget(currPos);
					targets[1].translateTarget(currPos);
					
					OnlinePlanner dynSolution = new OnlinePlanner(targets, r);
					
					dynSolution.run();
					
					energyMatrix[i][j][k] = dynSolution.getEnergyList().get(0) ;
					
					k++;
					
				}
				
				j++;
				
			}
			
		i++;
		
	}
	
		return energyMatrix;
		
	}

	private TransformGroup setPosition(double[] position){

		TransformGroup tg = new TransformGroup();

		Transform3D t3d = new Transform3D();

		t3d.setTranslation(new Vector3d(position[0], position[1], position[2]));

		tg.setTransform(t3d);

		return tg;

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
		
		double[] minValues = {0, 0, 0};
		double[] maxValues = {20, 20, 20};
		double resolution = 1;
		/*
		 * I'll have to find a way to be able to set Resolution < 1
		 * I'll have to change float to double and make things work
		 */
		Working_Envelope we = new Working_Envelope(minValues, maxValues, resolution);

		new map3d(directory, p, we);

	}

}
