import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.custom.CLabel;
import org.eclipse.swt.SWT;
//import org.eclipse.wb.swt.SWTResourceManager;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Text;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.widgets.ProgressBar;
import org.eclipse.swt.widgets.List;
import org.eclipse.swt.widgets.MessageBox;
import org.eclipse.swt.widgets.CoolBar;
import org.eclipse.swt.widgets.Combo;
import org.eclipse.swt.events.ModifyListener;
import org.eclipse.swt.events.ModifyEvent;

public class UserInterface {

	protected Shell shell;
	private Text textResolution;
	private Text textInPos00;
	private Text textInPos10;
	private Text textInPos20;
	private Text textInPos30;
	private Text textInPos01;
	private Text textInPos11;
	private Text textInPos21;
	private Text textInPos31;
	private Text textInPos02;
	private Text textInPos12;
	private Text textInPos22;
	private Text textInPos32;
	private Text textInPos03;
	private Text textInPos13;
	private Text textInPos23;
	private Text textInPos33;
	private Text textFinPos00;
	private Text textFinPos01;
	private Text textFinPos02;
	private Text textFinPos03;
	private Text textFinPos13;
	private Text textFinPos12;
	private Text textFinPos11;
	private Text textFinPos10;
	private Text textFinPos20;
	private Text textFinPos30;
	private Text textFinPos31;
	private Text textFinPos21;
	private Text textFinPos22;
	private Text textFinPos32;
	private Text textFinPos33;
	private Text textFinPos23;
	private Text textMax0;
	private Text textMax1;
	private Text textMax2;
	private Text textMin0;
	private Text textMin1;
	private Text textMin2;
	private Text textPos230;
	private Text textPos231;
	private Text textPos232;
	private Text textPos233;
	private Text textPos223;
	private Text textPos222;
	private Text textPos221;
	private Text textPos220;
	private Text textPos210;
	private Text textPos211;
	private Text textPos212;
	private Text textPos213;
	private Text textPos203;
	private Text textPos202;
	private Text textPos201;
	private Text textPos200;
	private Label lblndPosition;
	private Label lblrdPosition;
	private Text textPos300;
	private Text textPos310;
	private Text textPos320;
	private Text textPos330;
	private Text textPos331;
	private Text textPos321;
	private Text textPos311;
	private Text textPos301;
	private Text textPos302;
	private Text textPos303;
	private Text textPos313;
	private Text textPos312;
	private Text textPos322;
	private Text textPos323;
	private Text textPos333;
	private Text textPos332;
	private Text textVelocity;
	private Text textAcceleration;
	private Label lblNewLabel_2;
	private Label lblNewLabel_3;
	
	private boolean irb140;
	private boolean irb1600;
	
	private String trajectory;
	private String velocity;
	private String acceleration;
	
	private String inPos00;
	private String inPos10;
	private String inPos20;
	private String inPos01;
	private String inPos11;
	private String inPos21;
	private String inPos02;
	private String inPos12;
	private String inPos22;
	private String inPos03;
	private String inPos13;
	private String inPos23;
	private String inPos30;
	private String inPos31;
	private String inPos32;
	private String inPos33;
	
	private String pos200;
	private String pos210;
	private String pos220;
	private String pos201;
	private String pos211;
	private String pos221;
	private String pos202;
	private String pos212;
	private String pos222;
	private String pos203;
	private String pos213;
	private String pos223;
	private String pos230;
	private String pos231;
	private String pos232;
	private String pos233;
	
	private String pos300;
	private String pos310;
	private String pos320;
	private String pos301;
	private String pos311;
	private String pos321;
	private String pos302;
	private String pos312;
	private String pos322;
	private String pos303;
	private String pos313;
	private String pos323;
	private String pos330;
	private String pos331;
	private String pos332;
	private String pos333;

	private String finPos00;
	private String finPos01;
	private String finPos02;
	private String finPos03;
	private String finPos13;
	private String finPos12;
	private String finPos11;
	private String finPos10;
	private String finPos20;
	private String finPos21;
	private String finPos22;
	private String finPos23;
	private String finPos30;
	private String finPos31;
	private String finPos32;
	private String finPos33;
	
	private String max0;
	private String max1;
	private String max2;
	
	private String min0;
	private String min1;
	private String min2;
	
	private String resolution;
	

	/**
	 * Launch the application.
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			UserInterface window = new UserInterface();
			window.open();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Open the window.
	 */
	public void open() {
		Display display = Display.getDefault();
		createContents();
		shell.open();
		shell.layout();
		while (!shell.isDisposed()) {
			if (!display.readAndDispatch()) {
				display.sleep();
			}
		}
	}

	/**
	 * Create contents of the window.
	 */
	protected void createContents() {
		shell = new Shell();
		shell.setSize(1258, 819);
		shell.setText("Map3D");
		
		Label lblRobot = new Label(shell, SWT.NONE);
		lblRobot.setBounds(204, 10, 139, 25);
		lblRobot.setText("ROBOT");
		
		Label lblWorkingEnvelope = new Label(shell, SWT.NONE);
		lblWorkingEnvelope.setBounds(152, 159, 175, 25);
		lblWorkingEnvelope.setText("WORKING ENVELOPE");
		
		Label lblMinvalues = new Label(shell, SWT.NONE);
		lblMinvalues.setBounds(62, 227, 98, 25);
		lblMinvalues.setText("Min Values");
		
		Label lblMaxvalues = new Label(shell, SWT.NONE);
		lblMaxvalues.setBounds(62, 357, 98, 25);
		lblMaxvalues.setText("Max Values");
		
		Label lblResolution = new Label(shell, SWT.NONE);
		lblResolution.setBounds(62, 478, 98, 25);
		lblResolution.setText("Resolution");
		
		Button btnIrb140 = new Button(shell, SWT.RADIO);
		btnIrb140.setSelection(true);
		btnIrb140.setBounds(105, 90, 133, 25);
		btnIrb140.setText("IRB 140");
		
		Button btnIrb1600 = new Button(shell, SWT.RADIO);
		btnIrb1600.setBounds(244, 90, 133, 25);
		btnIrb1600.setText("IRB 1600");
		
		textResolution = new Text(shell, SWT.BORDER);
		textResolution.setBounds(243, 475, 39, 31);
		
		textInPos00 = new Text(shell, SWT.BORDER);
		textInPos00.setBounds(606, 185, 39, 31);
		
		textInPos10 = new Text(shell, SWT.BORDER);
		textInPos10.setBounds(606, 222, 39, 31);
		
		textInPos20 = new Text(shell, SWT.BORDER);
		textInPos20.setBounds(606, 259, 39, 31);
		
		Label lblInitialPosition = new Label(shell, SWT.NONE);
		lblInitialPosition.setBounds(464, 225, 116, 25);
		lblInitialPosition.setText("Initial Position");
		
		Label lbldefaultValues = new Label(shell, SWT.NONE);
		lbldefaultValues.setBounds(158, 264, 243, 25);
		lbldefaultValues.setText("( default: -1.0, -1.0, 0.0 [m] )");
		
		Label lbldefaultValues_1 = new Label(shell, SWT.NONE);
		lbldefaultValues_1.setText("( default: 1.0, 1.0, 1.0 [m] )");
		lbldefaultValues_1.setBounds(160, 394, 230, 25);
		
		Label lbldefault = new Label(shell, SWT.NONE);
		lbldefault.setText("( default: 0.05 [m] )");
		lbldefault.setBounds(174, 512, 173, 25);
		
		Label lblFinalPosition = new Label(shell, SWT.NONE);
		lblFinalPosition.setText("Final Position");
		lblFinalPosition.setBounds(464, 457, 122, 25);
		
		Combo comboTrajectory = new Combo(shell, SWT.READ_ONLY);
		comboTrajectory.addSelectionListener(new SelectionAdapter() {
			@Override
			public void widgetSelected(SelectionEvent e) {
				
				if(comboTrajectory.getText().contains("line")){
					
					textPos200.setEnabled(false);
					textPos201.setEnabled(false);
					textPos202.setEnabled(false);
					textPos203.setEnabled(false);
					textPos210.setEnabled(false);
					textPos211.setEnabled(false);
					textPos212.setEnabled(false);
					textPos213.setEnabled(false);
					textPos220.setEnabled(false);
					textPos221.setEnabled(false);
					textPos222.setEnabled(false);
					textPos223.setEnabled(false);
					textPos230.setEnabled(false);
					textPos231.setEnabled(false);
					textPos232.setEnabled(false);
					textPos233.setEnabled(false);
					
					textPos300.setEnabled(false);
					textPos301.setEnabled(false);
					textPos302.setEnabled(false);
					textPos303.setEnabled(false);
					textPos310.setEnabled(false);
					textPos311.setEnabled(false);
					textPos312.setEnabled(false);
					textPos313.setEnabled(false);
					textPos320.setEnabled(false);
					textPos321.setEnabled(false);
					textPos322.setEnabled(false);
					textPos323.setEnabled(false);
					textPos330.setEnabled(false);
					textPos331.setEnabled(false);
					textPos332.setEnabled(false);
					textPos333.setEnabled(false);
					
					lblndPosition.setEnabled(false);
					lblrdPosition.setEnabled(false);
					
				}else{
						
						textPos200.setEnabled(true);
						textPos201.setEnabled(true);
						textPos202.setEnabled(true);
						textPos203.setEnabled(true);
						textPos210.setEnabled(true);
						textPos211.setEnabled(true);
						textPos212.setEnabled(true);
						textPos213.setEnabled(true);
						textPos220.setEnabled(true);
						textPos221.setEnabled(true);
						textPos222.setEnabled(true);
						textPos223.setEnabled(true);
						textPos230.setEnabled(true);
						textPos231.setEnabled(true);
						textPos232.setEnabled(true);
						textPos233.setEnabled(true);
						
						textPos300.setEnabled(true);
						textPos301.setEnabled(true);
						textPos302.setEnabled(true);
						textPos303.setEnabled(true);
						textPos310.setEnabled(true);
						textPos311.setEnabled(true);
						textPos312.setEnabled(true);
						textPos313.setEnabled(true);
						textPos320.setEnabled(true);
						textPos321.setEnabled(true);
						textPos322.setEnabled(true);
						textPos323.setEnabled(true);
						textPos323.setEnabled(true);
						textPos330.setEnabled(true);
						textPos331.setEnabled(true);
						textPos332.setEnabled(true);
						textPos333.setEnabled(true);
						
						lblndPosition.setEnabled(true);
						lblrdPosition.setEnabled(true);
					
				}
				
			}
		});
		comboTrajectory.setItems(new String[] {"Straight line", "Squared path"});
		comboTrajectory.setBounds(596, 87, 143, 33);
		comboTrajectory.select(0);
		
		Label lblTrajectory = new Label(shell, SWT.NONE);
		lblTrajectory.setBounds(499, 90, 81, 25);
		lblTrajectory.setText("Trajectory");
		
		Label lblPath = new Label(shell, SWT.NONE);
		lblPath.setBounds(785, 10, 48, 25);
		lblPath.setText("PATH");
		
		Label lblwarningFinalPosition = new Label(shell, SWT.NONE);
		lblwarningFinalPosition.setText("( warning: final position must be initialized )");
		lblwarningFinalPosition.setBounds(651, 578, 342, 25);
		
		Label lblwarningInitialPosition = new Label(shell, SWT.NONE);
		lblwarningInitialPosition.setText("( warning: initial position must be initialized )");
		lblwarningInitialPosition.setBounds(651, 343, 357, 25);
		
		Button btnStart = new Button(shell, SWT.NONE);
		btnStart.addSelectionListener(new SelectionAdapter() {
			@Override
			public void widgetSelected(SelectionEvent e) {
				
				irb140 = btnIrb140.getSelection();
				irb1600 = btnIrb1600.getSelection();
				
				velocity = textVelocity.getText();
				acceleration = textAcceleration.getText();
				
				inPos00 = textInPos00.getText();
				inPos01 = textInPos01.getText();
				inPos02 = textInPos02.getText();
				inPos03 = textInPos03.getText();
				inPos10 = textInPos10.getText();
				inPos11 = textInPos11.getText();
				inPos12 = textInPos12.getText();
				inPos13 = textInPos13.getText();
				inPos20 = textInPos20.getText();
				inPos21 = textInPos21.getText();
				inPos22 = textInPos22.getText();
				inPos23 = textInPos23.getText();
				inPos30 = textInPos30.getText();
				inPos31 = textInPos31.getText();
				inPos32 = textInPos32.getText();
				inPos33 = textInPos33.getText();
				
				pos200 = textPos200.getText();
				pos201 = textPos201.getText();
				pos202 = textPos202.getText();
				pos203 = textPos203.getText();
				pos210 = textPos210.getText();
				pos211 = textPos211.getText();
				pos212 = textPos212.getText();
				pos213 = textPos213.getText();
				pos220 = textPos220.getText();
				pos221 = textPos221.getText();
				pos222 = textPos222.getText();
				pos223 = textPos223.getText();
				pos230 = textPos230.getText();
				pos231 = textPos231.getText();
				pos232 = textPos232.getText();
				pos233 = textPos233.getText();
				
				pos300 = textPos300.getText();
				pos301 = textPos301.getText();
				pos302 = textPos302.getText();
				pos303 = textPos303.getText();
				pos310 = textPos310.getText();
				pos311 = textPos311.getText();
				pos312 = textPos312.getText();
				pos313 = textPos313.getText();
				pos320 = textPos320.getText();
				pos321 = textPos321.getText();
				pos322 = textPos322.getText();
				pos323 = textPos323.getText();
				pos330 = textPos330.getText();
				pos331 = textPos331.getText();
				pos332 = textPos332.getText();
				pos333 = textPos333.getText();
				
				finPos00 = textFinPos00.getText();
				finPos01 = textFinPos01.getText();
				finPos02 = textFinPos02.getText();
				finPos03 = textFinPos03.getText();
				finPos10 = textFinPos10.getText();
				finPos11 = textFinPos11.getText();
				finPos12 = textFinPos12.getText();
				finPos13 = textFinPos13.getText();
				finPos20 = textFinPos20.getText();
				finPos21 = textFinPos21.getText();
				finPos22 = textFinPos22.getText();
				finPos23 = textFinPos23.getText();
				finPos30 = textInPos30.getText();
				finPos31 = textInPos31.getText();
				finPos32 = textInPos32.getText();
				finPos33 = textInPos33.getText();
				
				min0 = textMin0.getText();
				min1 = textMin1.getText();
				min2 = textMin2.getText();
				
				max0 = textMax0.getText();
				max1 = textMax1.getText();
				max2 = textMax2.getText();
				
				resolution = textResolution.getText();
				
				trajectory = comboTrajectory.getText();
				
				if 	(inPos00 == null || inPos00.isEmpty() || inPos01 == null || inPos01.isEmpty()
					|| inPos02 == null || inPos02.isEmpty() || inPos03 == null || inPos03.isEmpty()
					|| inPos10 == null || inPos10.isEmpty() || inPos11 == null || inPos11.isEmpty()
					|| inPos12 == null || inPos12.isEmpty() || inPos13 == null || inPos13.isEmpty()
					|| inPos20 == null || inPos20.isEmpty() || inPos21 == null || inPos21.isEmpty()
					|| inPos22 == null || inPos22.isEmpty() || inPos23 == null || inPos23.isEmpty()
					|| finPos00 == null || finPos00.isEmpty() || finPos01 == null || finPos01.isEmpty()
					|| finPos02 == null || finPos02.isEmpty() || finPos03 == null || finPos03.isEmpty()
					|| finPos10 == null || finPos10.isEmpty() || finPos11 == null || finPos11.isEmpty()
					|| finPos12 == null || finPos12.isEmpty() || finPos13 == null || finPos13.isEmpty()
					|| finPos20 == null || finPos20.isEmpty() || finPos21 == null || finPos21.isEmpty()
					|| finPos22 == null || finPos22.isEmpty() || finPos23 == null || finPos23.isEmpty())
				{
					
					String errorMsg = null;
					MessageBox messageBox = new MessageBox(shell, SWT.OK | SWT.ICON_ERROR);
					
					messageBox.setText("Warning");
					if (inPos00 == null || inPos00.isEmpty() || inPos01 == null || inPos01.isEmpty()
						|| inPos02 == null || inPos02.isEmpty() || inPos03 == null || inPos03.isEmpty()
						|| inPos10 == null || inPos10.isEmpty() || inPos11 == null || inPos11.isEmpty()
						|| inPos12 == null || inPos12.isEmpty() || inPos13 == null || inPos13.isEmpty()
						|| inPos20 == null || inPos20.isEmpty() || inPos21 == null || inPos21.isEmpty()
						|| inPos22 == null || inPos22.isEmpty() || inPos23 == null || inPos23.isEmpty()) {
						errorMsg = "Please enter initial position";
					} else if (finPos00 == null || finPos00.isEmpty() || finPos01 == null || finPos01.isEmpty()
							|| finPos02 == null || finPos02.isEmpty() || finPos03 == null || finPos03.isEmpty()
							|| finPos10 == null || finPos10.isEmpty() || finPos11 == null || finPos11.isEmpty()
							|| finPos12 == null || finPos12.isEmpty() || finPos13 == null || finPos13.isEmpty()
							|| finPos20 == null || finPos20.isEmpty() || finPos21 == null || finPos21.isEmpty()
							|| finPos22 == null || finPos22.isEmpty() || finPos23 == null || finPos23.isEmpty()) {
						errorMsg = "Please enter final position";
					}
					if (errorMsg != null) {
						messageBox.setMessage(errorMsg);
						messageBox.open();
					}
				} else {
					
					shell.close();
					
					/*
					 * Here I have to initialize Map3D
					 */
					
					// ROBOT
					
					String xmlFilePath;
					
					if (irb140){
						
						xmlFilePath = "src/Robot Constructors/XML ABB IRB 140.xml";
						
					}else{
						
						xmlFilePath = "src/Robot Constructors/XML ABB IRB 1600 TEST.xml";
						
					}
					
					// WORKING ENVELOPE
					
					double res;
					double[] min = new double[3];
					double[] max = new double[3];
					WorkingEnvelope we;
					
					if (resolution == null || resolution.isEmpty())
						
						res = 0.05;
					
					else{
						
						res = Double.parseDouble(resolution);
						
					}
					
					if (min0 == null || min0.isEmpty())
						
						min[0] = -1.0;
					
					else{
						
						min[0] = Double.parseDouble(min0);
						
					}
					
					if (min1 == null || min1.isEmpty())
						
						min[1] = -1.0;
					
					else{
						
						min[1] = Double.parseDouble(min1);
						
					}
					
					if (min2 == null || min2.isEmpty())
						
						min[2] = 0.0;
					
					else{
						
						min[2] = Double.parseDouble(min2);
						
					}
					
					if (max0 == null || max0.isEmpty())
						
						max[0] = 1.0;
					
					else{
						
						max[0] = Double.parseDouble(max0);
						
					}
					
					if (max1 == null || max1.isEmpty())
						
						max[1] = 1.0;
					
					else{
						
						max[1] = Double.parseDouble(max1);
						
					}
					
					if (max2 == null || max2.isEmpty())
						
						max[2] = 1.0;
					
					else{
						
						max[2] = Double.parseDouble(max2);
						
					}
					
					we = new WorkingEnvelope(min, max, res);
					
					// PATH
					
					double vel;
					double acc;
					double xSample = 1E-3;
					double tSample = 0.01;
					double[][] inPos = {
						{Double.parseDouble(inPos00), Double.parseDouble(inPos01), Double.parseDouble(inPos02), Double.parseDouble(inPos03)},
						{Double.parseDouble(inPos10), Double.parseDouble(inPos11), Double.parseDouble(inPos12), Double.parseDouble(inPos13)},
						{Double.parseDouble(inPos20), Double.parseDouble(inPos21), Double.parseDouble(inPos22), Double.parseDouble(inPos23)},
						{Double.parseDouble(inPos30), Double.parseDouble(inPos31), Double.parseDouble(inPos32), Double.parseDouble(inPos33)}
						};
					double[][] finPos = {
						{Double.parseDouble(finPos00), Double.parseDouble(finPos01), Double.parseDouble(finPos02), Double.parseDouble(finPos03)},
						{Double.parseDouble(finPos10), Double.parseDouble(finPos11), Double.parseDouble(finPos12), Double.parseDouble(finPos13)},
						{Double.parseDouble(finPos20), Double.parseDouble(finPos21), Double.parseDouble(finPos22), Double.parseDouble(finPos23)},
						{Double.parseDouble(finPos30), Double.parseDouble(finPos31), Double.parseDouble(finPos32), Double.parseDouble(finPos33)}
						};
					Target in = new Target(inPos);
					Target fin = new Target(finPos);
					Path p;
					
					if (trajectory.contains("Straight line"))
						p = new Path(in, fin, tSample, xSample);
					else{
						
						double[][] pos2 = {
								{Double.parseDouble(pos200), Double.parseDouble(pos201), Double.parseDouble(pos202), Double.parseDouble(pos203)},
								{Double.parseDouble(pos210), Double.parseDouble(pos211), Double.parseDouble(pos212), Double.parseDouble(pos213)},
								{Double.parseDouble(pos220), Double.parseDouble(pos221), Double.parseDouble(pos222), Double.parseDouble(pos223)},
								{Double.parseDouble(pos230), Double.parseDouble(pos231), Double.parseDouble(pos232), Double.parseDouble(pos233)}
								};
						double[][] pos3 = {
								{Double.parseDouble(pos300), Double.parseDouble(pos301), Double.parseDouble(pos302), Double.parseDouble(pos303)},
								{Double.parseDouble(pos310), Double.parseDouble(pos311), Double.parseDouble(pos312), Double.parseDouble(pos313)},
								{Double.parseDouble(pos320), Double.parseDouble(pos321), Double.parseDouble(pos322), Double.parseDouble(pos323)},
								{Double.parseDouble(pos330), Double.parseDouble(pos331), Double.parseDouble(pos332), Double.parseDouble(pos333)}
								};
						Target p2 = new Target(pos2);
						Target p3 = new Target(pos3);
						
						p = new SquarePath(in, p2, p3, fin, tSample, xSample);
						
					}
					
					if (velocity == null || velocity.isEmpty())
						
						vel = 0.2;
					
					else{
						
						vel = Double.parseDouble(velocity);
						
					}
					
					if (acceleration == null || acceleration.isEmpty())
						
						acc = 1;
					
					else{
						
						acc = Double.parseDouble(acceleration);
						
					}
					
					p.setMaxAcc(acc);
					p.setMaxVel(vel);
					
					new Map3D(xmlFilePath, p, we);
					
//					MessageBox messageBox = new MessageBox(shell, SWT.OK | SWT.ICON_WORKING);
//					messageBox.setText("Info");
//					messageBox.setMessage("Valid");
//					messageBox.open();
					
				}
			}
		});
		
		btnStart.setBounds(540, 666, 105, 35);
		btnStart.setText("Start");
		
		textInPos30 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textInPos30.setText("0.0");
		textInPos30.setBounds(606, 296, 39, 31);
		
		textInPos01 = new Text(shell, SWT.BORDER);
		textInPos01.setBounds(651, 185, 39, 31);
		
		textInPos11 = new Text(shell, SWT.BORDER);
		textInPos11.setBounds(651, 222, 39, 31);
		
		textInPos21 = new Text(shell, SWT.BORDER);
		textInPos21.setBounds(651, 259, 39, 31);
		
		textInPos31 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textInPos31.setText("0.0");
		textInPos31.setBounds(651, 296, 39, 31);
		
		textInPos02 = new Text(shell, SWT.BORDER);
		textInPos02.setBounds(696, 185, 39, 31);
		
		textInPos12 = new Text(shell, SWT.BORDER);
		textInPos12.setBounds(696, 222, 39, 31);
		
		textInPos22 = new Text(shell, SWT.BORDER);
		textInPos22.setBounds(696, 259, 39, 31);
		
		textInPos32 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textInPos32.setText("0.0");
		textInPos32.setBounds(696, 296, 39, 31);
		
		textInPos03 = new Text(shell, SWT.BORDER);
		textInPos03.setBounds(741, 185, 39, 31);
		
		textInPos13 = new Text(shell, SWT.BORDER);
		textInPos13.setBounds(741, 222, 39, 31);
		
		textInPos23 = new Text(shell, SWT.BORDER);
		textInPos23.setBounds(741, 259, 39, 31);
		
		textInPos33 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textInPos33.setText("1.0");
		textInPos33.setBounds(741, 296, 39, 31);
		
		textFinPos00 = new Text(shell, SWT.BORDER);
		textFinPos00.setBounds(606, 417, 39, 31);
		
		textFinPos01 = new Text(shell, SWT.BORDER);
		textFinPos01.setBounds(651, 417, 39, 31);
		
		textFinPos02 = new Text(shell, SWT.BORDER);
		textFinPos02.setBounds(696, 417, 39, 31);
		
		textFinPos03 = new Text(shell, SWT.BORDER);
		textFinPos03.setBounds(741, 417, 39, 31);
		
		textFinPos13 = new Text(shell, SWT.BORDER);
		textFinPos13.setBounds(741, 454, 39, 31);
		
		textFinPos12 = new Text(shell, SWT.BORDER);
		textFinPos12.setBounds(696, 454, 39, 31);
		
		textFinPos11 = new Text(shell, SWT.BORDER);
		textFinPos11.setBounds(651, 454, 39, 31);
		
		textFinPos10 = new Text(shell, SWT.BORDER);
		textFinPos10.setBounds(606, 454, 39, 31);
		
		textFinPos20 = new Text(shell, SWT.BORDER);
		textFinPos20.setBounds(606, 491, 39, 31);
		
		textFinPos30 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textFinPos30.setText("0.0");
		textFinPos30.setBounds(606, 528, 39, 31);
		
		textFinPos31 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textFinPos31.setText("0.0");
		textFinPos31.setBounds(651, 528, 39, 31);
		
		textFinPos21 = new Text(shell, SWT.BORDER);
		textFinPos21.setBounds(651, 491, 39, 31);
		
		textFinPos22 = new Text(shell, SWT.BORDER);
		textFinPos22.setBounds(696, 491, 39, 31);
		
		textFinPos32 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textFinPos32.setText("0.0");
		textFinPos32.setBounds(696, 528, 39, 31);
		
		textFinPos33 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textFinPos33.setText("1.0");
		textFinPos33.setBounds(741, 528, 39, 31);
		
		textFinPos23 = new Text(shell, SWT.BORDER);
		textFinPos23.setBounds(741, 491, 39, 31);
		
		textMax0 = new Text(shell, SWT.BORDER);
		textMax0.setBounds(198, 354, 39, 31);
		
		textMax1 = new Text(shell, SWT.BORDER);
		textMax1.setBounds(243, 354, 39, 31);
		
		textMax2 = new Text(shell, SWT.BORDER);
		textMax2.setBounds(288, 354, 39, 31);
		
		textMin0 = new Text(shell, SWT.BORDER);
		textMin0.setBounds(198, 224, 39, 31);
		
		textMin1 = new Text(shell, SWT.BORDER);
		textMin1.setBounds(243, 224, 39, 31);
		
		textMin2 = new Text(shell, SWT.BORDER);
		textMin2.setBounds(288, 224, 39, 31);
		
		textPos230 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos230.setEnabled(false);
		textPos230.setText("0.0");
		textPos230.setBounds(834, 296, 39, 31);
		
		textPos231 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos231.setEnabled(false);
		textPos231.setText("0.0");
		textPos231.setBounds(879, 296, 39, 31);
		
		textPos232 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos232.setEnabled(false);
		textPos232.setText("0.0");
		textPos232.setBounds(924, 296, 39, 31);
		
		textPos233 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos233.setEnabled(false);
		textPos233.setText("1.0");
		textPos233.setBounds(969, 296, 39, 31);
		
		textPos223 = new Text(shell, SWT.BORDER);
		textPos223.setEnabled(false);
		textPos223.setBounds(969, 259, 39, 31);
		
		textPos222 = new Text(shell, SWT.BORDER);
		textPos222.setEnabled(false);
		textPos222.setBounds(924, 259, 39, 31);
		
		textPos221 = new Text(shell, SWT.BORDER);
		textPos221.setEnabled(false);
		textPos221.setBounds(879, 259, 39, 31);
		
		textPos220 = new Text(shell, SWT.BORDER);
		textPos220.setEnabled(false);
		textPos220.setBounds(834, 259, 39, 31);
		
		textPos210 = new Text(shell, SWT.BORDER);
		textPos210.setEnabled(false);
		textPos210.setBounds(834, 222, 39, 31);
		
		textPos211 = new Text(shell, SWT.BORDER);
		textPos211.setEnabled(false);
		textPos211.setBounds(879, 222, 39, 31);
		
		textPos212 = new Text(shell, SWT.BORDER);
		textPos212.setEnabled(false);
		textPos212.setBounds(924, 222, 39, 31);
		
		textPos213 = new Text(shell, SWT.BORDER);
		textPos213.setEnabled(false);
		textPos213.setBounds(969, 222, 39, 31);
		
		textPos203 = new Text(shell, SWT.BORDER);
		textPos203.setEnabled(false);
		textPos203.setBounds(969, 185, 39, 31);
		
		textPos202 = new Text(shell, SWT.BORDER);
		textPos202.setEnabled(false);
		textPos202.setBounds(924, 185, 39, 31);
		
		textPos201 = new Text(shell, SWT.BORDER);
		textPos201.setEnabled(false);
		textPos201.setBounds(879, 185, 39, 31);
		
		textPos200 = new Text(shell, SWT.BORDER);
		textPos200.setEnabled(false);
		textPos200.setBounds(834, 185, 39, 31);
		
		lblndPosition = new Label(shell, SWT.NONE);
		lblndPosition.setEnabled(false);
		lblndPosition.setText("2nd Position");
		lblndPosition.setBounds(1047, 225, 116, 25);
		
		lblrdPosition = new Label(shell, SWT.NONE);
		lblrdPosition.setEnabled(false);
		lblrdPosition.setText("3rd Position");
		lblrdPosition.setBounds(1047, 457, 143, 25);
		
		textPos300 = new Text(shell, SWT.BORDER);
		textPos300.setEnabled(false);
		textPos300.setBounds(834, 417, 39, 31);
		
		textPos310 = new Text(shell, SWT.BORDER);
		textPos310.setEnabled(false);
		textPos310.setBounds(834, 454, 39, 31);
		
		textPos320 = new Text(shell, SWT.BORDER);
		textPos320.setEnabled(false);
		textPos320.setBounds(834, 491, 39, 31);
		
		textPos330 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos330.setEnabled(false);
		textPos330.setText("0.0");
		textPos330.setBounds(834, 528, 39, 31);
		
		textPos331 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos331.setEnabled(false);
		textPos331.setText("0.0");
		textPos331.setBounds(879, 528, 39, 31);
		
		textPos321 = new Text(shell, SWT.BORDER);
		textPos321.setEnabled(false);
		textPos321.setBounds(879, 491, 39, 31);
		
		textPos311 = new Text(shell, SWT.BORDER);
		textPos311.setEnabled(false);
		textPos311.setBounds(879, 454, 39, 31);
		
		textPos301 = new Text(shell, SWT.BORDER);
		textPos301.setEnabled(false);
		textPos301.setBounds(879, 417, 39, 31);
		
		textPos302 = new Text(shell, SWT.BORDER);
		textPos302.setEnabled(false);
		textPos302.setBounds(924, 417, 39, 31);
		
		textPos303 = new Text(shell, SWT.BORDER);
		textPos303.setEnabled(false);
		textPos303.setBounds(969, 417, 39, 31);
		
		textPos313 = new Text(shell, SWT.BORDER);
		textPos313.setEnabled(false);
		textPos313.setBounds(969, 454, 39, 31);
		
		textPos312 = new Text(shell, SWT.BORDER);
		textPos312.setEnabled(false);
		textPos312.setBounds(924, 454, 39, 31);
		
		textPos322 = new Text(shell, SWT.BORDER);
		textPos322.setEnabled(false);
		textPos322.setBounds(924, 491, 39, 31);
		
		textPos323 = new Text(shell, SWT.BORDER);
		textPos323.setEnabled(false);
		textPos323.setBounds(969, 491, 39, 31);
		
		textPos333 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos333.setEnabled(false);
		textPos333.setText("1.0");
		textPos333.setBounds(969, 528, 39, 31);
		
		textPos332 = new Text(shell, SWT.BORDER | SWT.READ_ONLY);
		textPos332.setEnabled(false);
		textPos332.setText("0.0");
		textPos332.setBounds(924, 528, 39, 31);
		
		textVelocity = new Text(shell, SWT.BORDER);
		textVelocity.setBounds(928, 58, 39, 31);
		
		textAcceleration = new Text(shell, SWT.BORDER);
		textAcceleration.setBounds(928, 119, 39, 31);
		
		Label lblNewLabel = new Label(shell, SWT.NONE);
		lblNewLabel.setBounds(817, 61, 81, 25);
		lblNewLabel.setText("Velocity");
		
		Label lblNewLabel_1 = new Label(shell, SWT.NONE);
		lblNewLabel_1.setBounds(804, 122, 105, 25);
		lblNewLabel_1.setText("Acceleration");
		
		lblNewLabel_2 = new Label(shell, SWT.NONE);
		lblNewLabel_2.setBounds(1001, 64, 159, 25);
		lblNewLabel_2.setText("( default: 0.2 [m/s] )");
		
		lblNewLabel_3 = new Label(shell, SWT.NONE);
		lblNewLabel_3.setBounds(1001, 125, 162, 25);
		lblNewLabel_3.setText("( default: 1 [m/s^2] )");

	}
}
