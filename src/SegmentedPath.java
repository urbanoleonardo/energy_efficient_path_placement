import java.util.List;

public class SegmentedPath extends Path{
	private List<Target> points = null;
	
	public SegmentedPath(){
		super();
	}
	
	public SegmentedPath(List<Target> points){
		super(points.get(0), points.get(points.size() - 1));
		this.points.addAll(points);
	}
	
	public SegmentedPath(List<Target> points, double tSample, double xSample){
		super(points.get(0), points.get(points.size() - 1), tSample, xSample );
		this.points.addAll(points);
	}
	
	public SegmentedPath(Target[] points){
		super(points[0], points[points.length - 1]);
		for(int i = 0; i < points.length; i++){
			this.points.add(points[i]);
		}
	}
	
	public SegmentedPath(Target[] points, double tSample, double xSample){
		super(points[0], points[points.length - 1], tSample, xSample);
		for(int i = 0; i < points.length; i++){
			this.points.add(points[i]);
		}
	}

	public void setPoints (List<Target> points){
		this.points.clear();
		this.points.addAll(points);
	}
	
	public void addPoints (List<Target> points){
		this.points.addAll(points);
		this.finalPosition = this.points.get(this.points.size() - 1);
	}
	
	public void addPoints (int index, List<Target> points){
		this.points.addAll(index, points);
		this.finalPosition = this.points.get(this.points.size() - 1);
	}
	
	@Override
	public Target[] getPathPositions(){
		Target[] pathPositions = new Target[this.points.size()];
		for(int i = 0; i < pathPositions.length; i++){
			pathPositions[i] = this.points.get(i);
		}
		
		return pathPositions;
	}
	
	@Override 
	public void translatePath(double[] vector){
		for(Target point : points){
			point.translateTarget(vector);
		}
		
		if(this.interpolatedPath != null){
			this.translateTrajectory(vector);
		}
	}
	
	/**
	 * This method is interpolating through the points in the list.
	 *  
 	 * @return the trajectory that will be assigned to the Path (Parent) attribute.
	 */
	
	public Trajectory interpolate(){
		Trajectory interpolatedPath = new Trajectory();
		
		for(int i = 0; i < this.points.size() - 1; i++){
			Path tempPath = new Path(points.get(i), points.get(i+1), this.tSample, this.xSample);
			tempPath.setMaxAcc(super.getMaxAcc());
			tempPath.setMaxVel(super.getMaxVel());
			Trajectory interpolatedSide = tempPath.interpolate();
			
			if(i != 0){
				double time = interpolatedPath.timeInstants.get(interpolatedPath.timeInstants.size() - 1);
				
				interpolatedSide.points.remove(0);
				interpolatedSide.timeInstants.remove(0);
				interpolatedSide.shiftTime(time);
				
			}
			
			//add the lists of points, times, forces and torques to the complete Trajectory
			interpolatedPath.points.addAll(interpolatedSide.points);
			interpolatedPath.timeInstants.addAll(interpolatedSide.timeInstants);
			interpolatedPath.extForces.addAll(interpolatedSide.extForces);
			interpolatedPath.extTorques.addAll(interpolatedSide.extTorques);
		}
		
		this.interpolatedPath = interpolatedPath;
		return interpolatedPath;
	}
	
}
