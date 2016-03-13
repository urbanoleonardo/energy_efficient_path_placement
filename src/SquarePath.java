
public class SquarePath extends Path{
	/*
	 * This class is used when a square path is given in input, that is with a total of 4 points in the space.
	 * The trajectory will go along all 4 points and finish back to the start (initial position).
	 * The interpolation method o SquarePath should call 4 times Path's interpolate().
	 */
	
	private Target vertex1;
	private Target vertex2;
	
	public SquarePath(){
		super();
	}
	
	public SquarePath(Target[] squarePoints){
		super(squarePoints[0], squarePoints[3]);
		this.vertex1 = squarePoints[1];
		this.vertex2 = squarePoints[2];
	}
	
	public SquarePath(Target[] squarePoints, double tSample, double xSample){
		super(squarePoints[0], squarePoints[3], tSample, xSample);
		this.vertex1 = squarePoints[1];
		this.vertex2 = squarePoints[2];
	}
	
	public SquarePath(Target point1, Target point2, Target point3, Target point4,  double tSample, double xSample){
		super(point1, point2, tSample, xSample);
		this.vertex1 = point3;
		this.vertex2 = point4;
	}
	
	public SquarePath(Target point1, Target point2, Target point3, Target point4){
		super(point1, point2);
		this.vertex1 = point3;
		this.vertex2 = point4;
	}
	
	public Target getVertex1() {
		return vertex1;
	}

	public void setVertex1(Target vertex1) {
		this.vertex1 = vertex1;
	}

	public Target getVertex2() {
		return vertex2;
	}

	public void setVertex2(Target vertex2) {
		this.vertex2 = vertex2;
	}

	public Trajectory interpolate(){
		/*
		 * This method is interpolating through the 4 points specified in the constructor. 
		 * It returns the trajectory that will be the Path (Parent) attribute.
		 * The Parent method interpolate will be called for each couple of points and then merged.
		 */
		Trajectory interpolatedPath = new Trajectory();
		Target[] points = {super.getInitialPosition(), this.vertex1, this.vertex2, super.getFinalPosition(), super.getInitialPosition()};
		
		
		for(int i = 0 ; i < (points.length - 1); i++){
			Path tempPath = new Path(points[i], points[i+1], super.getTimeSample(), super.getSpaceSample());
			tempPath.setMaxAcc(super.getMaxAcc());
			tempPath.setMaxVel(super.getMaxVel());
			Trajectory interpolatedSide = tempPath.interpolate();
			
//			double time = 0;
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
		
		super.interpolatedPath = interpolatedPath;
		return interpolatedPath;
	}
	
	
}
