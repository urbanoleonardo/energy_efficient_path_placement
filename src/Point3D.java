
public class Point3D {
	
	private final double x;
	private final double y;
	private final double z;
	
	public Point3D(double x, double y, double z){
		
		this.x = x;
		this.y = y;
		this.z = z;
		
	}
	
	public Point3D(double[] pos){
		
		if(pos == null || pos.length != 3)
			throw new IllegalArgumentException();
		
		this.x = pos[0];
		this.y = pos[1];
		this.z = pos[2];
		
	}

	public double getX() {
		return x;
	}

	public double getY() {
		return y;
	}

	public double getZ() {
		return z;
	}
	
	public double[] getPosition() {
		
		return new double[] {x, y, z};
		
	}

}
