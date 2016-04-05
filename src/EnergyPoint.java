import javax.vecmath.Color3f;

public class EnergyPoint implements Comparable<EnergyPoint>{

	private double x;
	private double y;
	private double z;
	
	private double energy;
	private Color3f color;
	
	public EnergyPoint () {}
	
	public EnergyPoint (double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public EnergyPoint (double energy){
		this.energy = energy;
	}
	
	public double[] getPosition() {
		return new double[] {x, y, z};
	}
	
	public void setPosition(double[] pos) {
		
		x = pos[0];
		y = pos[1];
		z = pos[2];
		
	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	public double getEnergy() {
		return energy;
	}

	public void setEnergy(double energy) {
		this.energy = energy;
	}

	public Color3f getColor() {
		return color;
	}

	public void setColor(Color3f color) {
		
		if(color == null)
			throw new IllegalArgumentException(x + " " + y + " " + z);
		
		this.color = color;
	}

	@Override
	public int compareTo(EnergyPoint e) {

		if(this.energy < e.energy)
			return -1;
		
		if(this.energy > e.energy)
			return 1;
		
		return 0;
		
	}
}

