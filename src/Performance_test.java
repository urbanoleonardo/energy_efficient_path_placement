import java.lang.Math;


public class Performance_test {

	public static int TIMES = 1000000;
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		Test_Math();

	}
	
	public static void Test_Math(){
		
		
		double[] theta_5 = new double[2];
		double[][] R_zyz = {
				{0.45,  0.0,  0.20},
				{0.20, -0.90, 0.10},
				{0.30, -0.10, 1.00}
		};
		long time1 = System.nanoTime();
		for(int i = 0; i < TIMES; i++)
		{
		theta_5[0] = Math.atan2(Math.sqrt(R_zyz[2][0]*R_zyz[2][0] + R_zyz[2][1]*R_zyz[2][1]), R_zyz[2][2]);
		theta_5[1] = Math.atan2(- Math.sqrt(R_zyz[2][0]*R_zyz[2][0] + R_zyz[2][1]*R_zyz[2][1]), R_zyz[2][2]);
		}
		long time2 = System.nanoTime();
		
		long time3 = System.nanoTime();
		for(int i = 0; i < TIMES; i++)
		{
			
		double anglex = Math.sqrt(R_zyz[2][0]*R_zyz[2][0] + R_zyz[2][1]*R_zyz[2][1]);
		theta_5[0] = Math.atan2(anglex, R_zyz[2][2]);
		theta_5[1] = Math.atan2(-anglex, R_zyz[2][2]);
		
		}
		long time4 = System.nanoTime();
		
		System.out.println("Time for the first loop: " + (time2-time1)/1E9);
		System.out.println("Time for the second loop: " + (time4-time3)/1E9);
		
	}

}
