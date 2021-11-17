/**
 *
 * @author  Javier
 * @author  Nicolas
 * @author  Andy
 * @author  Harvy
 * @author  Mileidy
 * @author  Maria Alejandra
 * @author  Juan David
 */

public class Point implements Cloneable {
	public float[] coords;
	
	public Point()
	{
		coords = new float[2];
	}
	public Point(float x, float y)
	{
		coords = new float[2];
		coords[0]=x;
		coords[1]=y;
	}

	public Point clone()
	{
		Point pointret = new Point(coords[0],coords[1]);
		return pointret;
	}
	
	public float getX() {
		return coords[0];
	}
	
	public float getY() {
		return coords[1];
	}
	public float distance(float x, float y)
	{
		return (float) Math.sqrt(Math.pow(x-coords[0], 2)+Math.pow(y-coords[1], 2));
	}
}