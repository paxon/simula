package simula;
//gesture
/**
 * Created by HP on 04.10.2016.
 */
public class Tests {
    public static void testLJPotentialValues(int steps, double step)
    {
        double x = 0;
        for (int i=0; i<steps; i++)
        {
            System.out.println(PotentialsForces.potenLJ(x*x));
            x += step;
        }
    }
}
