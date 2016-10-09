package simula.routines;

/**
 * Created by sash on 29.09.16.
 */
public class PotentialsForces {
    public static double potenLJ(double rSqrd) {
        double rTwo = 1.0/rSqrd;
        double rSix = rTwo*rTwo*rTwo;
        return 4.0*rSix*(rSix-1.0);
    }

    public static double forceLJ(double rSqrd) {
        double ro2 = 1.0/rSqrd;
        double ro4 = ro2*ro2;
        double ro8 = ro4*ro4;
        return 48.0*ro8*(ro4*ro2-0.5);
    }
}
