package simula;

/**
 * Created by sash on 29.09.16.
 */
public class PotentialsForces {
    public static double potenLJ(double rSqrd) {
        double rTwo = 1/rSqrd;
        double rSix = rTwo*rTwo*rTwo;
        return 4*(rSix*rSix-rSix);
    }

    public static double forceLJ(double rSqrd) {
        double ro2 = 1/rSqrd;
        double ro4 = ro2*ro2;
        double ro8 = ro4*ro4;
        return 48*ro8*(ro4*ro2-0.5);
    }
}
