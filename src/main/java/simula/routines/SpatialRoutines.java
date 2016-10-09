package simula.routines;

/**
 * Created by HP on 09.10.2016.
 */
public class SpatialRoutines {

    public static int[] getTupleFromIndex(int index, int[] tuple, int[] tupleSize)
    {
        if (tupleSize.length>0 && (tuple.length >= tupleSize.length))
        {
            tuple[0] = index % tupleSize[0];

            if (tupleSize.length>1) {
                int m = tupleSize[0];
                for (int i=1; i<tupleSize.length-1; i++)
                {
                    tuple[i] = (index / m) % tupleSize[i];
                    m *= tupleSize[i];
                }
                tuple[tupleSize.length-1]= index / m;
            }
        }
        return tuple;
    }

    public static int[] getTupleFromIndex(int index, int[] tupleSize)
    {
        if (tupleSize.length>0)
        {
            int[] tuple = new int[tupleSize.length];
            tuple[0] = index % tupleSize[0];

            if (tupleSize.length>1) {
                int m = tupleSize[0];
                for (int i=1; i<tupleSize.length-1; i++)
                {
                    tuple[i] = (index / m) % tupleSize[i];
                    m *= tupleSize[i];
                }
                tuple[tupleSize.length-1]= index / m;
            }
            return tuple;
        } else {
            int[] nulltuple = {};
            return nulltuple;
        }
    }
    public static double getSquaredDistance(double[] vector)
    {
        double res = 0;
        for (int dim=0; dim<vector.length; dim++)
        {
            double distance = vector[dim];
            res += distance*distance;
        }
        return res;
    }

    public static double getSquaredVectorLength(double[] vector)
    {
        double res = 0;
        for (int dim=0; dim<vector.length; dim++)
        {
            double distance = vector[dim];
            res += distance*distance;
        }
        return res;
    }

    public static double periodicPositionShift(double position, double period)
    {
        while (position>period) position -= period;
        while (position<0) position += period;
        return position;
    }

    public static double periodicDistanceShift(double distance, double period)
    {
        while (distance > period*0.5) distance -= period;
        while (distance < -period*0.5) distance += period;
        return distance;
    }
}
