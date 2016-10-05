package simula;
//gesture

import com.sun.org.apache.xerces.internal.impl.xpath.regex.RegularExpression;

import java.util.concurrent.ThreadLocalRandom;


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

    public interface Tupling {
        String tellName();

        void doTuple(int cycles, int bound, int[] tuple, int[] tupleSize);
    }

    public static Tupling[] tuplings = new Tupling[] {
            new Tupling() {
                public String tellName () {return "Unreceived Single Creation";}
                public void doTuple(int cycles, int bound, int[] tuple, int[] tupleSize)
                {
                    indexTuplerUnreceivedSingleCreation(cycles, bound, tuple, tupleSize);
                }
            },
            new Tupling() {
                public String tellName () {return "Received Single Creation";}
                public void doTuple(int cycles, int bound, int[] tuple, int[] tupleSize)
                {
                    indexTuplerReceivedSingleCreation(cycles, bound, tuple, tupleSize);
                }
            },
            new Tupling() {
                public String tellName () {return "Unreceived Multiple Creation";}
                public void doTuple(int cycles, int bound, int[] tuple, int[] tupleSize)
                {
                    indexTuplerUnreceivedMultipleCreation(cycles, bound, tupleSize);
                }
            },
            new Tupling() {
                public String tellName () {return "Received Multiple Creation";}
                public void doTuple(int cycles, int bound, int[] tuple, int[] tupleSize)
                {
                    indexTuplerReceivedMultipleCreation(cycles, bound, tuple, tupleSize);
                }
            },

    };

    private static void indexTuplerUnreceivedSingleCreation(int cycles, int bound, int[] tuple, int[] tupleSize) {
        for (int i=0; i<cycles; i++)
        SparseEnsemble.getTupleFromIndex(ThreadLocalRandom.current().nextInt(bound),tuple,tupleSize);
    }

    private static  void indexTuplerReceivedSingleCreation(int cycles, int bound, int[] tuple, int[] tupleSize) {
        for (int i=0; i<cycles; i++)
            tuple = SparseEnsemble.getTupleFromIndex(ThreadLocalRandom.current().nextInt(bound),tuple,tupleSize);
    }

    private static void indexTuplerUnreceivedMultipleCreation(int cycles, int bound, int[] tupleSize) {
        for (int i=0; i<cycles; i++)
            SparseEnsemble.getTupleFromIndex(ThreadLocalRandom.current().nextInt(bound),tupleSize);
    }

    private static void indexTuplerReceivedMultipleCreation(int cycles, int bound, int[] tuple, int[] tupleSize) {
        for (int i=0; i<cycles; i++)
            tuple = SparseEnsemble.getTupleFromIndex(ThreadLocalRandom.current().nextInt(bound),tupleSize);
    }

    public static void doTuple(int i, int cycles, int bound, int[] tuple, int[] tupleSize)
    {
        tuplings[i].doTuple(cycles, bound, tuple, tupleSize);
    }

    public static void indexTuplerCompare (int warmupRuns, int regularRuns, int cyclesPerCount)
    {
        long nanoT;
        int[] tupleSize = {4,7,11};
        int[] tuple = new int[tupleSize.length];
        int bound = 300;


        for (int h=0; h<4; h++) {
            System.out.println(tuplings[h].tellName());
            for (int i=0; i<warmupRuns+regularRuns; i++)
            {
                if (i<warmupRuns) System.out.print("Warm-up run #" + i);
                else System.out.print("Regular run #");
                nanoT = System.nanoTime();
                doTuple(h,cyclesPerCount,bound,tuple,tupleSize);
                nanoT = System.nanoTime()-nanoT;
                System.out.println(" - " + nanoT/cyclesPerCount + " ns per call");

            }
        }
    }
}
