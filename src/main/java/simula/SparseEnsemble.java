package simula;

import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Created by sash on 30.09.16.
 */
public class SparseEnsemble{
    private final int DIMENSIONS = 3;
    private final double MINDISTANCE = 1.05;
    private final double CALCDISTANCE = 2;
    private final double TIMESTEP = 0.01;

    private final double dt = TIMESTEP;
    private final double dtByTwo = dt /2.0;
    private final double dtSqrdByTwo = dtByTwo * dt;

    private double[] size = new double[DIMENSIONS];
    private int partQty;
    private double[][] state;

    private double[][] pos;
    private double[][] vel;
    private double[][] acc;

    private double potenEnergyTotal =0.0;
    private double virialTotal = 0.0;
    private int steps = 0;

    private void setSize(int xi, double size)
    {
        this.size[xi] = size;
    }

    private void setPartQty (int partQty)
    {
        this.partQty = partQty;
    }

    public boolean init () {
        Scanner sc = new Scanner(System.in);
        System.out.println("Please enter size of cell:");
        for (int i=0; i<DIMENSIONS; i++) setSize(i,sc.nextDouble());
        sc.nextLine();
        System.out.println("Enter quantity of particles:");
        setPartQty(sc.nextInt());
        sc.nextLine();
        state = new double[partQty][DIMENSIONS*3];
        return rectAlign(MINDISTANCE);
    }

    private boolean rectAlign (double distance)
    {
        /**
         * Checks whether there is enough space to lay particles upon rectangular lattice. Aligns particles spaced by
         * var distance. If there is not enough particles to fill the mesh some nodes are skipped in random order.
         */
        int[] axisFlat = new int[DIMENSIONS];
        int latticeCapacity = 1;
        for (int i=0; i<DIMENSIONS; i++) {
            int intervals = (int) (size[i] / distance);
            axisFlat[i] = intervals;
            latticeCapacity *= intervals;
        }

        if (latticeCapacity>=partQty) {
            Random nodeRand = new Random();
            boolean[] lattice = new boolean[latticeCapacity];
            for (int partNum = 0; partNum<partQty; partNum++)
            {
                int num;
                boolean undone = true;
                while (undone) {
                    num = nodeRand.nextInt(latticeCapacity);
                    if (!lattice[num])
                    {
                        lattice[num] = true;
                        state[partNum][0] = (num % axisFlat[0]) * distance;
                        for (int i=1; i<DIMENSIONS; i++){
                        state[partNum][i] = ((num / axisFlat[i-1]) % axisFlat[i]) * distance;
                    }
                    undone = false;
                    }
                }
            }
            return true;
        }
        else {
            System.out.println("Not enough nodes to place all particles.");
            return false;
        }
    }



    public void initParticles(double initKinetEnergy, double minVelocity, double maxVelocity) {
        double[] generalVelocity = new double[DIMENSIONS];
        double velocitySqrSum = 0.0;
        for (int j=0; j<DIMENSIONS; j++)
        {
            generalVelocity[j] = 0.0;

            for (int i=0; i<partQty; i++)
            {
                double vel = ThreadLocalRandom.current().nextDouble(minVelocity,maxVelocity);
                state[i][j+DIMENSIONS] = vel;
                generalVelocity[j] += vel;
            }

            generalVelocity[j] /= partQty;

            for (int i=0; i<partQty; i++)
            {
                double vel = state[i][j+DIMENSIONS];
                vel -= generalVelocity[j];
                state[i][j+DIMENSIONS] = vel;
                velocitySqrSum = vel*vel;
            }
        }

        double meanParticleKinetEnergy = velocitySqrSum/(2*partQty);
        double energyAdjCoef = Math.sqrt(initKinetEnergy/meanParticleKinetEnergy);

        for (int j=0; j<DIMENSIONS; j++)
        {
            for (int i=0; i<partQty; i++)
            {
                state[i][j+DIMENSIONS] *= energyAdjCoef;
                state[i][j+DIMENSIONS*2] = 0.0;
            }
        }
    }

    public void listParticlesCoord(){
        for(int i=0; i<partQty; i++)
        {
            System.out.print("Prtcl #"+(i+1)+": ");
            for (int j=0; j<DIMENSIONS; j++){
                System.out.print("x" + (j+1) + ": " + state[i][j] + " ");
            }
            System.out.println();
        }
    }


    public void start()
    {
        initParticles(0,-0.5,0.5);
        calc();
        doSteps(10,1);
    }

    private void step()
    {
        for (int partNum=0; partNum<partQty; partNum++)
        {
            for (int dim=0; dim<DIMENSIONS; dim++)
            {
                state[partNum][dim] = state[partNum][DIMENSIONS+dim]* dtSqrdByTwo;
                state[partNum][DIMENSIONS+dim] =  state[partNum][DIMENSIONS*2+dim]* dtByTwo;

                calc();

                state[partNum][DIMENSIONS+dim] =  state[partNum][DIMENSIONS*2+dim]* dtByTwo;
            }
        }
        steps++;
    }

    private void doSteps(int stepsPerRun, int runsToDo)
    {
        for(int i=0; i<runsToDo; i++)
        {
            for (int j=0; j<stepsPerRun; j++)
            {
                step();
            }
            System.out.println(potenEnergyTotal);
        }
    }

    public void calc()
    {
        double[] vector = new double[DIMENSIONS];
        double[] forceVector = new double[DIMENSIONS];
        for (int partOne=0; partOne<partQty-1; partOne++)
        {
            for (int partTwo=partOne+1; partTwo<partQty; partTwo++)
            {
                for (int dim = 0; dim<DIMENSIONS; dim++)
                {
                    vector[dim] = periodicDistanceAdjust(state[partOne][dim]-state[partTwo][dim], size[dim]);
                }
                double sqrdDistance = getSqrdDistance(vector);
                double force = PotentialsForces.forceLJ(sqrdDistance);
                for (int dim = 0; dim<DIMENSIONS; dim++)
                {
                    forceVector[dim] = force*vector[dim];
                    state[partOne][dim+DIMENSIONS*2] += forceVector[dim];
                    state[partTwo][dim+DIMENSIONS*2] -= forceVector[dim];
                    virialTotal += vector[dim]*forceVector[dim];
                }
                potenEnergyTotal += PotentialsForces.potenLJ(sqrdDistance);

            }
        }
        System.out.println(potenEnergyTotal);
    }

    public double getSqrdDistance(double[] vector)
    {
        double res = 0;
        for (int dim=0; dim<vector.length; dim++)
        {
            double distance = periodicDistanceAdjust(vector[dim], size[dim]);
            res += distance*distance;
        }
        return res;
    }

    public double periodicDistanceAdjust(double distance, double axisPeriod)
    {
        if (distance > axisPeriod/2.0) distance -= axisPeriod;
        else if (distance < -axisPeriod/2.0) distance += axisPeriod;
        return distance;
    }


}
