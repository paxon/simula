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

    private double[][] position;
    private double[][] velocity;
    private double[][] acceleration;

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

    public boolean getParamsFromConsole()
    {
        double[] size = new double[DIMENSIONS];
        int partQty = 0;
        Scanner sc = new Scanner(System.in);
        System.out.println("Please enter size of cell:");
        for (int i=0; i<DIMENSIONS; i++) size[i] = sc.nextDouble();
        sc.nextLine();
        System.out.println("Enter quantity of particles:");
        partQty = sc.nextInt();
        sc.nextLine();
        return init(size, partQty);
    }

    private boolean init (double[] size, int partQty) {
        this.size = size;
        this.partQty = partQty;
        position = new double[partQty][DIMENSIONS];
        velocity = new double[partQty][DIMENSIONS];
        acceleration = new double[partQty][DIMENSIONS];
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
                        position[partNum][0] = (num % axisFlat[0]) * distance;
                        for (int i=1; i<DIMENSIONS; i++){
                        position[partNum][i] = ((num / axisFlat[i-1]) % axisFlat[i]) * distance;
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
                velocity[i][j] = vel;
                generalVelocity[j] += vel;
            }

            generalVelocity[j] /= partQty;

            for (int i=0; i<partQty; i++)
            {
                double vel = velocity[i][j];
                vel -= generalVelocity[j];
                velocity[i][j] = vel;
                velocitySqrSum = vel*vel;
            }
        }

        double meanParticleKineticEnergy = velocitySqrSum/(2*partQty);
        double velocityAdjustment = Math.sqrt(initKinetEnergy/meanParticleKineticEnergy);

        for (int j=0; j<DIMENSIONS; j++)
        {
            for (int i=0; i<partQty; i++)
            {
                velocity[i][j] *= velocityAdjustment;
                acceleration[i][j] = 0.0;
            }
        }
    }

    public void listParticlesCoord(){
        for(int i=0; i<partQty; i++)
        {
            System.out.print("Prtcl #"+(i+1)+": ");
            for (int j=0; j<DIMENSIONS; j++){
                System.out.print("x" + (j+1) + ": " + position[i][j] + " ");
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
                position[partNum][dim] = velocity[partNum][dim] * dtSqrdByTwo;
                velocity[partNum][dim] = acceleration[partNum][dim] * dtByTwo;

                calc();

                velocity[partNum][dim] = acceleration[partNum][dim] * dtByTwo;
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
                    vector[dim] = periodicDistanceAdjust(position[partOne][dim]-position[partTwo][dim], size[dim]);
                }
                double sqrdDistance = getSqrdDistance(vector);
                double force = PotentialsForces.forceLJ(sqrdDistance);
                for (int dim = 0; dim<DIMENSIONS; dim++)
                {
                    forceVector[dim] = force*vector[dim];
                    acceleration[partOne][dim] += forceVector[dim];
                    acceleration[partTwo][dim] -= forceVector[dim];
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

    private double periodicPositionAdjust(double distance, double axisPeriod)
    {
        if (distance > axisPeriod) distance -= axisPeriod;
        else if (distance < 0) distance += axisPeriod;
        return distance;
    }

    private double periodicDistanceAdjust(double distance, double axisPeriod)
    {
        if (distance > axisPeriod/2.0) distance -= axisPeriod;
        else if (distance < -axisPeriod/2.0) distance += axisPeriod;
        return distance;
    }
}
