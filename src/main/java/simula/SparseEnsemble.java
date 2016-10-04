package simula;

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

    private double[][] pos;
    private double[][] vel;
    private double[][] acc;

    private double potentialEnergyAccumulator = 0.0;
    private double virialTotal = 0.0;
    private double kineticEnergyAccumulator = 0.0;
    private double kineticEnergySquaredAccumulator = 0.0;
    private int steps = 0;


    public boolean getParamsFromConsole()
    {
        double[] size = new double[DIMENSIONS];
        int partQty = 0;

        Scanner sc = new Scanner(System.in);
        System.out.println("Please enter size of cell:");
        for (int i=0; i<size.length; i++) size[i] = sc.nextDouble();
        sc.nextLine();

        System.out.println("Enter quantity of particles:");
        partQty = sc.nextInt();
        sc.nextLine();

        return init(size, partQty);
    }

    private boolean init (double[] size, int partQty)
    {
        this.size = size;
        this.partQty = partQty;

        pos = new double[partQty][DIMENSIONS];
        vel = new double[partQty][DIMENSIONS];
        acc = new double[partQty][DIMENSIONS];

        return rectAlign(MINDISTANCE);
    }


    private boolean rectAlign (double distance)
    {
        /**
         * Checks whether there is enough space to lay particles upon rectangular lattice. Aligns particles spaced by
         * var distance. If there is not enough particles to fill the mesh some nodes are skipped in random order.
         */
        int[] nodesPerAxis = new int[DIMENSIONS];
        int latticeCapacity = 1;
        for (int i=0; i<DIMENSIONS; i++) {
            int intervals = (int) (size[i] / distance);
            nodesPerAxis[i] = intervals;
            latticeCapacity *= intervals;
        }

        if (latticeCapacity>=partQty) {
            boolean[] lattice = new boolean[latticeCapacity];
            for (int partNum = 0; partNum<partQty; partNum++)
            {
                int num;
                boolean undone = true;
                while (undone) {
                    num = ThreadLocalRandom.current().nextInt(latticeCapacity);
                    if (!lattice[num])
                    {
                        lattice[num] = true;
                        pos[partNum][0] = ( (num+1) % nodesPerAxis[0]) * distance;
                        for (int i=1; i<DIMENSIONS; i++){
                        pos[partNum][i] = (((num+1) / nodesPerAxis[i-1]) % nodesPerAxis[i]) * distance;
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

    public void initParticles(double initKinetEnergy, double minVelocity, double maxVelocity)
    {
        double velocity = 0;
        double[] generalVelocity = new double[DIMENSIONS];
        double velocitySqrSum = 0.0;
        for (int j=0; j<DIMENSIONS; j++)
        {
            generalVelocity[j] = 0.0;

            for (int i=0; i<partQty; i++)
            {
                velocity = ThreadLocalRandom.current().nextDouble(minVelocity,maxVelocity);
                vel[i][j] = velocity;
                generalVelocity[j] += velocity;
            }

            generalVelocity[j] /= partQty;

            for (int i=0; i<partQty; i++)
            {
                velocity = vel[i][j];
                velocity -= generalVelocity[j];
                vel[i][j] = velocity;
                velocitySqrSum += velocity*velocity;
            }
        }

        double particialKineticEnergy = velocitySqrSum/(2*partQty);
        double energyScalingRatio = Math.sqrt(initKinetEnergy/particialKineticEnergy);

        for (int j=0; j<DIMENSIONS; j++)
        {
            for (int i=0; i<partQty; i++)
            {
                vel[i][j] *= energyScalingRatio;
                acc[i][j] = 0.0;
            }
        }
    }

    public void listParticlesCoord(){
        for(int i=0; i<partQty; i++)
        {
            System.out.print("Prtcl #"+(i+1)+": ");
            for (int j=0; j<DIMENSIONS; j++){
                System.out.print("x" + (j+1) + ": " + pos[i][j] + " ");
            }
            System.out.println();
        }
    }


    public void start()
    {
        initParticles(0,-0.5,0.5);
        calculateForces();
        doSteps(1,10);
    }

    private void step()
    {
        double kineticEnergyTotal = 0;
        for (int partNum=0; partNum<partQty; partNum++)
        {
            for (int dim = 0; dim < DIMENSIONS; dim++)
            {
                pos[partNum][dim] = periodicPositionShift(vel[partNum][dim]*dtSqrdByTwo, size[dim]);
                vel[partNum][dim] = acc[partNum][dim] * dtByTwo;
            }
        }

        calculateForces();

        for (int partNum=0; partNum<partQty; partNum++)
        {
            for (int dim=0; dim<DIMENSIONS; dim++)
            {
                vel[partNum][dim] = acc[partNum][dim]* dtByTwo;
                kineticEnergyTotal += getSquaredVectorLength(vel[partNum]);
            }
        }
        kineticEnergyTotal /= 2.0;
        kineticEnergyAccumulator += kineticEnergyTotal;
        kineticEnergySquaredAccumulator += kineticEnergyTotal*kineticEnergyTotal;

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
            System.out.println(potentialEnergyAccumulator);
        }
    }

    public void calculateForces()
    {
        double[] distance = new double[DIMENSIONS];
        double[] forceVector = new double[DIMENSIONS];
        for (int p1=0; p1<partQty-1; p1++)
        {
            for (int p2=p1+1; p2<partQty; p2++)
            {
                for (int dim = 0; dim<DIMENSIONS; dim++)
                {
                    distance[dim] = periodicDistanceShift(pos[p1][dim]-pos[p2][dim], size[dim]);
                }
                double sqrdDistance = getSquaredDistance(distance);
                double force = PotentialsForces.forceLJ(sqrdDistance);
                for (int dim = 0; dim<DIMENSIONS; dim++)
                {
                    forceVector[dim] = force*distance[dim];
                    acc[p1][dim] += forceVector[dim];
                    acc[p2][dim] -= forceVector[dim];
                    virialTotal += distance[dim]*forceVector[dim];
                }
                potentialEnergyAccumulator += PotentialsForces.potenLJ(sqrdDistance);
            }
        }
        System.out.println(potentialEnergyAccumulator);
    }

    public double getSquaredDistance(double[] vector)
    {
        double res = 0;
        for (int dim=0; dim<DIMENSIONS; dim++)
        {
            double distance = periodicDistanceShift(vector[dim], size[dim]);
            res += distance*distance;
        }
        return res;
    }

    private double getSquaredVectorLength(double[] vector)
    {
        double res = 0;
        for (int dim=0; dim<DIMENSIONS; dim++)
        {
            double distance = vector[dim];
            res += distance*distance;
        }
        return res;
    }

    public double periodicPositionShift(double position, double axisPeriod)
    {
        if (position > 0) while (position>axisPeriod) position -= axisPeriod;
        else while (position<0) position += axisPeriod;
        return position;
    }

    public double periodicDistanceShift(double distance, double axisPeriod)
    {
        if (distance > 0) while (distance > axisPeriod/2.0) distance -= axisPeriod;
        else while (distance < -axisPeriod/2.0) distance += axisPeriod;
        return distance;
    }


}
