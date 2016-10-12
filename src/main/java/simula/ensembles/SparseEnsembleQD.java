package simula.ensembles;

import simula.routines.PotentialsForces;

import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;

import static simula.routines.SpatialRoutines.*;
import static simula.routines.SpatialRoutines.getSquaredDistance;

/**
 * Created by HP on 11.10.2016.
 */
public class SparseEnsembleQD implements Runnable{


    private final int DIMENSIONS = 3;
    private final double MINDISTANCE = 1.05;
    private final double TIMESTEP = 0.01;
    private final double INITKE = 1.0;


    private final double dt = TIMESTEP;
    private final double dtByTwo = dt /2.0;
    private final double dtSqrdByTwo = dtByTwo * dt;

    private double threshold = 3;
    private double thresholdSqrd = threshold*threshold;
    private double[] size = new double[DIMENSIONS];
    private int partQty = 0;

    private double[][] pos;
    private double[][] vel;
    private double[][] acc;

    private double potentialEnergyAccumulator = 0.0;
    private double virialTotal = 0.0;
    private double kineticEnergyAccumulator = 0.0;
    private double kineticEnergySquaredAccumulator = 0.0;
    private int steps = 0;



// gesture

    public boolean getParamsFromConsole()
    {
        Scanner sc = new Scanner(System.in);
        System.out.println("Please enter size of cell:");
        for (int i=0; i<size.length; i++) size[i] = sc.nextDouble();
        sc.nextLine();

        System.out.println("Enter quantity of particles:");
        partQty = sc.nextInt();
        sc.nextLine();

        return init(size, partQty);
    }

    public boolean init (double[] size, int partQty)
    {
        this.size = size;
        this.partQty = partQty;

        pos = new double[partQty][DIMENSIONS];
        vel = new double[partQty][DIMENSIONS];
        acc = new double[partQty][DIMENSIONS];

        return rectAlign(MINDISTANCE, size, partQty);
    }


    private boolean rectAlign (double distance, double[] size, int partQty)
    {
        /**
         * Particles alignment routine.
         * Checks whether there is enough space to lay particles upon rectangular lattice. Aligns particles spaced by
         * var distance. If there is not enough particles to fill the mesh some nodes are skipped in random order.
         */
        int[] nodesPerAxis = new int[size.length];
        int[] posTuple = new int[nodesPerAxis.length];
        int latticeCapacity = 1;
        for (int i=0; i<nodesPerAxis.length; i++) {
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
                        posTuple = getTupleFromIndex(num,posTuple,nodesPerAxis);
                        for (int i=0; i<pos[partNum].length; i++){
                            pos[partNum][i] = posTuple[i] * distance;
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

        /**
         * Particles are supplied with initial kinetic energy. Random generated velocities are normalized in order to
         * fit declared kinetic energy.
         */
        double velocity;
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

    public void listSystemState () {
        System.out.print("Mean Temp: ");
        System.out.printf("%.4f", getMeanTemperature());
        System.out.print(" Mean Pressure: ");
        System.out.printf("%.4f", getMeanPressure());
        System.out.print(" Mean Energy: ");
        System.out.printf("%.4f", getMeanEnergy());
        System.out.print(" ");
    }


    public void run()
    {
        initParticles(INITKE,-0.5,0.5);
        calculateForces();
        doSteps(10,10000);
    }

    private void step()
    {
        /**
         * Simulation step use simplified form of Verlet method of dynamics calculation. Velocity change is split to
         * actual and previous acceleration.
         */
        double kineticEnergyTotal = 0;
        for (int partNum=0; partNum<partQty; partNum++)
        {
            for (int dim = 0; dim < DIMENSIONS; dim++)
            {
                pos[partNum][dim] = periodicPositionShift(pos[partNum][dim]+vel[partNum][dim]*dt+acc[partNum][dim]*dtSqrdByTwo, size[dim]);
                vel[partNum][dim] += acc[partNum][dim] * dtByTwo;
            }
        }

        calculateForces();

        for (int partNum=0; partNum<partQty; partNum++) for (int dim=0; dim<DIMENSIONS; dim++)
        {
            vel[partNum][dim] += acc[partNum][dim] * dtByTwo;
            kineticEnergyTotal += getSquaredVectorLength(vel[partNum]);
        }

        kineticEnergyTotal /= 2.0;
        kineticEnergyAccumulator += kineticEnergyTotal;
        kineticEnergySquaredAccumulator += kineticEnergyTotal*kineticEnergyTotal;
        steps++;
    }

    private void doSteps(int stepsPerRun, int runsToDo)
    {
        Long nanoT;


        for(int i=0; i<runsToDo; i++)
        {
            nanoT = System.nanoTime();
            for (int j=0; j<stepsPerRun; j++) { calculateForces(); step();   }

            listSystemState();
            System.out.println((System.nanoTime()-nanoT)/stepsPerRun + "ns per step (QD)");
        }
    }

    private void calculateForces()
    {
        double[] distance = new double[DIMENSIONS];
        double[] forceVector = new double[DIMENSIONS];
        double force = 0;
        for (int i=0; i<acc.length; i++) for (int j=0; j<acc[i].length; j++) acc[i][j] = 0.0;
        for (int p1=0; p1<partQty-1; p1++)
        {
            for (int p2=p1+1; p2<partQty; p2++)
            {
                for (int dim = 0; dim<distance.length; dim++)
                {
                    distance[dim] = periodicDistanceShift(pos[p1][dim]-pos[p2][dim], size[dim]);
                }
                double sqrdDistance = getSquaredDistance(distance);
                if (sqrdDistance<thresholdSqrd)
                {
                    force = PotentialsForces.forceLJ(sqrdDistance);

                    for (int dim = 0; dim<distance.length; dim++)
                    {
                        forceVector[dim] = force*distance[dim];
                        acc[p1][dim] += forceVector[dim];
                        acc[p2][dim] -= forceVector[dim];
                        virialTotal += distance[dim]*forceVector[dim];
                    }
                    potentialEnergyAccumulator += PotentialsForces.potenLJ(sqrdDistance);
                }
            }
        }
    }

    public double getMeanTemperature() {
        return kineticEnergyAccumulator/(partQty*steps);
    }
    public double getMeanEnergy() { return (kineticEnergyAccumulator+potentialEnergyAccumulator)/steps; }
    public double getMeanPressure() { return 1.0+0.5*virialTotal/(steps*partQty*getMeanTemperature());}
}
