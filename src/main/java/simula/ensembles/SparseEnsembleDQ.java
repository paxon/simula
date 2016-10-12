package simula.ensembles;

import simula.routines.PotentialsForces;

import java.util.Scanner;
import java.util.concurrent.*;

import static simula.routines.SpatialRoutines.*;

/**
 * Created by sash on 30.09.16.
 */
public class SparseEnsembleDQ implements Runnable{

    /**
     * Sparse Ensemble uses Lennard-Jones modelling potential to simulate a small system of interacting particles.
     * MINDISTANCE is a distance where LJ turns to zero (so-called sigma for non-macroscopic for of LJ). System is
     * initialized as set of particles put in nodes of rectangular lattice.
     *
     * This class is used to simulate dynamics of such system by providing macroscopic values expressed in "molecular"
     * units.
     */

    private final int DIMENSIONS = 3;
    private final double MINDISTANCE = 1.05;
    private final double TIMESTEP = 0.01;
    private final double INITKE = 1.0;

    private final double dt = TIMESTEP;
    private final double dtByTwo = dt /2.0;
    private final double dtSqrdByTwo = dtByTwo * dt;

    private double[] size = new double[DIMENSIONS];
    private int partQty = 0;
    private double threshold = 3;
    private double thresholdSqrd = threshold*threshold;

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

        pos = new double[size.length][partQty];
        vel = new double[size.length][partQty];
        acc = new double[size.length][partQty];

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
                        for (int i=0; i<pos.length; i++){
                        pos[i][partNum] = posTuple[i] * distance;
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
        double[] generalVelocity = new double[size.length];
        double velocitySqrSum = 0.0;
        for (int dim=0; dim<size.length; dim++)
        {
            generalVelocity[dim] = 0.0;

            for (int partNum=0; partNum<partQty; partNum++)
            {
                velocity = ThreadLocalRandom.current().nextDouble(minVelocity,maxVelocity);
                vel[dim][partNum] = velocity;
                generalVelocity[dim] += velocity;
            }

            generalVelocity[dim] /= partQty;

            for (int partNum=0; partNum<partQty; partNum++)
            {
                velocity = vel[dim][partNum];
                velocity -= generalVelocity[dim];
                vel[dim][partNum] = velocity;
                velocitySqrSum += velocity*velocity;
            }
        }

        double particialKineticEnergy = velocitySqrSum/(2*partQty);

        double energyScalingRatio = Math.sqrt(initKinetEnergy/particialKineticEnergy);

        for (int dim=0; dim<size.length; dim++)
        {
            for (int partNum=0; partNum<partQty; partNum++)
            {
                vel[dim][partNum] *= energyScalingRatio;
                acc[dim][partNum] = 0.0;
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
        doSteps(1,100);
    }

    private void step()
    {
        /**
         * Simulation step use simplified form of Verlet method of dynamics calculation. Velocity change is split to
         * actual and previous acceleration.
         */
        double kineticEnergyTotal = 0;
        for (int dim = 0; dim < DIMENSIONS; dim++) for (int partNum=0; partNum<partQty; partNum++)
        {
            pos[dim][partNum] = periodicPositionShift(pos[dim][partNum]+vel[dim][partNum]*dt+acc[dim][partNum]*dtSqrdByTwo, size[dim]);
            vel[dim][partNum] += acc[dim][partNum] * dtByTwo;
        }


        calculateForces();

        for (int dim=0; dim<DIMENSIONS; dim++) for (int partNum=0; partNum<partQty; partNum++)
        {
                vel[dim][partNum] += acc[dim][partNum] * dtByTwo;
                kineticEnergyTotal += getSquaredVectorLength(vel[dim]);
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

            //listSystemState();
            System.out.println((System.nanoTime()-nanoT)/stepsPerRun + "ns per step (DQ)");
        }
    }

    private void calculateForces()
    {
        double[] distance = new double[size.length];
        double[] forceVector = new double[size.length];
        double force = 0;
        for (int i=0; i<acc.length; i++) for (int j=0; j<acc[i].length; j++) acc[i][j] = 0.0;
        for (int p1=0; p1<partQty-1; p1++)
        {
            for (int p2=p1+1; p2<partQty; p2++)
            {
                for (int dim = 0; dim<distance.length; dim++)
                {
                    distance[dim] = periodicDistanceShift(pos[dim][p1]-pos[dim][p2], size[dim]);
                }
                double sqrdDistance = getSquaredDistance(distance);

                if (sqrdDistance<thresholdSqrd) {
                    force = PotentialsForces.forceLJ(sqrdDistance);
                    potentialEnergyAccumulator += PotentialsForces.potenLJ(sqrdDistance);

                    for (int dim = 0; dim < distance.length; dim++) {
                        forceVector[dim] = force * distance[dim];
                        acc[dim][p1] += forceVector[dim];
                        acc[dim][p2] -= forceVector[dim];
                        virialTotal += distance[dim] * forceVector[dim];
                    }
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
