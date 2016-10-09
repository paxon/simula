package simula.workers;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveTask;

/**
 * Created by HP on 09.10.2016.
 */
public class FJVelocityUpdater extends RecursiveTask<double[][]>{
    private double[][] vel;
    private double[][] acc;
    private double dtByTwo;


    public FJVelocityUpdater(double[][] vel, double[][] acc, double dtByTwo) {
        this.vel = vel;
        this.acc = acc;
        this.dtByTwo = dtByTwo;
    }

    protected double[][] compute() {
        double[][] res = vel.clone();
        List<ForkJoinTask> children = new ArrayList();
        for (int partNum=0; partNum<vel.length; partNum++) children.add( new FJParticleVelocityUpdater(vel,acc,res,partNum,dtByTwo));
        invokeAll(children);
        return res;
    }
}
