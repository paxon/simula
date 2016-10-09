package simula.workers;

import java.util.concurrent.RecursiveTask;

/**
 * Created by HP on 10.10.2016.
 */
public class FJParticleVelocityUpdater extends RecursiveTask{
    private double[][] vel;
    private double[][] acc;
    private double[][] res;
    private double dtByTwo;
    private int partNum;



        public FJParticleVelocityUpdater(double[][] vel, double[][] acc, double[][] res, int partNum, double dtByTwo) {
            this.vel = vel;
            this.acc = acc;
            this.res = res;
            this.partNum = partNum;
            this.dtByTwo = dtByTwo;
        }

        protected double[][] compute() {
            double[][] res = vel.clone();
            for (int dim=0; dim<vel[partNum].length; dim++) res[partNum][dim] += acc[partNum][dim] * dtByTwo ;
            return res;
        }
}
