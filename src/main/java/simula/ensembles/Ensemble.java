package simula.ensembles;

/**
 * Created by HP on 10.10.2016.
 */
public interface Ensemble {

    boolean init (double[] size, int partQty);
    void initParticles(double initKinetEnergy, double minVelocity, double maxVelocity);
    double getMeanTemperature();
    double getMeanEnergy();
    double getMeanPressure();

}
