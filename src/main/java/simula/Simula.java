package simula;

/**
 * Created by sash on 23.09.16.
 */
public class Simula {
//gesture

    public static void main(String[] args) {
        //Boxie b = new Boxie();
        //b.initSize();
        //b.initPar(BigInteger.valueOf(100));
        //b.listPar();
        //Sy4stem.out.println(b.getPotEnergy());

        SparseEnsemble spe = new SparseEnsemble();
        System.out.println(PotentialsForces.forceLJ(1));
        System.out.println(PotentialsForces.potenLJ(1));
        if (spe.getParamsFromConsole()) {
            //spe.listParticlesState();
            spe.start();
        }
    }


}
