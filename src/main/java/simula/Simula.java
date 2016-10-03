package simula;

/**
 * Created by sash on 23.09.16.
 */
public class Simula {

    public static void main(String[] args) {
        //Boxie b = new Boxie();
        //b.initSize();
        //b.initPar(BigInteger.valueOf(100));
        //b.listPar();
        //System.out.println(b.getPotEnergy());

        SparseEnsemble spe = new SparseEnsemble();
        spe.getParamsFromConsole();
        spe.start();



    }


}
