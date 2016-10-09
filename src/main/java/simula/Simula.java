package simula;

import simula.ensembles.SparseEnsemble;

/**
 * Created by sash on 23.09.16.
 */
public class Simula {
//gesture

    public static void main(String[] args) {
        SparseEnsemble spe = new SparseEnsemble();
        if (spe.getParamsFromConsole()) {
            spe.start();
        }
    }
}
