package simula;

import simula.ensembles.SparseEnsembleDQ;
import simula.ensembles.SparseEnsembleQD;

/**
 * Created by sash on 23.09.16.
 */
public class Simula {
//gesture

    public static void main(String[] args) {
        SparseEnsembleDQ spe = new SparseEnsembleDQ();
        if (spe.getParamsFromConsole()) {
            spe.start();
        }
    }
}
