package simula;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Created by sash on 24.09.16.
 */

public class Partie {
    private Map<String,Double> coo = new HashMap<String, Double>();

    public void init(Map<String,Double> coo){
        for (Map.Entry<String,Double> ord:coo.entrySet())
        {
            this.coo.put(ord.getKey(), ThreadLocalRandom.current().nextDouble(0,ord.getValue()));
        }
    }

    public void listCoo() {
        System.out.print("Particle "+this.toString()+" ( ");
        for (Map.Entry<String, Double> ord: coo.entrySet())
        {
            System.out.print(ord.getKey()+": "+ord.getValue()+ " ");
        }
        System.out.println(")");
    }

    public double getOrd(String ordStr) {
        return coo.get(ordStr);
    }
}
