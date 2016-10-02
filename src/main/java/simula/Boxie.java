package simula;

import java.math.BigInteger;
import java.util.*;

/**
 * Created by sash on 24.09.16.
 */
public class Boxie {
    private Map<String,Double> size = new HashMap<String,Double>();
    private BigInteger partQty;
    private List<Partie> parties = new ArrayList<>();

    private double[] cellSize = new double[3];
    


    private void setSize(Map<String, Double> size) {
        this.size = size;
    }

    public Map<String, Double> getSize() {
        return size;
    }

    public void initSize(){
        Scanner sc = new Scanner(System.in);
        int numOrd = 0;
        System.out.println("Enter number of coords: ");
        numOrd=sc.nextInt();
        sc.nextLine();
        for (int i=0; i<numOrd; i++){
            System.out.println("Enter bound for x"+(int) (i+1)+": ");
            this.size.put("x"+(int) (i+1),sc.nextDouble());
            sc.nextLine();
        }
    }

    public void initPar(BigInteger qty){
        this.partQty = qty;
        for (BigInteger bigI=BigInteger.valueOf(0); bigI.compareTo(partQty)<0; bigI=bigI.add(BigInteger.ONE))
        {
            Partie p = new Partie();
            p.init(size);
            parties.add(p);
        }

    }

    public void listPar(){
        System.out.println("list");
        for (Partie p : this.parties
                ) {
            p.listCoo();
        }
    }

    public double getPotEnergy() {
        double res = 0;
        for (Partie p1 : this.parties ) {
            for (Partie p2 : this.parties ) {
                if (p1!=p2) res += actingBinaryPotentials(p1,p2);
            }
        }
        return res;
    }

    private double actingBinaryPotentials (Partie p1, Partie p2) {
        return PotentialsForces.potenLJ(getSqrdDistance(p1,p2));
    }

    private double getSqrdDistance(Partie p1, Partie p2) {
        double res = 0;
        for (Map.Entry<String,Double> ord :
                size.entrySet()) {
            res += p1.getOrd(ord.getKey())-p2.getOrd(ord.getKey());
            res *= res;
        }
        return res;
    }

}

