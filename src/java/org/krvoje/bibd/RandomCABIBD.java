package org.krvoje.bibd;

import java.util.Random;

public class RandomCABIBD {

    private static final Random rnd = new java.util.Random(System.currentTimeMillis());

    public static void main(String[] args) throws Exception {
        int v,k,lambda;
        if(args.length != 3) {
            System.out.println("Usage: ");
            System.out.println("java -jar ca-bibd.jar v k lambda");
        }
        else {
            v = Integer.parseInt(args[0]);
            k = Integer.parseInt(args[1]);
            lambda = Integer.parseInt(args[2]);

            IncidenceStructure bibd = new IncidenceStructure(v,k,lambda);
            new RandomCABIBD().randomCaBibd(bibd);

            System.out.println("An incidence matrix for the given parameters found!");
            System.out.print(bibd);
        }
    }

    private void randomCaBibd(IncidenceStructure is){
        int maxChangeFactor = (is.v-1)*Math.abs(is.r-is.lambda) + is.b;
        for(int i=0; i<is.v; i++) {
            for (int j=0; j<is.r; j++) {
                is.incidences[i][j]=1;
            }
        }
        is.updateCache();


        while(true)
        {
            is.generations++;
            if(is.isBIBD()) return;

            int row = rnd.nextInt(is.v);
            int activeCol = randomActiveIn(is, row);
            int dormantCol = randomDormantIn(is, row);

            int cfa = changeFactor(is,row,activeCol);
            int cfd = changeFactor(is,row,dormantCol);

            int doWeChange = rnd.nextInt(maxChangeFactor);
            if(doWeChange < cfa && doWeChange < cfd) {
                is.flip(row,activeCol);
                is.flip(row,dormantCol);
            }
        }
    }

    private static int changeFactor(IncidenceStructure is, int row, int col) {
        int changeFactor = 0;
        for(int otherRow=0; otherRow<is.v; otherRow++){
            if(otherRow==row) continue;
            int delta = is.lambda - is.rowIntersection[row][otherRow];
            if(delta < 0) {
                if(is.active(row,col) && is.active(otherRow, col))
                    changeFactor -= delta;
                else if(is.active(row,col) && is.dormant(otherRow, col))
                    changeFactor += delta;
            }
            if(delta < 0) {
                if(is.dormant(row, col) && is.active(otherRow, col))
                    changeFactor+=delta;
                else if(is.dormant(row, col) && is.dormant(otherRow, col))
                    changeFactor-=delta;
            }
        }

        int delta = is.k - is.sumInCol[col];
        if(is.active(row,col))
            changeFactor -= delta;
        else if(is.dormant(row,col))
            changeFactor += delta;

        return changeFactor >= 0 ? changeFactor : 0;
    }

    public int randomActiveIn(IncidenceStructure is, int row) {
        int col = rnd.nextInt(is.b);
        while(!is.active(row, col))
            col = rnd.nextInt(is.b);
        return col;
    }

    public int randomDormantIn(IncidenceStructure is, int row) {
        int col = rnd.nextInt(is.b);
        while(!is.dormant(row, col))
            col = rnd.nextInt(is.b);
        return col;
    }
}
