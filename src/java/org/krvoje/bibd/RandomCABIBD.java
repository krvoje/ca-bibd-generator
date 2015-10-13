package org.krvoje.bibd;
import static java.lang.Math.abs;
import static java.lang.Math.max;

import java.util.Random;

public class RandomCABIBD {

    private static final Random rnd = new java.util.Random(System.currentTimeMillis());

    private final int[][] changeFactor;
    private final IncidenceStructure is;

    public static void main(String[] args) throws Exception {
        int vertices, blocksPerVertex,lambda;
        if(args.length != 3) {
            System.out.println("Usage: ");
            System.out.println("java -jar ca-bibd.jar v k lambda");
        }
        else {
            vertices = Integer.parseInt(args[0]);
            blocksPerVertex = Integer.parseInt(args[1]);
            lambda = Integer.parseInt(args[2]);

            IncidenceStructure is = new RandomCABIBD(vertices, blocksPerVertex, lambda).getBibd();

            System.out.println("An incidence matrix for the given parameters found!");
            System.out.print(is);
        }
    }

    public RandomCABIBD(int v, int k, int lambda) throws Exception
    {
        this.is = new MatrixIncidenceStructure(v, k, lambda);
        this.changeFactor = new int[is.v()][is.b()];
        for(int i = 0; i < is.k(); i++) {
            for (int col =0; col < is.b(); col++) {
                int row = rnd.nextInt(is.v());
                while(is.active(row,col))
                    row = rnd.nextInt(is.v());
                is.setIncidence(row, col, true);
            }
        }

        for(int row=0; row<is.v(); row++)
            for(int col =0; col <is.b(); col++)
                changeFactor[row][col] = changeFactor(is,row,col);
    }

    private IncidenceStructure getBibd(){
        int maxChangeFactor = (is.v()-1)*abs(is.r() - is.lambda()) + (is.b() - is.k()) + (is.v() - is.r());
        int generations = 0;
        int iterations = 0;
        int unchanged = 0;
        int maxUnchanged = 0;

        while(true)
        {
            iterations++;

            if(unchanged > 10) {
                //System.out.println(unchanged + ", " + maxUnchanged);
                //System.out.println(is);
            }

            if(is.isBIBD()) {
                System.out.println("Generations: " + generations);
                System.out.println("Iterations: " + iterations);
                System.out.println("Max unchanged number of iterations: " + maxUnchanged);
                return this.is; // Sparta
            }

            int col = rnd.nextInt(is.b());
            int activeRow = randomActiveInCol(is, col);
            int dormantRow = randomDormantInCol(is, col);

            int cfa = changeFactor(is, activeRow, col);
            int cfd = changeFactor(is, dormantRow, col);

            int doWeChange = rnd.nextInt(maxChangeFactor);
            if(doWeChange < cfa && doWeChange < cfd) {
                is.flip(activeRow, col);
                is.flip(dormantRow, col);
                maxUnchanged=max(unchanged,maxUnchanged);
                unchanged = 0;
                generations++;
            }
            else {
                //System.out.println(doWeChange + ", " + cfa+ ", " + cfd);
                unchanged ++;
            }
        }
    }

    private static int changeFactor(IncidenceStructure is, int row, int col) {
        int changeFactor = 0;
        int delta = 0;

        // Max increment of: (v-1) * r - lambda
        for(int otherRow=0; otherRow<is.v(); otherRow++){
            if(otherRow==row) continue;
            delta = is.lambda() - is.rowIntersection(row, otherRow);
            if(delta < 0) {
                if(is.active(row,col) && is.active(otherRow, col))
                    changeFactor -= delta;
                else if(is.active(row,col) && is.dormant(otherRow, col))
                    changeFactor += delta;
            }
            if(delta > 0) {
                if(is.dormant(row, col) && is.active(otherRow, col))
                    changeFactor+=delta;
                else if(is.dormant(row, col) && is.dormant(otherRow, col))
                    changeFactor-=delta;
            }
        }

        // Max increment of: (b - k)
        delta = is.k() - is.sumInCol(col);
        if(is.active(row,col))
            changeFactor -= delta;
        else if(is.dormant(row,col))
            changeFactor += delta;

        // Max increment of: (v - r)
        delta = is.r() - is.sumInRow(row);
        if(is.active(row,col))
            changeFactor -= delta;
        else if(is.dormant(row,col))
            changeFactor += delta;

        return changeFactor;
    }

    public int randomActiveInRow(IncidenceStructure is, int row) {
        int col = rnd.nextInt(is.b());
        while(!is.active(row, col))
            col = rnd.nextInt(is.b());
        return col;
    }

    public int randomActiveInCol(IncidenceStructure is, int col) {
        int row = rnd.nextInt(is.v());
        while(!is.active(row, col))
            row = rnd.nextInt(is.v());
        return row;
    }

    public int randomDormantInRow(IncidenceStructure is, int row) {
        int col = rnd.nextInt(is.b());
        while(!is.dormant(row, col))
            col = rnd.nextInt(is.b());
        return col;
    }

    public int randomDormantInCol(IncidenceStructure is, int col) {
        int row = rnd.nextInt(is.v());
        while(!is.dormant(row, col))
            row = rnd.nextInt(is.v());
        return row;
    }
}
