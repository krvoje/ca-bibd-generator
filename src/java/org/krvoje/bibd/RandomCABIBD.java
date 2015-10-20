package org.krvoje.bibd;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;

import java.util.Random;

public class RandomCABIBD {

    private static final Random rnd = new java.util.Random(System.currentTimeMillis());

    private final int[][] changeFactor;
    private final IncidenceStructure is;
    private int currentMaxChangeFactor=0;

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
        randomize();
    }

    private void randomize() {
        for(int i = 0; i < is.k(); i++) {
            for (int col =0; col < is.b(); col++) {
                int row = rnd.nextInt(is.v());
                while(is.active(row,col))
                    row = rnd.nextInt(is.v());
                is.setIncidence(row, col, true);
            }
        }
        calculateChangeFactors();
    }

    private IncidenceStructure getBibd(){
        int generations = 1;
        int iterations = 0;
        int unchanged = 0;

        int col, activeRow, dormantRow;
        int cfa, cfd, randomChangeFactor;
        boolean doWeChange, stale;

        col = -1;
        while(true)
        {
            col = (++col) % is.b();
            activeRow = randomActiveInCol(is, col);
            dormantRow = randomDormantInCol(is, col);

            if(is.isBIBD()) {
                System.out.println("Generations: " + generations);
                System.out.println("Iterations: " + iterations);
                return this.is; // Sparta
            }

            iterations++;

            cfa = changeFactor[activeRow][col];
            cfd = changeFactor[dormantRow][col];

            randomChangeFactor = currentMaxChangeFactor > 1 ?
                    rnd.nextInt(currentMaxChangeFactor)
                    : 1;

            if(unchanged > is.v()*is.b()*10) {
                System.out.println("Stale, randomizing...");
                this.randomize();
                unchanged = 0;
                continue;
            }
            doWeChange = (cfa > randomChangeFactor && cfd > randomChangeFactor)
                    || currentMaxChangeFactor <= 0
                    || unchanged > is.v()*is.b()*3 && (cfa>0 || cfd>0)
                    ;

            if (doWeChange) {
                is.flip(activeRow, col);
                is.flip(dormantRow, col);
                calculateChangeFactors();
                //System.out.println(unchanged + ", " + currentMaxChangeFactor);
                unchanged = 0;
                generations++;
            }
            else {
                unchanged ++;
            }
        }
    }

    private void calculateChangeFactors() {
        currentMaxChangeFactor = 0;
        for(int row = 0; row < is.v(); row++) {
            for (int col =0; col < is.b(); col++) {
                changeFactor[row][col] = changeFactor(row,col);
                currentMaxChangeFactor = max(changeFactor[row][col], currentMaxChangeFactor);
            }
        }
    }

    private int changeFactor(int row, int col) {
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
                    changeFactor += delta;
                else if(is.dormant(row, col) && is.dormant(otherRow, col))
                    changeFactor -= delta;
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

    public int randomActiveInCol(IncidenceStructure is, int col) {
        int row = rnd.nextInt(is.v());
        while(!is.active(row, col))
            row = rnd.nextInt(is.v());
        return row;
    }

    public int randomDormantInCol(IncidenceStructure is, int col) {
        int row = rnd.nextInt(is.v());
        while(!is.dormant(row, col))
            row = rnd.nextInt(is.v());
        return row;
    }
}
