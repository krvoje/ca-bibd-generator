package org.krvoje.bibd;

/**
 * Created by krvoje on 10/8/15.
 */
public interface IncidenceStructure {

    public boolean isBIBD();
    public void flip(int row, int col);
    public boolean active(int row, int col);
    public boolean dormant(int row, int col);

    public int incidences(int row, int col);
    public void setIncidence(int row, int col, boolean active);

    public int rowIntersection(int row, int otherRow);
    public int colIntersection(int col, int otherCol);

    public int sumInRow(int row);
    public int sumInCol(int col);
}
