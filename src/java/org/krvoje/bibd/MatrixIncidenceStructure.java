package org.krvoje.bibd;

public class MatrixIncidenceStructure implements IncidenceStructure {

	private int[][] incidences;
	private int[][] rowIntersection;
	private int[][] colIntersection;
	private int[] sumInRow;
	private int[] sumInCol;
	
	private final int v, r, b, k, lambda;
	private int sumTotal;
	private int sumIdeal;
	
	private int heuristicDistance;
	private int maxHeuristicDistance;
	
	public MatrixIncidenceStructure(int v, int k, int lambda) throws Exception
	{
		this.v = v;
		this.k = k;
		this.lambda = lambda;
		
		Double r_ = 1.0*lambda*(v-1)/(k-1);
		Double b_ = r_*v/k;
		
		if(r_.doubleValue() != r_.intValue()
				|| b_.doubleValue() != b_.intValue())
		{
			throw new Exception("Invalid BIBD parameters.");
		}
		
		this.r = r_.intValue();
		this.b = b_.intValue();
		
		this.incidences = new int[v][b];
		this.rowIntersection = new int[v][v];
		this.colIntersection = new int[b][b];
		this.sumInRow = new int[v];
		this.sumInCol = new int[b];
		this.sumTotal = 0;
		this.sumIdeal = v*r;

		for(int row=0; row<v; row++) {
			for(int col =0; col <b; col++) {
				this.incidences[row][col]=0;
			}
		}

		this.updateCache();
		System.out.println(String.format("(v=%d, k=%d, λ=%d, b=%d, r=%d)", v, k, lambda, b, r));
	}

	@Override
	public int v() {
		return this.v;
	}

	@Override
	public int r() {
		return this.r;
	}

	@Override
	public int b() {
		return this.b;
	}

	@Override
	public int k() {
		return this.k;
	}

	@Override
	public int vertices() {
		return this.v;
	}

	@Override
	public int verticesInBlock() {
		return this.r;
	}

	@Override
	public int blocks() {
		return this.b;
	}

	@Override
	public int blocksPerVertex() {
		return this.k;
	}

	@Override
	public int lambda() {
		return this.lambda;
	}

	public boolean isBIBD() {
		if(sumTotal != sumIdeal) {
			return false;
		}

		for(int row =0; row<v; row++)
		{
			if(sumInRow[row] != r) {
				return false;
			}
		}

		for(int col=0; col<b; col++)
		{
			if(sumInCol[col] != k) {
				return false;
			}
		}

		for(int row=0; row<v-1; row++)
		for(int otherRow=row+1; otherRow<v; otherRow++)
		{
			if(rowIntersection[row][otherRow] != lambda)
				return false;
		}

		return true;
	}

	public boolean active(int row, int col)
	{
		return this.incidences[row][col] != 0;
	}

	public boolean dormant(int row, int col)
	{
		return this.incidences[row][col] == 0;
	}

	@Override
	public int incidences(int row, int col) {
		return this.incidences[row][col];
	}

	@Override
	public void setIncidence(int row, int col, boolean active) {
		//this.incidences[row][col] = active ? 1 : 0;
		if(!(this.active(row,col) && active))
			this.flip(row,col);
	}

	@Override
	public int rowIntersection(int row, int otherRow) {
		return this.rowIntersection[row][otherRow];
	}

	@Override
	public int colIntersection(int col, int otherCol) {
		return this.colIntersection[col][otherCol];
	}

	@Override
	public int sumInRow(int row) {
		return this.sumInRow[row];
	}

	@Override
	public int sumInCol(int col) {
		return this.sumInCol[col];
	}

	public void flip(int row, int col) {
		int newVal = Math.abs(incidences[row][col] - 1);
		incidences[row][col]=newVal;
		int increment = new Double(Math.pow(-1, newVal + 1)).intValue();

		this.sumTotal += increment;
		this.sumInRow[row] += increment;
		this.sumInCol[col] += increment;
		this.rowIntersection[row][row] += increment;
		this.colIntersection[col][col] += increment;

		for(int otherRow=0; otherRow<this.v; otherRow++)
		{
			if(this.incidences[otherRow][col]==1)
			{
				rowIntersection[otherRow][row] += increment;
				rowIntersection[row][otherRow] += increment;
			}
		}

		for(int otherCol=0; otherCol <this.b; otherCol++)
		{
			if(this.incidences[row][otherCol]==1)
			{
				colIntersection[otherCol][col] += increment;
				colIntersection[col][otherCol] += increment;
			}
		}
		this.calculateHeuristicDistance();
	}

	private void calculateHeuristicDistance()
	{
		this.heuristicDistance=0;
		for(int col =0; col <this.b; col++)
		{
			heuristicDistance+= Math.abs(sumInCol[col]-this.k);
		}
		for(int row1=0; row1<this.v; row1++)
			for(int row2=0; row2<this.v; row2++)
			{
				if(row1==row2) continue;
				heuristicDistance+=Math.abs(rowIntersection[row1][row2]-this.lambda);
			}
		this.maxHeuristicDistance=b*(v-k) + (v*v - v)*(r-lambda);
	}
	
	public void updateCache() {

		this.sumTotal = 0;
		for(int row=0; row<this.v; row++) {
			this.sumInRow[row] = 0;
			for(int col=0; col<this.b; col++) {
				this.sumInRow[row] += incidences[row][col];
			}
			this.sumTotal += this.sumInRow[row];
		}
		for(int col=0; col<this.b; col++) {
			this.sumInCol[col] = 0;
			for(int row=0; row<this.v; row++) {
				this.sumInCol[col] += incidences[row][col];
			}
		}
		
		for(int row1=0; row1<this.v; row1++) {
			for(int row2=0; row2<this.v; row2++) {
				this.rowIntersection[row1][row2] = 0;
				for(int col=0; col<this.b; col++) {
					rowIntersection[row1][row2] += this.incidences[row1][col]*this.incidences[row2][col];
				}
			}
		}

		for(int col1=0; col1<this.b; col1++) {
			for(int col2=0; col2<this.b; col2++) {
				this.colIntersection[col1][col2] = 0;
				for(int row =0; row<this.v; row++) {
					colIntersection[col1][col2] += incidences[row][col1]*incidences[row][col2];
				}
			}
		}

		this.calculateHeuristicDistance();
	}

	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<v; i++) {
			for(int j=0; j <b; j++) {
				sb.append(incidences[i][j]);
			}
			sb.append("\n");
		}
		return sb.toString();
	}
}
