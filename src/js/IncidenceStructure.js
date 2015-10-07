function IncidenceStructure(v,k,lambda) {
    var self = this;

    self.v = v;
    self.k = k;
    self.lambda = lambda;

    var r_ = 1.0*lambda*(v-1)/(k-1);
    var b_ = r_*v/k;

    if(r != parseInt(r_, 10)
                || b_ != parseInt(b_, 10))
    {
        throw "Invalid BIBD parameters.";
    }

    self.r = parseInt(r_, 10);
    self.b = parseInt(b_, 10);

    self.incidences = [[]];
    self.rowIntersection = [[]];
    self.colIntersection = [[]];
    self.sumInRow = [];
    self.sumInCol = [];
    self.sumTotal = 0;
    self.generations = 0;

    self.updateCache();
    console.log(String.format("(v=%d, k=%d, λ=%d, b=%d, r=%d)", v, k, lambda, b, r));
}

IncidenceStructure.prototype.isBIBD = function() {
            var self = this;
            if(sumTotal != sum_ideal) {
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

    		console.log("Generations: " + self.generations);

    		return true;
}

IncidenceStructure.prototype.active = function(row,col) {
    var self = this;
    return self.incidences[row][col] != 0;
}

IncidenceStructure.prototype.dormant = function(row,col) {
    var self = this;
    return self.incidences[row][col] == 0;
}

IncidenceStructure.prototype.flip = function(row, col) {
    var self = this;

	var int newVal = abs(incidences[row][col] - 1);
	incidences[row][col]=newVal;
	var increment = pow(-1, newVal + 1);

	self.sumTotal += increment;
	self.sumInRow[row] += increment;
	self.sumInCol[col] += increment;
	self.rowIntersection[row][row] += increment;
	self.colIntersection[col][col] += increment;

	for(otherRow=0; otherRow<self.v; otherRow++)
	{
		if(self.incidences[otherRow][col]==1)
		{
			rowIntersection[otherRow][row] += increment;
			rowIntersection[row][otherRow] += increment;
		}
	}

	for(int otherCol=0; otherCol <self.b; otherCol++)
	{
		if(self.incidences[row][otherCol]==1)
		{
			colIntersection[otherCol][col] += increment;
			colIntersection[col][otherCol] += increment;
		}
	}
	//self.calculateHeuristicDistance();
}