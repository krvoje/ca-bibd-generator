package org.krvoje.bibd

class MatrixIncidenceStructure(val v: Int, val k: Int, val lambda: Int) extends IncidenceStructure{

  val (r, b) = computeParams(v, k, lambda)
  println(s"(v=$v, k=$k, λ=$lambda, b=$b, r=$r)")

  private val _incidences = Array.ofDim[Int](v, b)
  private val _rowIntersection = Array.ofDim[Int](v, b)
  private val _colIntersection = Array.ofDim[Int](v, b)
  private val _sumInRow = Array.ofDim[Int](v)
  private val _sumInCol = Array.ofDim[Int](b)
  private var _sumTotal = 0
  private val _sumIdeal = v*r
  private var _heuristicDistance = 0
  private var _maxHeuristicDistance = 0

  updateCache()

  override def rowIntersection(row: Int, col: Int): Int = _rowIntersection(row)(col)

  override def colIntersection(row: Int, col: Int): Int = _colIntersection(row)(col)

  override def sumInRow(row: Int): Int = _sumInRow(row)

  override def sumInCol(col: Int): Int = _sumInCol(col)

  override def dormant(row: Int, col: Int): Boolean = incidences(row, col) == 0

  override def active(row: Int, col: Int): Boolean = incidences(row, col) != 0

  override def setIncidence(row: Int, col: Int, active: Boolean): Unit = _incidences(row)(col) = if(active) 1 else 0

  override def incidences(row: Int, col: Int): Int = _incidences(row)(col)

  override def isBIBD: Boolean = {
    if(_sumTotal != _sumIdeal) {
      return false
    }

    forIndex(0, v) { row =>
      if(sumInRow(row) != r) return false
    }

    forIndex(0, b) { col =>
      if(sumInCol(col) != k) return false
    }

    forIndex(0, v) { row =>
      forIndex(row + 1, v) { otherRow =>
        if(_rowIntersection(row)(otherRow) != lambda) return false
      }
    }

    return true
  }

  override def flip(row: Int, col: Int): Unit = {
    val newVal: Int = math.abs(_incidences(row)(col) - 1)
    _incidences(row)(col) = newVal
    val increment: Int = math.pow(-1, newVal.toDouble + 1).intValue

    _sumTotal += increment
    _sumInRow(row) += increment
    _sumInCol(col) += increment
    _rowIntersection(row)(row) += increment
    _colIntersection(col)(col) += increment

    forIndex(0, v) { otherRow =>
      if (_incidences(otherRow)(col) == 1) {
        _rowIntersection(otherRow)(row) += increment
        _rowIntersection(row)(otherRow) += increment
      }
    }

    forIndex(0, b) {otherCol =>
      if (_incidences(row)(otherCol) == 1) {
        _colIntersection(otherCol)(col) += increment
        _colIntersection(col)(otherCol) += increment
      }
    }
    calculateHeuristicDistance
  }

  override def heuristicDistance: Int = _heuristicDistance


  private def computeParams(v: Int, k: Int, lambda: Int) = {
    val rDouble = lambda.toDouble * (v-1) / (k-1)
    val bDouble = rDouble * v / k

    if(!rDouble.isValidInt || !bDouble.isValidInt)
      throw new RuntimeException("Invalid BIBD parameters.")

    (rDouble.toInt, bDouble.toInt)
  }

  private def updateCache(): Unit = {
    _sumTotal = 0

    // Row sum
    forIndex(0, v) { row =>
      _sumInRow(row) = 0
      _sumInCol(row) = 0
      forIndex(0, b) { col =>
        _sumInRow(row) += _incidences(row)(col)
        _sumInCol(col) += _incidences(row)(col)
      }
      _sumTotal += sumInRow(row)
    }

    forIndex(0, v) { row1 =>
      forIndex(0, v) { row2 =>
        _rowIntersection(row1)(row2) = 0
        forIndex(0, b) { col =>
          _rowIntersection(row1)(row2) += _incidences(row1)(col) * _incidences(row2)(col)
        }
      }
    }

    calculateHeuristicDistance()
  }

  private def calculateHeuristicDistance(): Unit = {
    _heuristicDistance = 0
    forIndex(0, v) { row =>
      _heuristicDistance += math.abs(sumInRow(row) - r)
    }
    forIndex(0, b) { col =>
      _heuristicDistance += math.abs(sumInCol(col) - k)
    }
    forIndex(0, v) { row1 =>
      forIndex(0, v) { row2 =>
        if(row1 != row2) {
          _heuristicDistance += _rowIntersection(row1)(row2) - lambda
        }
      }
    }
    _maxHeuristicDistance = b * (v - k) + (v * v - v) * (r - lambda)
  }

  override def toString: String = {
    val sb = StringBuilder.newBuilder
    forIndex(0, v) {
      row => sb.append("\n")
        forIndex(0, b) {
          col => sb.append(_incidences(row)(col))
        }
    }
    sb.toString()
  }
}
