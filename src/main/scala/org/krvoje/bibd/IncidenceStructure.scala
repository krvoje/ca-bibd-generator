package org.krvoje.bibd

trait IncidenceStructure {

  val v: Int
  val r: Int
  val b: Int
  val k: Int
  val lambda: Int

  def heuristicDistance: Int

  def isBIBD: Boolean

  def flip(row: Int, Col: Int): Unit
  def active(row: Int, Col: Int): Boolean
  def dormant(row: Int, Col: Int): Boolean

  def incidences(row: Int, Col: Int): Int
  def setIncidence(row: Int, Col: Int, active: Boolean): Unit

  def rowIntersection(row: Int, Col: Int): Int
  def colIntersection(row: Int, Col: Int): Int

  def sumInRow(row: Int): Int
  def sumInCol(col: Int): Int
}
