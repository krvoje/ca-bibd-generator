package org.krvoje.bibd

case class ReferenceFrame(
  var lastChangeTimestamp: Long,
  var currentIteration: BigInt = 0,
  var lastChangeIteration: BigInt = 0,
  var generation: BigInt = 0
) {

  def unchanged = currentIteration - lastChangeIteration
  def iterate = currentIteration += 1
  def changed = {
    lastChangeTimestamp = System.currentTimeMillis()
    lastChangeIteration = currentIteration
    generation += 1
  }

  def iterationsSinceLastChange = (currentIteration - lastChangeIteration)
  def lastChangeHowLongAgo = System.currentTimeMillis() - lastChangeTimestamp
}
