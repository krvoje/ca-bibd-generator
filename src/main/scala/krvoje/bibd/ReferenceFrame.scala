package krvoje.bibd

case class ReferenceFrame(
  var lastChangeTimestamp: Long,
  var currentIteration: BigInt = 0,
  var lastChangeIteration: BigInt = 0,
  var generation: BigInt = 0,
  var staleResets: BigInt = 0
) {

  def iterate = currentIteration += 1
  def changed(staleReset: Boolean) = {
    if(staleReset) staleResets += 1
    lastChangeTimestamp = System.currentTimeMillis()
    lastChangeIteration = currentIteration
    generation += 1
  }

  def unchanged = (currentIteration - lastChangeIteration)
  def lastChangeHowLongAgo = System.currentTimeMillis() - lastChangeTimestamp
}
