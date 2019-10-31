package krvoje.bibd

case class ReferenceFrame(
  var lastChangeTimestamp: Long,
  var currentIteration: Int = 0,
  var lastChangeIteration: Int = 0,
  var generation: Int = 0,
  var staleResets: Int = 0
) {

  def iterate(): Unit = currentIteration += 1
  def changed(staleReset: Boolean): Unit = {
    if(staleReset) staleResets += 1
    lastChangeTimestamp = System.currentTimeMillis()
    lastChangeIteration = currentIteration
    generation += 1
  }

  def unchanged: Int = (currentIteration - lastChangeIteration)
  def lastChangeHowLongAgo: Long = System.currentTimeMillis() - lastChangeTimestamp
}
