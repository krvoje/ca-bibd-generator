package org.krvoje

package object bibd {

  def now = System.currentTimeMillis()

  def forIndex(from: Int, to: Int)(fn: Int => Unit) = {
    var index = from
    while(index < to) {
      fn(index)
      index += 1
    }
  }

  def min(left: BigInt, right: BigInt): BigInt = if(left < right) left else right
  def max(left: BigInt, right: BigInt): BigInt = if(left > right) left else right

  implicit def anyOpt[T](e: T): Option[T] = Some(e)
  implicit def anyMethod[T](e: T): () => T = () => e
  implicit def anySto(e: Int)(implicit lastChange: ReferenceFrame): Stochasticity = Stochasticity(e, e+1)
  implicit def anyStoOpt(e: Int)(implicit lastChange: ReferenceFrame): Option[Stochasticity] = Stochasticity(e, e+1)
  implicit def sto2bigint(sto: Stochasticity): BigInt = sto.value()
}
