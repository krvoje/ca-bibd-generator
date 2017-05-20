package org.krvoje

import scala.util.Random

package object bibd {

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
}
