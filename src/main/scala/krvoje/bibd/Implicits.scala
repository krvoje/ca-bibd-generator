package krvoje.bibd

import scala.language.implicitConversions

object Implicits {

  def now(): Long = System.currentTimeMillis()

  def forIndex(from: Int, to: Int)(fn: Int => Unit): Unit = {
    var index = from
    while(index < to) {
      fn(index)
      index += 1
    }
  }

  def min(left: Int, right: Int): Int = if(left < right) left else right
  def max(left: Int, right: Int): Int = if(left > right) left else right

  implicit def anyOpt[A, B](a: A)(implicit f: A => B): Option[B] = Some(f(a))
  implicit def anyMethod[T](e: T): () => T = () => e
  implicit def doubleToInt(number: Double): Int = number.intValue()
  //implicit def leftLift[X,Y](value : X): Either[X,Y] = Left(value)
  //implicit def rightLift[X,Y](value : Y): Either[X,Y] = Right(value)
}
