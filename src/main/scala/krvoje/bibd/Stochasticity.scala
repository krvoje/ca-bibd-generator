package krvoje.bibd

import krvoje.bibd.util.Implicits._

import scala.util.Random

case class Stochasticity(
  min: Int,
  max: Int,
  changeTreshold: Option[Stochasticity] = None,
  increment: Option[Stochasticity] = None,
  initialValue: Option[Int] = None,
)(implicit rf: ReferenceFrame) {

  private var rate = initialValue.getOrElse {
    min + Random.nextInt((max - min).abs + 1)
  }

  private var up = true

  def value: () => Int = {

    val doWeChange = {
      val treshold: Int = changeTreshold
      rf.unchanged > treshold
    }

    if(doWeChange) {
      up = Random.nextBoolean()
      if(rate >= max) up = false
      if(rate <= min) up = true

      if (up) rate += increment
      else rate -= increment
    }
    rate = math.min(rate, max)
    rate = math.max(rate, min)

    rate
  }
}

object Stochasticity {

  implicit def intToOptSto(int: Int)(implicit rf: ReferenceFrame): Option[Stochasticity] = Some(intToSto(int))

  implicit def intToSto(int: Int)(implicit rf: ReferenceFrame): Stochasticity = new Stochasticity(
    int, int, None, None, None
  )(rf)

  implicit def optStoToInt(optSto: Option[Stochasticity]): Int = optSto.map(_.value()).getOrElse(1)
}