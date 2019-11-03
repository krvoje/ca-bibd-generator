package krvoje.bibd

import Implicits._

import scala.language.implicitConversions
import scala.util.Random

case class Stochasticity(
  min: Int,
  max: Int,
  increment: Option[Stochasticity] = None,
  changeThreshold: Option[Stochasticity] = None,
  initialValue: Option[Int] = None,
  isRandom: Boolean = false,
  isProbabilistic: Boolean = false,
)(implicit rf: ReferenceFrame) {

  import Stochasticity._

  private var rate = initialValue.getOrElse {
    randomValue
  }

  private var up = true

  def value: () => Int = if (isRandom) {
    randomValue
  } else{

    val isChange: Boolean = rf.unchanged >= changeThreshold

    if(isChange) {
      up = Random.nextBoolean()

      if (up) rate += increment
      else rate -= increment
    }
    rate = math.min(rate, max)
    rate = math.max(rate, min)


    if (isProbabilistic) {
      min + Random.nextInt(rate)
    } else {
      rate
    }
  }

  private def randomValue = min + Random.nextInt((max - min).abs + 1)
}

object Stochasticity {

  implicit def intToOptSto(int: Int)(implicit rf: ReferenceFrame): Option[Stochasticity] = Some(intToSto(int))

  implicit def intToSto(int: Int)(implicit rf: ReferenceFrame): Stochasticity = new Stochasticity(
    int, int, None, None, None
  )(rf)

  implicit def stochasticityToInt(sto: Stochasticity): Int = sto.value()

  implicit def optStochasticityToInt(optSto: Option[Stochasticity]): Int = optSto.map(_.value()).getOrElse(1)
}
