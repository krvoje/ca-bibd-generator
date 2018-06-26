package org.krvoje.bibd
import scala.util.Random

case class Stochasticity(
  min: Int,
  max: Int,
  changeInterval: Either[Stochasticity, Int],
  initialValue: Option[Int] = None,
)(implicit referenceFrame: ReferenceFrame) {

  private var rate = initialValue.getOrElse {
    min + Random.nextInt((max - min).abs + 1)
  }

  private var up = true

  def value: () => BigInt = {

    val doWeChange = changeInterval match {
      case Left(stochasticity) => referenceFrame.unchanged > stochasticity.value()
      case Right(const) => referenceFrame.unchanged > const
    }

    if(doWeChange) {
      up = Random.nextBoolean()
      if(rate >= max) up = false
      if(rate <= min) up = true

      rate += (if(up) 1 else -1)
    }
    rate = math.min(rate, max)
    rate = math.max(rate, min)

    BigInt(rate)
  }
}