package org.krvoje.bibd
import scala.util.Random

case class Stochasticity(
  min: Int,
  max: Int,
  changeInterval: Option[Stochasticity] = None,
  initialValue: Option[Int] = None)(implicit referenceFrame: ReferenceFrame) {

  private val DefaultSto = BigInt(1)

  private var rate = initialValue.getOrElse {
    min + Random.nextInt((max - min).abs + 1)
  }

  private var up = true

  def value: () => BigInt = {
    if(referenceFrame.iterationsSinceLastChange > changeInterval.map(_.value()).getOrElse(DefaultSto)) {
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