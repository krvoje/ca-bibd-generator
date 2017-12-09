package org.krvoje.bibd
import scala.util.Random

case class Stochasticity(min: Int,
  max: Int,
  increment: () => Int = () => 1,
  stochasticity: Option[Stochasticity] = None,
  initialRate: Option[Int] = None)(implicit fieldLastChange: LastChange) {
  def DEFAULT_CHANGE_INTERVAL: Int = 5

  private var rate = initialRate.getOrElse {
    min + Random.nextInt((max - min).abs)
  }

  private var up = true

  def changeInterval: Int = stochasticity.map(_.value()).getOrElse(DEFAULT_CHANGE_INTERVAL)

  def value: () => Int = {
    if((System.currentTimeMillis() - fieldLastChange.timestamp > changeInterval)) {
      up = Random.nextBoolean()
      if(rate >= max) up = false
      if(rate <= min) up = true
      if(up) rate += increment()
      else rate -= increment()
    }
    if(rate > max) max
    else if(rate < min) min
    else rate
  }
}