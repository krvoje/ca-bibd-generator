package org.krvoje.bibd
import scala.util.Random

case class Stochasticity(min: Int,
  max: Int,
  increment: () => Int,
  stochacity: Option[Stochasticity] = None,
  initialRate: Option[Int] = None)(implicit fieldLastChange: Moment) {
  def DEFAULT_CHANGE_INTERVAL: Int = 100

  private var rate = initialRate.getOrElse(min + (max - min) / 2)

  private var up = true

  def changeInterval: Int = stochacity.map(_.value()).getOrElse(DEFAULT_CHANGE_INTERVAL)

  def value: () => Int = {
    if((System.currentTimeMillis() - fieldLastChange.when > changeInterval)) {
      up = Random.nextBoolean()
      if(rate >= max) up = false
      if(rate <= min) up = true
      if(up) rate += increment()
      else rate -= increment()
    }
    rate
  }
}