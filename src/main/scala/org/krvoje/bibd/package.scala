package org.krvoje

package object bibd {
  def forIndex(from: Int, to: Int)(fn: Int => Unit) {
    var index = from
    while(index < to) {
      fn(index)
      index += 1
    }
  }
}
