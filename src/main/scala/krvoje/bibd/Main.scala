package krvoje.bibd

import Implicits._

object Main {
  def main(args: Array[String]): Unit = {
    if(args.length != 3) {
      System.out.println("Usage: ")
      System.out.println("java -jar ca-bibd.jar v k lambda")
      System.exit(0)
    }

    val vertices = Integer.parseInt(args(0))
    val blocksPerVertex = Integer.parseInt(args(1))
    val lambda = Integer.parseInt(args(2))

    implicit val lastChange: ReferenceFrame = ReferenceFrame(now())

    val randomCABIBD = new MutationCABIBD(vertices, blocksPerVertex, lambda)
    randomCABIBD.findBIBD
  }
}
