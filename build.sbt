name := "Cellular Automaton BIBD generator"
version := "1.0"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"

mainClass in assembly := Some("org.krvoje.bibd.RandomCABIBD")

assemblyMergeStrategy in assembly := {
  case PathList("META-INF", xs @ _*) => MergeStrategy.discard
  case _ => MergeStrategy.first
}