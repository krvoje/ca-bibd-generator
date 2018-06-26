name := "Cellular Automaton BIBD generator"
version := "1.0"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"

mainClass in assembly := Some("org.krvoje.bibd.Program")

licenses += "MIT License" -> url("https://opensource.org/licenses/MIT")

assemblyMergeStrategy in assembly := {
  case PathList("META-INF", xs @ _*) => MergeStrategy.discard
  case _ => MergeStrategy.first
}

test in assembly := {}