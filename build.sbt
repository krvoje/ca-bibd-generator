name := "Cellular Automaton BIBD generator"
version := "1.0"

libraryDependencies += "org.scalactic" %% "scalactic" % "3.0.1"
libraryDependencies +="com.lihaoyi" %% "utest" % "0.4.5" % "test"
testFrameworks += new TestFramework("utest.runner.Framework")

mainClass in assembly := Some("org.krvoje.bibd.RandomCABIBD")

assemblyMergeStrategy in assembly := {
  case PathList("META-INF", xs @ _*) => MergeStrategy.discard
  case _ => MergeStrategy.first
}