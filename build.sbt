enablePlugins(ScalaJSPlugin)

name := "Cellular Automaton BIBD generator"
version := "1.0"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"

licenses += "MIT License" -> url("https://opensource.org/licenses/MIT")

scalaJSUseMainModuleInitializer := true