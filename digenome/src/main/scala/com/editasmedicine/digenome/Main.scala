package com.editasmedicine.digenome

import com.editasmedicine.commons.clp.ClpMain

/** Main object that is run from the command line. */
object Main {
  /** The main method */
  def main(args: Array[String]): Unit = new ClpMain("grna-qc").makeItSoAndExit(args)
}
