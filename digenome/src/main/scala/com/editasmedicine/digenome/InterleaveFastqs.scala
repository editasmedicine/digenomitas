package com.editasmedicine.digenome

import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.fastq.{FastqRecord, FastqSource, FastqWriter}
import com.fulcrumgenomics.sopt.{arg, clp}

@clp(group=ClpGroups.Digenome, description=
  """
    |Takes in one or more pairs of fastq files and produces interleaved fastq. The same number of R1
    |and R2 files must be provided.  Pairing is auto-detected based on the name of the first read in
    |each file.
    |
    |If at any point the reads coming from a pair of fastq files do not have matching read names an
    |error will be raised.
  """)
class InterleaveFastqs
( @arg(flag='1', doc="One or more read 1 fastq files, optionally gzipped.") val readOne: Seq[PathToFastq],
  @arg(flag='2', doc="One or more read 2 fastq files, optionally gzipped.") val readTwo: Seq[PathToFastq],
  @arg(flag='o', doc="Output interleaved fastq file.") val output: PathToFastq = Io.StdOut
) extends EditasTool {

  validate(readOne.size == readTwo.size, "The same number of R1 and R2 files must be provided.")

  override def execute(): Unit = {
    val r1s = readOne.map(FastqSource(_).bufferBetter)
    val r2s = readTwo.map(FastqSource(_).bufferBetter)
    val out = FastqWriter(output)

    val iter = new FqPairIterator(r1s, r2s)
    iter.foreach { case (r1, r2) =>
      out += r1
      out += r2
    }

    out.close()
  }
}

/**
  * Iterator that takes in a collection of R1 and R2 fastqs and generates an iterator over pairs of fastqs.
  */
class FqPairIterator(r1s: Seq[BetterBufferedIterator[FastqRecord]], r2s: Seq[BetterBufferedIterator[FastqRecord]])
  extends Iterator[(FastqRecord,FastqRecord)] {
  require(r1s.size == r2s.size, "Must provide same number of R1 and R2 iterators.")

  // Match up the R1 and R2 iterators
  private val (r1Iterator, r2Iterator) = {
    val tmp = r2s.toBuffer
    val pairs = r1s.map { r1 =>
      val name = r1.head.name
      val r2 = tmp.indexWhere(_.head.name == name) match {
        case -1 => throw new IllegalArgumentException("Could not find an R2 iterator starting with read name $name")
        case i => tmp.remove(i)
      }

      (r1, r2)
    }

    require(tmp.isEmpty, s"Still had R2 iterators after matching all R2 iterators: ${tmp.map(_.head.name)}")
    (pairs.iterator.flatMap(_._1).bufferBetter, pairs.iterator.flatMap(_._2).bufferBetter)
  }

  /** Throws an exception if one iterator is exhausted and the other is not, or if both iterators have records
    * but their names differ.
    */
  def assertValid(): Unit = (r1Iterator.headOption, r2Iterator.headOption) match {
    case (None, None)     => Unit
    case (Some(r1), None) => throw new IllegalStateException(s"R2 exhausted, but R1 still has records. Next record: ${r1.name}")
    case (None, Some(r2)) => throw new IllegalStateException(s"R1 exhausted, but R2 still has records. Next record: ${r2.name}")
    case (Some(r1), Some(r2)) =>
      require(r1.name == r2.name, s"R1 and R2 out of sync. Next R1: ${r1.name}. Next R2: ${r2.name}.")
  }

  /** Returns true if either iterator has more records. */
  override def hasNext: Boolean = r1Iterator.hasNext || r2Iterator.hasNext

  /** Validates that the R1 and R2 iterators are sync'd then returns the next pair of records. */
  override def next(): (FastqRecord,FastqRecord) = {
    require(hasNext, "next() called on empty iterator.")
    assertValid()
    (r1Iterator.next(), r2Iterator.next())
  }
}
