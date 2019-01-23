package com.editasmedicine.aligner

import java.io.PrintStream

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.LinearMatrix
import com.fulcrumgenomics.alignment.Mode.Glocal
import htsjdk.samtools.reference.ReferenceSequenceFile
import htsjdk.samtools.util.SequenceUtil

import scala.math.{abs, max, min}

/** Represents the alignment of a guide to a section of a reference sequence / chromosome. The alignment
  * is always returned in the same orientation as the guide.  I.e. if the guide matched the negative
  * strand of the genome, the target/genome sequence will be reverse complemented, and the guide returned
  * as is.
  *
  * @param guide the guide sequence (including the PAM if given)
  * @param chrom the chromosome to which the guide was aligned
  * @param start the first base on the chromosome to which the guide is aligned (always < end)
  * @param end the last base on the chromosome to which the guide is aligned
  * @param strand the strand of the chromosome to which the guide is aligned
  * @param score the score of the alignment
  * @param paddedGuide the padded guide sequence
  * @param paddedAlignment a string the same length as the paddedGuide and paddedTarget which represents
  *                        whether each base is a match/mismatch or part of a gap
  * @param paddedTarget the padded target sequence
  */
case class GuideAlignment(guide: String,
                          chrom: String,
                          start: Int,
                          end: Int,
                          strand: Char,
                          score: Int,
                          paddedGuide: String,
                          paddedAlignment: String,
                          paddedTarget: String
                         ) {
  require(paddedGuide.length  == paddedAlignment.length, "Padded guide and alignment string are different lengths.")
  require(paddedTarget.length == paddedAlignment.length, "Padded target and alignment string are different lengths.")

  /* For debugging only - print the alignment to stdout. */
  def print(out: PrintStream = System.out): Unit = Seq(paddedGuide, paddedAlignment, paddedTarget).foreach(out.println)

  def mismatches: Int = paddedAlignment.count(_ == '.')
  def gapBases:   Int = paddedAlignment.count(_ == ' ')
  def edits:      Int = paddedAlignment.count(ch => ch == '.' || ch == ' ')
}


object Aligner {
  /** EDNAfull matrix extended to also include 'U' as equivalent to T. */
  private val EdnaFull =
    """
      |    A   T    U    G    C    S   W   R   Y   K   M   B   V   H   D   N
      |A   5  -4   -4   -4   -4   -4   1   1  -4  -4   1  -4  -1  -1  -1  -2
      |T  -4   5    5   -4   -4   -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
      |U  -4   5    5   -4   -4   -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
      |G  -4  -4   -4    5   -4    1  -4   1  -4   1  -4  -1  -1  -4  -1  -2
      |C  -4  -4   -4   -4    5    1  -4  -4   1  -4   1  -1  -1  -1  -4  -2
      |S  -4  -4   -4    1    1   -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1
      |W   1   1    1   -4   -4   -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1
      |R   1  -4   -4    1   -4   -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1
      |Y  -4   1    1   -4    1   -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1
      |K  -4   1    1    1   -4   -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1
      |M   1  -4   -4   -4    1   -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1
      |B  -4  -1   -1   -1   -1   -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1
      |V  -1  -4   -4   -1   -1   -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1
      |H  -1  -1   -1   -4   -1   -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1
      |D  -1  -1   -1   -1   -4   -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1
      |N  -2  -2   -2   -2   -2   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
      |a  12 -12  -12  -12  -12   -4   1   1  -4  -4   1  -4  -1  -1  -1  -2
      |t -12  12   12  -12  -12   -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
      |u -12  12   12  -12  -12   -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
      |g -12 -12  -12   12  -12    1  -4   1  -4   1  -4  -1  -1  -4  -1  -2
      |c -12 -12  -12  -12   12    1  -4  -4   1  -4   1  -1  -1  -1  -4  -2
      |s -12 -12  -12   12   12   -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1
      |w  12  12   12  -12  -12   -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1
      |r  12 -12  -12   12  -12   -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1
      |y -12  12   12  -12   12   -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1
      |k -12  12   12   12  -12   -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1
      |m  12 -12  -12  -12   12   -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1
      |b -12  12   12   12   12   -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1
      |v  12 -12  -12   12   12   -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1
      |h  12  12   12  -12   12   -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1
      |d  12  12   12   12  -12   -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1
      |n  12  12   12   12   12   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
    """.stripMargin.lines.map(_.trim).filter(_.nonEmpty).toSeq

  /** Penalty charged to open a gap (extend penalty is also charged for first base in the gap). */
  val DefaultGapOpenPenalty: Int = 2

  /** Penalty charged to extend a gap, including the first base in the gap. */
  val DefaultGapExtendPenalty: Int = 5

  @inline private final def toOffset(ch: Char): Int = ch - 'A'
  @inline private final def toOffset(ch: Byte): Int = ch - 'A'

  // Parse the above matrix into a query-able 2D matrix using the offset from 'A' as the indices
  private val matrix = {
    val cols = EdnaFull.head.filter(_.isLetter).map(toOffset)
    val targetDim = cols.max + 1
    val queryDim  = EdnaFull.map(_.charAt(0)).filter(_.isLetter).map(toOffset).max + 1

    val m = new LinearMatrix[Int](queryDim, targetDim)
    for (x <- 0 until queryDim; y <- 0 until targetDim) m(x,y) = Int.MinValue

    EdnaFull.tail.foreach { line =>
      val xs     = line.split("""\s+""")
      val query  = toOffset(xs.head.charAt(0))
      val scores = xs.tail.map(_.toInt)

      forloop(from=0, until=scores.length) { i =>
        val target = cols(i)
        val score  = scores(i)
        m(query, target) = score
      }
    }

    m
  }

  /** Returns the score from the EDNA full matrix, or throws an exception if invalid bases are passed. */
  def score(query: Byte, target: Byte): Int = {
    val x =this.matrix(toOffset(query), toOffset(target))
    if (x == Int.MinValue) throw new IllegalArgumentException(s"${query.toChar}/${target.toChar} not present in scoring matrix.")
    else x
  }
}

/**
  * Class that performs alignments of guide sequences (or really any short sequence) against a small section
  * of a chromosome from a reference genome.  The alignment is performed using a modified Needleman-Wunch aligner
  * which performs a Glocal (i.e. semi-local or global/local) alignment.
  *
  * @param refFile the reference genome which target sequences are to be drawn from
  */
class Aligner(refFile: ReferenceSequenceFile,
              gapOpen: Int = Aligner.DefaultGapOpenPenalty,
              gapExtend: Int = Aligner.DefaultGapExtendPenalty) {
  require(refFile.isIndexed, s"Cannot work with a non-indexed reference: $refFile")

  // The actual aligner
  private val nw = new com.fulcrumgenomics.alignment.Aligner(Aligner.score, gapOpen= -abs(gapOpen), gapExtend= -abs(gapExtend), mode=Glocal)

  /**
    * Aligns a guide sequence to a region around the given position.  Attempts alignment on both the F and R
    * strands and returns the best scoring alignment of the pair.
    *
    * @param guide the guide sequence; must be only A/C/G/T and N.
    * @param chrom the chromosome on which the putative hit is found
    * @param pos the approximate location at which to align the guide
    * @return the best alignment of the guide to either the F or R strand
    */
  def align(guide: String, chrom: String, pos: Int): GuideAlignment = {
    val ref = refFile.getSequenceDictionary.getSequence(chrom)
    require(ref != null, s"Unknown chromosome: $chrom")

    val padding = guide.length * 2
    val (regionStart, regionEnd) = (max(pos-padding, 1), min(pos+padding, ref.getSequenceLength))
    val targetFwd = refFile.getSubsequenceAt(chrom, regionStart, regionEnd).getBases
    SequenceUtil.upperCase(targetFwd)
    val targetRev = targetFwd.clone()
    SequenceUtil.reverseComplement(targetRev)

    val query = guide.getBytes

    val fwd = this.nw.align(query, targetFwd)
    val rev = this.nw.align(query, targetRev)
    val best = if (fwd.score >= rev.score) fwd else rev
    val (start, end, strand) =
      if (best eq fwd) (regionStart + fwd.targetStart - 1, regionStart + fwd.targetEnd - 1, 'F')
      else             (regionEnd - rev.targetEnd + 1,     regionEnd - rev.targetStart + 1, 'R')

    val Seq(paddedGuide, alignString, paddedTarget) = best.paddedString()

    // Regenerate the align string to account for Us
    val buffer = new Array[Byte](paddedGuide.length)
    forloop (from=0, until=buffer.length) { i =>
      val q = paddedGuide(i)
      val t = paddedTarget(i)
      buffer(i) =
        if (q == '-' || t == '-') ' '
        else {
          val q2 = if (q == 'U') 'T' else q
          val t2 = if (t == 'U') 'T' else t
          if (SequenceUtil.readBaseMatchesRefBaseWithAmbiguity(t2.toByte, q2.toByte)) '|' else '.'
        }
    }

    GuideAlignment(
      guide           = guide,
      chrom           = chrom,
      start           = start,
      end             = end,
      strand          = strand,
      score           = best.score,
      paddedGuide     = paddedGuide,
      paddedAlignment = new String(buffer),
      paddedTarget    = paddedTarget
    )
  }
}
