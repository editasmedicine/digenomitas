package com.editasmedicine.digenome

import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.fastq.{FastqRecord, FastqSource, FastqWriter}

class InterleaveFastqsTest extends UnitSpec {

  def fq(name: String, readNumber: Int, length: Int = 10): FastqRecord = {
    FastqRecord(name, "A"*length, "#"*length, comment = None, readNumber = Some(readNumber))
  }

  // A whole bunch of small Seqs of FastqRecord
  private val r1a = (1 to 5).map(i => fq("a-" + i, 1))
  private val r2a = (1 to 5).map(i => fq("a-" + i, 2))
  private val r1b = (1 to 4).map(i => fq("b-" + i, 1))
  private val r2b = (1 to 4).map(i => fq("b-" + i, 2))
  private val r1c = (1 to 6).map(i => fq("c-" + i, 1))
  private val r2c = (1 to 6).map(i => fq("c-" + i, 2))
  private val r1d = (1 to 5).map(i => fq("d-" + i, 1))
  private val r2d = (1 to 5).map(i => fq("d-" + i, 2))

  // Handy implicit to build a buffered iterator from a Seq
  private implicit def seqToBufferedIterator(fs: Seq[FastqRecord]): BetterBufferedIterator[FastqRecord] = fs.iterator.bufferBetter

  "FqPairIterator" should "work correctly with just a single pair of iterators" in {
    val iterator = new FqPairIterator(Seq(r1a), Seq(r2a))
    val pairs = iterator.toSeq
    pairs should have size 5
    pairs.foreach { case (r1, r2) =>
      r1.name shouldBe r2.name
      r1.readNumber == r2.readNumber shouldBe false
    }
  }

  it should "work with several iterators in order" in {
    val iterator = new FqPairIterator(Seq(r1a, r1b, r1c), Seq(r2a, r2b, r2c))
    val pairs = iterator.toSeq
    pairs should have size 5+4+6
    pairs.foreach { case (r1, r2) =>
      r1.name shouldBe r2.name
      r1.readNumber == r2.readNumber shouldBe false
    }
  }

  it should "work with several iterators out of order" in {
    val iterator = new FqPairIterator(Seq(r1a, r1b, r1c), Seq(r2b, r2c, r2a))
    val pairs = iterator.toSeq
    pairs should have size 5+4+6
    pairs.foreach { case (r1, r2) =>
      r1.name shouldBe r2.name
      r1.readNumber == r2.readNumber shouldBe false
    }
  }

  it should "fail if given different numbers of iterators" in {
    an[Exception] shouldBe thrownBy { new FqPairIterator(Seq(r1a), Seq(r2a, r2b)).toList }
  }

  it should "fail if given the same number, but mismatched, iterators" in {
    an[Exception] shouldBe thrownBy { new FqPairIterator(Seq(r1a), Seq(r2d)).toList }
  }

  it should "fail if given two iterators in different orders" in {
    an[Exception] shouldBe thrownBy { new FqPairIterator(Seq(r1a), Seq(r2a.reverse)).toList }
  }

  it should "fail if given iterators of different lengths" in {
    an[Exception] shouldBe thrownBy { new FqPairIterator(Seq(r1a), Seq(r2a.take(3))).toList }
  }

  "InterleaveFastqs" should "run end to end and generate valid output" in {
    def fqToFile(fs: Seq[FastqRecord]): FilePath = {
      val tmp = makeTempFile("test.", ".fastq")
      val out = FastqWriter(tmp)
      out ++= fs
      out.close()
      tmp
    }

    val r1s = Seq(r1a, r1b, r1c, r1d).map(fqToFile)
    val r2s = Seq(r2a, r2c, r2d, r2b).map(fqToFile)
    val out = makeTempFile("result.", ".fastq")
    new InterleaveFastqs(readOne=r1s, readTwo=r2s, output=out).execute()
    val results  = FastqSource(out).toIndexedSeq
    val expected = Seq(r1a, r1b, r1c, r1d).flatten
      .zip(Seq(r2a, r2b, r2c, r2d).flatten)
      .toIndexedSeq
      .flatMap { case (r1, r2) => IndexedSeq(r1, r2) }

    results should contain theSameElementsInOrderAs expected
  }
}
