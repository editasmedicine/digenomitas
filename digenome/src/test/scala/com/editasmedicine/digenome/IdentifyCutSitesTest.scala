package com.editasmedicine.digenome

import com.editasmedicine.commons.testing.UnitSpec
import com.editasmedicine.digenome.IdentifyCutSites.{CutSiteAccumulator, CutSiteParams}
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder}
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.Interval

import scala.util.Random

class IdentifyCutSitesTest extends UnitSpec {
  /** A set of default parameters/filters for calling cut sites, which can be adjusted easliy using `copy()`. */
  private val DefaultParams = CutSiteParams(
    readLength            = 100,
    insertSize            = 300,
    guide                 = Some("CTATTTCTCGATCGATCGAT"),
    enzyme                = Some("SpyCas9"),
    pamFivePrime          = None,
    pamThreePrime         = Some("NGG"),
    overhang              = 0,
    maxOffset             = 2,
    minDepth              = 10,
    maxDepth              = 300,
    maxLowMapqFraction    = 0.2,
    minForwardReads       = 2,
    minReverseReads       = 2,
    minSupportingReads    = 5,
    minSupportingFraction = 0.2
  )
  
  // Generate a reference fasta, index, etc. with known sequence.
  private val (refFile, ref, dict) = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1").add("TACGTCGATTGTGCTAGCTGATATTTCGGCTATTTCTCGATCGATCGATCGGCGCTATAAAGCTCTTTAGCGCTAGCTCGATTCTAGCTCTTAGATATCG", 100)
    builder.add("chr2").add("CGATATCTAAGAGCTAGAATCGAGCTAGCGCTAAAGAGCTTTATAGCGCCGATCGATCGATCGAGAAATAGCCGAAATATCAGCTAGCACAATCGACGTA", 100)
    builder.add("chr3")
      .add("N", 1000) //1000
      .add("A", 500)
      .add("N", 100)
      .add("T", 400)
      .add("N", 600)  // 2000
      .add("C", 400)
      .add("N", 1000) // 3000
    val path = builder.toTempFile()
    val ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(path)
    (path, ref, ref.getSequenceDictionary)
  }

  /** Generates a SamBuilder for building up SamRecords for testing. */
  private def newBuilder(readLength: Int): SamBuilder = {
    new SamBuilder(readLength=readLength, sd=Some(dict), sort=Some(SamOrder.Coordinate))
  } 

  "IdentifyCutSites.filtered" should "filter out reads from non-FR mapped pairs" in {
    val builder = newBuilder(readLength=25)
    builder.addFrag("frag1", start=1000)
    builder.addFrag("frag2", unmapped=true)
    builder.addPair("badpair1", start1=1000, start2=1200, strand1=Plus, strand2=Plus)
    builder.addPair("badpair2", start1=1000, start2=1000, unmapped2=true)
    builder.addPair("badpair3", start1=1000, start2=1200) match { case Seq(r1, r2) =>
        r1.mateRefIndex = 1
        r2.refIndex = 1
    }
    builder.addPair("badpair4", start1=1000, start2=9000)
    builder.addPair("ok1", start1=1000, start2=1400)

    val filtered = IdentifyCutSites.filtered(builder.iterator).toIndexedSeq
    filtered should have size 2
    filtered.map(_.name).distinct shouldBe Seq("ok1")
  }

  it should "filter out duplicate reads" in {
    val builder = newBuilder(readLength=25)
    builder.addPair(start1=100, start2=200)

    IdentifyCutSites.filtered(builder.iterator).toSeq should have size 2
    IdentifyCutSites.filtered(builder.iterator.map{r => r.duplicate = true; r}).toSeq should have size 0
  }

  it should "filter out secondary and supplementary reads" in {
    val builder = newBuilder(readLength=25)
    builder.addPair(start1=100, start2=200)

    IdentifyCutSites.filtered(builder.iterator).toSeq should have size 2
    IdentifyCutSites.filtered(builder.iterator.map{r => r.secondary = true; r}).toSeq should have size 0
    IdentifyCutSites.filtered(builder.iterator.map{r => r.supplementary= true; r}).toSeq should have size 0
  }

  it should "work on an empty iterator" in {
    val filtered = IdentifyCutSites.filtered(Iterator.empty)
    filtered.hasNext shouldBe false
  }

  "IdentifyCutSites.determineReadLengthAndInsertSize" should "correctly calculate values when all reads are the same length" in {
    val builder = newBuilder(readLength=100)
    builder.addPair(start1=1000, start2=1200)
    builder.addPair(start1=1000, start2=1300)
    builder.addPair(start1=1000, start2=1400)
    builder.addPair(start1=1000, start2=8000) // should be filtered out

    val (rlen, isize) = IdentifyCutSites.determineReadLengthAndInsertSize(Seq(new SelfClosingIterator(builder.iterator, () => Unit)))
    rlen shouldBe 100
    isize shouldBe 400
  }

  it should "calculate read length correctly when there is a mix of read lengths" in {
    val builder1 = newBuilder(readLength=100)
    val builder2 = newBuilder(readLength=200)
    builder1.addPair(start1=1000, start2=1200) // isize = 300
    builder1.addPair(start1=1000, start2=1300) // isize = 400
    builder2.addPair(start1=1000, start2=1200) // isize = 400
    builder2.addPair(start1=1000, start2=1300) // isize = 500
    builder2.addPair(start1=1100, start2=1400) // isize = 500
    builder2.addPair(start1=1100, start2=1500) // isize = 600

    val iter = new SelfClosingIterator(builder1.iterator ++ builder2.iterator, () => Unit)
    val (rlen, isize) = IdentifyCutSites.determineReadLengthAndInsertSize(Seq(iter))
    rlen shouldBe (4*100 + 8*200) / 12
    isize shouldBe 450
  }

  "IdentifyCutSites.wholeGenomeIntervalList" should "cut on large blocks of Ns" in {
    val ilist = IdentifyCutSites.wholeGenomeIntervalList(this.refFile, 500)
    val regions  = ilist.getIntervals.toIndexedSeq.map(i => new Interval(i.getContig, i.getStart, i.getEnd)) // drop names for ease of comparison
    val expected = IndexedSeq(
      new Interval("chr1", 1, 10000),
      new Interval("chr2", 1, 10000),
      new Interval("chr3", 1001, 2000),
      new Interval("chr3", 2601, 3000)
    )

    regions should contain theSameElementsInOrderAs expected
  }

  "CutSiteAccumulator" should "ignore unmapped reads" in {
    val builder = newBuilder(readLength=20)
    builder.addPair(unmapped1=true, unmapped2=true)
    builder.addPair(unmapped1=true, start2=1100)
    builder.addPair(start1=1100, unmapped2=true)
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    acc ++= builder.iterator.filter(_.unmapped)

    acc.totalFwdStarts shouldBe 0
    acc.totalRevStarts shouldBe 0
    acc.totalReadDepth shouldBe 0
    acc.totalTemplateDepth shouldBe 0
    acc.totalLowMapqReads shouldBe 0
  }

  it should "accumulate a valid FR read pair as expected" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    acc ++= newBuilder(readLength=20).addPair(start1=1100, start2=1100, strand1=Plus, strand2=Minus, mapq1=60, mapq2=60)

    acc.totalFwdStarts shouldBe 1
    acc.totalRevStarts shouldBe 1
    acc.fStarts(1100) shouldBe 1
    acc.rStarts(1119) shouldBe 1

    Range.inclusive(1000, 2000).foreach { pos =>
      acc.lowMapqFraction(pos) shouldBe 0.0
      if (pos >= 1100 && pos <= 1119) {
        acc.readDepth(pos) shouldBe 2
        acc.templateDepth(pos) shouldBe 1
      }
      else {
        acc.readDepth(pos) shouldBe 0
        acc.templateDepth(pos) shouldBe 0
      }
    }
  }

  it should "correctly count read depth and template depth based on aligned start/end and include gaps" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    acc ++= newBuilder(readLength=20).addPair(start1=1100, start2=1271, strand1=Plus, strand2=Minus, cigar1="20M", cigar2="10M10D10M")
    val r1Span       = Range.inclusive(1100, 1119)
    val r2Span       = Range.inclusive(1271, 1300)
    val templateSpan = Range.inclusive(1100, 1300)

    acc.totalFwdStarts shouldBe 1
    acc.totalRevStarts shouldBe 1
    acc.fStarts(1100) shouldBe 1
    acc.rStarts(1300) shouldBe 1

    Range.inclusive(1000, 2000).foreach { pos =>
      acc.lowMapqFraction(pos) shouldBe 0.0
      if (r1Span.contains(pos) || r2Span.contains(pos))
        acc.readDepth(pos) shouldBe 1
      else
        acc.readDepth(pos) shouldBe 0

      if (templateSpan.contains(pos))
        acc.templateDepth(pos) shouldBe 1
      else
        acc.templateDepth(pos) shouldBe 0
    }
  }

  it should "count reads with clipping at the 5' end of the read for depth but not starts" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    acc ++=   newBuilder(readLength=20).addPair(start1=1100, start2=1180, strand1=Plus, strand2=Minus, cigar1="5S15M", cigar2="18M2S")
    val r1Span       = Range.inclusive(1100, 1114)
    val r2Span       = Range.inclusive(1180, 1197)
    val templateSpan = Range.inclusive(1100, 1197)

    acc.totalFwdStarts shouldBe 0
    acc.totalRevStarts shouldBe 0

    Range.inclusive(1000, 2000).foreach { pos =>
      acc.lowMapqFraction(pos) shouldBe 0.0
      if (r1Span.contains(pos) || r2Span.contains(pos))
        acc.readDepth(pos) shouldBe 1
      else
        acc.readDepth(pos) shouldBe 0

      if (templateSpan.contains(pos))
        acc.templateDepth(pos) shouldBe 1
      else
        acc.templateDepth(pos) shouldBe 0
    }
  }

  it should "handle reads that map entirely outside of the region" in {
    val acc     = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=20)
    acc ++= builder.addPair(start1=500, start2=980)   // last base of R read is at 999
    acc ++= builder.addPair(start1=2001, start2=2100) // first base of F read is at 2001
    acc.totalFwdStarts shouldBe 0
    acc.totalRevStarts shouldBe 0
    acc.totalReadDepth shouldBe 0
    acc.totalTemplateDepth shouldBe 0

    acc ++= builder.addPair(start1=500, start2=2500) // template straddle region but reads don't overlap
    acc.totalReadDepth shouldBe 0
    acc.totalTemplateDepth shouldBe 1001 // 1 at each position
  }

  it should "handle reads that only partially overlap the region" in {
    val acc     = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=20)
    acc ++= builder.addPair(start1=990, start2=2000)
    acc.totalFwdStarts shouldBe 0
    acc.totalRevStarts shouldBe 0
    acc.totalReadDepth shouldBe 11 // 10bp overlap from R1 and 1bp overlap from R2

    acc ++= builder.addPair(start1=1999, start2=2200)
    acc ++= builder.addPair(start1=950,  start2=991)
    acc.totalFwdStarts shouldBe 1
    acc.totalRevStarts shouldBe 1
    acc.fStarts(1999) shouldBe 1
    acc.rStarts(1010) shouldBe 1
  }

  it should "count low-mapq reads towards depth, but not read starts" in {
    val builder = newBuilder(readLength=20)
    builder.addPair(start1=1100, start2=1100, mapq1=19, mapq2=19)
    builder.addPair(start1=1100, start2=1100, mapq1=0,  mapq2=0)

    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    acc ++= builder

    acc.totalFwdStarts shouldBe 0
    acc.totalRevStarts shouldBe 0
    Range.inclusive(1000, 2000).foreach { pos =>
      acc.readDepth(pos)       shouldBe 0
      acc.templateDepth(pos)   shouldBe 0
      if (pos >= 1100 && pos <= 1119)
        acc.lowMapqFraction(pos) shouldBe 1.0
      else
        acc.lowMapqFraction(pos) shouldBe 0.0
    }

    acc ++= builder.addPair(start1=1200, start2=1200, mapq1=20, mapq2=20)
    acc.fStarts(1200) shouldBe 1
    acc.rStarts(1219) shouldBe 1
  }


  "CutSiteAccumulator.median" should "return the index of the median bin" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 1100, minMapQ=20, ref=ref)
    acc.median(Seq(5,5,5), 15) shouldBe 1
    acc.median(Seq(6,5), 11) shouldBe 0
    acc.median(Seq(5,5), 10) shouldBe 0 // we don't do half-bins, so it'll round down
    acc.median(Seq(1,2,3,4,5), 15) shouldBe 3
  }

  it should "return the middle index when the total is 0" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 1100, minMapQ=20, ref=ref)
    acc.median(Seq(0), 0) shouldBe 0
    acc.median(Seq(0,0,0), 0) shouldBe 1
    acc.median(Seq(0,0,0,0,0), 0) shouldBe 2
  }


  "CutSiteAccumulator.ok" should "correctly filter" in {
    val params = DefaultParams.copy(
      minDepth              = 10,
      maxDepth              = 300,
      maxLowMapqFraction    = 0.1,
      minForwardReads       = 2,
      minReverseReads       = 2,
      minSupportingReads    = 5,
      minSupportingFraction = 0.2
    )

    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 1100, minMapQ=20, ref=ref)
    acc.ok(fStarts=5, rStarts=5, depth=50,  filters=params) shouldBe true
    acc.ok(fStarts=2, rStarts=2, depth=20,  filters=params) shouldBe false // fails for minSupportingReads < 5
    acc.ok(fStarts=3, rStarts=3, depth=100, filters=params) shouldBe false // fails for minSupportingFraction < 0.2
    acc.ok(fStarts=1, rStarts=10, depth=15, filters=params) shouldBe false // fails for fStarts < minForwardReads
    acc.ok(fStarts=1, rStarts=1, depth=15, filters=params) shouldBe false // fails for rStarts < minReverseReads
  }


  "CutSiteAccumulator.pValue" should "calculate meaningful p-values in most cases" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 1100, minMapQ=20, ref=ref)
    acc.pValue(starts=1, depth=100, p=1/150d) shouldBe 0.48 +- 0.01
    acc.pValue(starts=2, depth=100, p=1/150d) shouldBe 0.14 +- 0.01
    acc.pValue(starts=3, depth=100, p=1/150d) shouldBe 0.03 +- 0.01

    Range.inclusive(0, 50).map(acc.pValue(_, 50, p=5/150d)).sliding(2).foreach { case Seq(lhs, rhs) =>
      lhs > 0 shouldBe true
      rhs > 0 shouldBe true
      lhs should be >= rhs
    }
  }


  "CutSiteAccumulator.build" should "generate no cut sites when no sites have sufficient starts/stops" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=100)
    Range.inclusive(900, 2100, step=10).foreach { p => builder.addPair(start1=p, start2=p+200) }
    acc ++= builder

    Range.inclusive(1000, 2000).foreach { pos => acc.build(pos, DefaultParams) shouldBe None }
  }

  it should "generate a single cut site on the F strand at position 1500 with overhang=0 and maxOffset=1" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=100)
    forloop(from=0, until=20) { i => builder.addPair(start1=1200, start2=1401) }
    forloop(from=0, until=7)  { i => builder.addPair(start1=1500, start2=1700) }
    forloop(from=0, until=8)  { i => builder.addPair(start1=1501, start2=1704) }
    forloop(from=0, until=9)  { i => builder.addPair(start1=1502, start2=1711) }

    acc ++= builder
    acc.build(1500, DefaultParams.copy(overhang=0, maxOffset=1)) match {
      case None => fail("Failed to generate a cut site at pos=1500.")
      case Some(cut) =>
        cut.strand shouldBe "F" // guide "should" match F strand
        cut.forward_starts        shouldBe 24
        cut.reverse_starts        shouldBe 20
        cut.read_depth            shouldBe 44
        cut.template_depth        shouldBe 44
        cut.read_fraction_cut     shouldBe 1.0
        cut.template_fraction_cut shouldBe 1.0
        cut.median_overhang       shouldBe 0
        cut.overhang_distribution shouldBe "7;8;9"
    }
  }

  it should "generate a single cut site on the R strand at position 1500 with overhang=0 and maxOffset=1" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=100)
    forloop(from=0, until= 9) { i => builder.addPair(start1=1200, start2=1399) } // end at 1498
    forloop(from=0, until=10) { i => builder.addPair(start1=1200, start2=1400) } // end at 1499
    forloop(from=0, until=11) { i => builder.addPair(start1=1200, start2=1401) } // end at 1500
    forloop(from=0, until=31) { i => builder.addPair(start1=1500, start2=1700) }

    acc ++= builder
    acc.build(1500, DefaultParams.copy(overhang=0, maxOffset=1)) match {
      case None => fail("Failed to generate a cut site at pos=1500.")
      case Some(cut) =>
        cut.strand shouldBe "R" // guide "should" match R strand
        cut.forward_starts        shouldBe 31
        cut.reverse_starts        shouldBe 30
        cut.read_depth            shouldBe 61
        cut.template_depth        shouldBe 61
        cut.read_fraction_cut     shouldBe 1.0
        cut.template_fraction_cut shouldBe 1.0
        cut.median_overhang       shouldBe 0
        cut.overhang_distribution shouldBe "11;10;9"
    }
  }

  it should "respect the max-offset for F strand cut sites" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=100)
    forloop(from=0, until=20) { i => builder.addPair(start1=1200, start2=1401) }

    // Add 8 reads each at the exact spot and at each position +/- 3bp
    for (off <- -3 to 3; idx <- 1 to 8) { builder.addPair(start1=1501+off, start2=1700+off) }

    acc ++= builder
    (0 to 3) foreach { offset =>
      acc.build(1500, DefaultParams.copy(overhang=0, maxOffset=offset)) match {
        case None      => fail("Failed to generate a cut site at pos=1500.")
        case Some(cut) =>
          cut.strand shouldBe "F" // guide "should" match F strand
          cut.forward_starts        shouldBe (8 + (8*2*offset))
          cut.reverse_starts        shouldBe 20
      }
    }
  }

  it should "respect the max-offset for R strand cut sites" in {
    val acc = new CutSiteAccumulator("s1", "chr1", 1000, 2000, minMapQ=20, ref=ref)
    val builder = newBuilder(readLength=100)

    // Add 8 reads each at the exact spot and at each position +/- 3bp
    for (off <- -3 to 3; idx <- 1 to 8) { builder.addPair(start1=1200+off, start2=1400+off) } // end = 1499 + offset

    forloop(from=0, until=20) { i => builder.addPair(start1=1500, start2=1708) }

    acc ++= builder
    (0 to 3) foreach { offset =>
      acc.build(1500, DefaultParams.copy(overhang=0, maxOffset=offset)) match {
        case None      => fail("Failed to generate a cut site at pos=1500.")
        case Some(cut) =>
          cut.strand shouldBe "R" // guide "should" match R strand
          cut.forward_starts        shouldBe 20
          cut.reverse_starts        shouldBe (8 + (8*2*offset))
      }
    }
  }
  
  "IdentifyCutSites" should "fail if all the min-reads parameters are set to 0" in {
    an[ValidationException] shouldBe thrownBy {
      val path = makeTempFile("foo.", ".bam")
      new IdentifyCutSites(input=Seq(path), output=path, guide=Some("CTATTTCTCGATCGATCGAT"), pamThreePrime=Some("NGG"), enzyme=Some("SpyCas9"),
        ref=refFile, minForwardReads=0, minReverseReads=0, minSupportingReads=0)
    }
  }

  Range.inclusive(1, 2).foreach { bamCount =>
   it should s"run end to end on a 'whole genome' across $bamCount BAM files and generate some cut sites" in {
     val builders = (1 to bamCount) map { _ => newBuilder(readLength=100) }
     val random   = new Random(42)
     def builder  = builders(random.nextInt(bamCount))

     // Evenly spaced reads up to around 2000
     for (pos <- Range.inclusive(start=200,end=2000,step=5)) {
       builder.addPair(start1=pos-199, start2=pos-105)
     }

     // Evenly spaced reads from around 2200
     for (pos <- Range.inclusive(start=2300,end=4000,step=5)) {
       builder.addPair(start1=pos, start2=pos+110)
     }

     // A bunch of stuff for an F-strand cut site at 2150
     //  - Nearly all R reads on site, but with a tiny wobble
     //  - A good spread of F reads around the expected site (with overhang = -2)
     builder.addPair(start1=1960, start2=2050)
     forloop (from=0, until=50) { _ => builder.addPair(start1=1950, start2=2051) }
     builder.addPair(start1=1960, start2=2052)

     forloop (from=0, until= 2) { _ => builder.addPair(start1=2151, start2=2267) } // -2
     forloop (from=0, until=12) { _ => builder.addPair(start1=2152, start2=2271) } // -1
     forloop (from=0, until=50) { _ => builder.addPair(start1=2153, start2=2264) } // on
     forloop (from=0, until=15) { _ => builder.addPair(start1=2154, start2=2283) } // +1
     forloop (from=0, until= 3) { _ => builder.addPair(start1=2155, start2=2277) } // +2

     // Double check that all BAMs got some reads
     builders.size shouldBe bamCount
     builders.forall(_.nonEmpty) shouldBe true

     val bams = builders.map(_.toTempFile())
     val out = makeTempFile("metrics.", ".txt")
     new IdentifyCutSites(input=bams, output=out, guide=Some("CTATTTCTCGATCGATCGAT"), pamThreePrime=Some("NGG"), enzyme=Some("SpyCas9"),
       ref=refFile, overhang= -2, maxOffset=1, minDepth=5).execute()
     val metrics = Metric.read[CutSiteInfo](out)

     metrics.size shouldBe 1
     val cut = metrics.head
     cut.chrom shouldBe "chr1"
     cut.pos   shouldBe 2149 // output is 0-based
     cut.strand shouldBe "F"
     cut.forward_starts shouldBe 12+50+15
     cut.reverse_starts shouldBe 50
     cut.read_depth shouldBe 128 // the starts, plus one rogue F read that's one base too far away
     cut.read_fraction_cut shouldBe     127/128d +- 0.0001
     cut.template_fraction_cut shouldBe 127/128d +- 0.0001
     cut.median_overhang shouldBe -2
     cut.overhang_distribution shouldBe "12;50;15"
   }
  }
}
