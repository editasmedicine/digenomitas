package com.editasmedicine.digenome

import java.io.Closeable

import com.editasmedicine.aligner.Aligner
import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{QueryType, SamRecord, SamSource}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.NumericCounter
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Metric, ProgressLogger}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory, ReferenceSequenceFileWalker}
import htsjdk.samtools.util.{CoordMath, Interval, IntervalList}
import org.apache.commons.math3.distribution.PoissonDistribution

import scala.collection.mutable
import scala.math.{log10, max, min}

object IdentifyCutSites {
  /**
    * Simple case class to encapsulate all the parameters and filters used when building/filtering cut sites
    * from simple counts.
    *
    * @param readLength the average read length in the data
    * @param insertSize the median insert size of the data
    * @param guide the optional sequence of the guide used
    * @param pamFivePrime the optional sequence of any PAM at the 5' end of the guide
    * @param pamThreePrime the optional sequence of any PAM at the 3' end of the guide
    * @param enzyme optionally the enzyme used in the cutting reaction
    * @param overhang the expected overhang at the cut sites (positive for overhang, negative for gapped)
    * @param maxOffset the maximum offset from the expected cut site at which to count read starts of the opposite strand
    * @param minDepth the minimum number of reads covering a site in order to emit
    * @param maxDepth the maximum number of reads covering a site in order to emit
    * @param maxLowMapqFraction the maximum fraction of low-mapq reads at a site in order to emit
    * @param minForwardReads the minimum number of forward reads starting at the cut site in order to emit
    * @param minReverseReads the minimum number of reverse reads starting at the cut site in order to emit
    * @param minSupportingReads the minimum total number of reads starting at the cut site in order to emit
    * @param minSupportingFraction the minimum fraction of reads covering the site which start at the site
    */
  case class CutSiteParams(readLength: Int,
                           insertSize: Int,
                           guide: Option[String],
                           pamFivePrime: Option[String],
                           pamThreePrime: Option[String],
                           enzyme: Option[String],
                           overhang: Int,
                           maxOffset: Int,
                           minDepth: Int,
                           maxDepth: Int,
                           maxLowMapqFraction: Double,
                           minForwardReads: Int,
                           minReverseReads: Int,
                           minSupportingReads: Int,
                           minSupportingFraction: Double
                           ) {

    /** The guide with the PAM concatenated to it. Guide sequence is upper case.  PAM sequence is lower case. */
    val guidePlusPam: Option[String] = guide.map(pamFivePrime.getOrElse("").toLowerCase + _.toUpperCase + pamThreePrime.getOrElse("").toLowerCase)

    /** The probability of a read starting at a given site given that it covers the site. */
    val readStartProbability: Double  = (1 + Range.inclusive(-maxOffset, maxOffset).length) / (2*readLength).toDouble

    /** The probability of a read starting at a cut site given that the template covers the cut site. */
    val templateStartProbability:Double = (1 + Range.inclusive(-maxOffset, maxOffset).length) / insertSize.toDouble
  }

  /**
    * Class which uses a number of long arrays to track counts of various things (read starts, coverage)
    * efficiently.  Counted values are as follows:
    *
    * - fStarts: the number of reads on the forward strand starting at that position, which do not have clipping
    *            at the beginning of the read/alignment
    * - rStarts: the number of reads on the reverse strand starting at that position (i.e. the first base read is
    *            at that position, so the last base in the alignment is at that position) which do not have
    *            clipping at the end of the alignment
    * - readDepth: the number of high mapq reads that cover this position, where cover is defined as having the
    *              position between the start and end of the read alignment, even if the position is gapped in
    *              the alignment
    * - templateDepth: the number of templates (i.e. read pairs) which cover this position, where cover is defined
    *                  as having the position between the first base in the template and the last, including any
    *                  un-sequenced portion in the middle. Excludes templates with low mapping quality.
    * - lowMapqReads: the number of reads covering the position with low mapping quality
    *
    * Important Notes:
    * 1. fStarts/rStarts and readDepth/templateDepth are NOT mutually exclusive.  I.e. at a given position where
    *    a forward strand read starts, we will add one to fStarts(pos) and one to readDepth(pos).
    * 2. readDepth/templateDepth and lowMapqReads ARE mutually exclusive.  Either a reason is of sufficient mapping
    *    quality and is counted in readDepth/templateDepth OR it is counted in lowMapqReads.
    *
    * @param sample the sample being accumulated for
    * @param chrom the chromosome on which the region being accumulated exists
    * @param start the start of the region over which to track counts
    * @param end the end of the region over which to track counts
    * @param minMapQ the minimum mapping quality for a read to be considered well mapped
    */
  class CutSiteAccumulator(val sample: String,
                           val chrom: String,
                           val start: Int,
                           val end: Int,
                           val minMapQ: Int,
                           val ref : ReferenceSequenceFile) {
    private val length: Int = end - start + 1
    private val aligner = new Aligner(ref)

    // A whole bunch of arrays the length of the interval for tracking various counts
    private val _fStarts       = new Array[Short](length)
    private val _rStarts       = new Array[Short](length)
    private val _readDepth     = new Array[Short](length)
    private val _templateDepth = new Array[Short](length)
    private val _lowMapqReads  = new Array[Short](length)

    /** Increments a Short value while capping instead of overflowing. */
    @inline private def safeIncrement(xs: Array[Short], pos: Int): Unit = {
      if (pos >= start && pos <= end && xs(pos-start) < Short.MaxValue) xs(pos-start) = (xs(pos-start) + 1).toShort
    }

    // Private methods for incrementing counts in the various arrays, being careful not to overflow Short's limits
    @inline private def incrementForwardStarts(pos: Int): Unit = safeIncrement(_fStarts, pos)
    @inline private def incrementReverseStarts(pos: Int): Unit = safeIncrement(_rStarts, pos)
    @inline private def incrementReadDepth(pos: Int): Unit     = safeIncrement(_readDepth, pos)
    @inline private def incrementTemplateDepth(pos: Int): Unit = safeIncrement(_templateDepth, pos)
    @inline private def incrementLowMapqReads(pos: Int): Unit  = safeIncrement(_lowMapqReads, pos)

    /** Returns the value from the short array if the position is within (start,end), else returns 0. */
    @inline private def safeFetch(xs: Array[Short], pos: Int): Short = {
      if (pos >= start && pos <= end) xs(pos - start) else 0
    }

    // Methods for fetching counts by genomic position within the range.
    @inline final def fStarts(pos: Int): Int       = safeFetch(_fStarts, pos)
    @inline final def rStarts(pos: Int): Int       = safeFetch(_rStarts, pos)
    @inline final def readDepth(pos: Int): Int     = safeFetch(_readDepth, pos)
    @inline final def templateDepth(pos: Int): Int = safeFetch(_templateDepth, pos)
    @inline final def poorlyMapped(pos: Int): Int  = safeFetch(_lowMapqReads, pos)
    final def lowMapqFraction(pos: Int): Double = {
      val lowMapq = poorlyMapped(pos)
      val total   = poorlyMapped(pos) + readDepth(pos)
      if (total == 0) 0 else lowMapq / total.toDouble
    }

    // Methods, mostly for testing, to get totals of various counts
    private[digenome] def totalFwdStarts:     Long = _fStarts.foldLeft(0L)(_ + _)
    private[digenome] def totalRevStarts:     Long = _rStarts.foldLeft(0L)(_ + _)
    private[digenome] def totalReadDepth:     Long = _readDepth.foldLeft(0L)(_ + _)
    private[digenome] def totalTemplateDepth: Long = _templateDepth.foldLeft(0L)(_ + _)
    private[digenome] def totalLowMapqReads:  Long = _lowMapqReads.foldLeft(0L)(_ + _)

    /** Method to accumulate counts from a single SamRecord.  Method should only be passed
      * "filtered" reads, excluding any reads that are not from an FR mapped pair or whose
      * insert sizes are anomolous.
      *
      * @param rec a SamRecord that is part of a reasonable FR pair
      */
    def accumulate(rec: SamRecord): Unit = if (rec.mapped) {
        // Determine if the read is overclipped, and what the read start site is based on strand
        val clipped      = if (rec.positiveStrand) rec.cigar.head.operator.isClipping else rec.cigar.last.operator.isClipping
        val fivePrimeEnd = if (rec.positiveStrand) rec.start else rec.end

        // Count the read coverage appropriately
        val lowMapQ  = rec.mapq < this.minMapQ

        forloop (from=max(start, rec.start), until=min(end, rec.end)+1) { i =>
          if (lowMapQ) incrementLowMapqReads(i) else incrementReadDepth(i)
        }

        // Count the template coverage
        if (!lowMapQ && rec.firstOfPair) {
          val (insStart, insEnd) = if (rec.positiveStrand) (rec.start, rec.start + rec.insertSize - 1) else (rec.end + rec.insertSize + 1, rec.end)
          forloop (from=max(start, insStart), until=min(end, insEnd)+1) { i => incrementTemplateDepth(i) }
        }

        // Count the read start
        if (!clipped && rec.mapq >= this.minMapQ) {
          if      (rec.positiveStrand && rec.start >= start) incrementForwardStarts(rec.start)
          else if (rec.negativeStrand && rec.end   <= end)   incrementReverseStarts(rec.end)
        }
    }

    /** Synonym for accumulate(). */
    def +=(rec: SamRecord): Unit = accumulate(rec)

    /** Accumulate multiple records. */
    def ++=(recs: TraversableOnce[SamRecord]): Unit = recs.foreach(accumulate)

    /**
      * Builds a cut site from the data accumulated at the given position.  If the cut site fails any of the
      * filtering critera, None will be returned, otherwise Some(CutSiteInfo).
      *
      * @param pos the position at which to build the [[CutSiteInfo]]
      * @param params the parameters and filters to use when building the cut site
      * @return either None, or Some(CutSiteInfo)
      */
    def build(pos: Int, params: CutSiteParams): Option[CutSiteInfo] = {
      if (readDepth(pos) < params.minDepth || readDepth(pos) > params.maxDepth || lowMapqFraction(pos) > params.maxLowMapqFraction) {
        None
      }
      else {
        // Construct a neighborhood measure of start site frequency vs. depth. The measure will be approximately
        // 1 when the rate of start sites vs. depth is expected/random, and approximately 0 when there is coverage
        // but no read starts.
        lazy val nearby = {
          var totalDepth  = 0
          var totalStarts = 0
          val startDistance = math.abs(params.overhang) + params.maxOffset + 2

          forloop(from=startDistance, until=startDistance+75) { i =>
            totalDepth  += templateDepth(pos+i)
            totalStarts += fStarts(pos+i)
            totalDepth  += templateDepth(pos-i)
            totalStarts += rStarts(pos-i)
          }

          totalStarts / totalDepth.toDouble * params.insertSize
        }

        // Lazily generate an alignment of the guide to the genome at this position
        lazy val alignment = params.guidePlusPam.map(query => this.aligner.align(query, chrom, pos))

        // fwd = Guide matches fwd sequence, binds to rev strand, yielding clean cut on rev and variable cut on fwd strand
        val fwd = {
          val fPos           = pos + 1 - params.overhang
          val fRange         = Range.inclusive(fPos - params.maxOffset, fPos + params.maxOffset)
          val fs             = fRange.map(fStarts)
          val fStartCount    = fs.sum
          val medianOverhang = median(fs.reverse, fStartCount) - params.maxOffset + params.overhang
          val rStartCount    = rStarts(pos)
          val bonusDepth     = fRange.iterator.filter(_ > pos).map(this.fStarts).sum
          val depth          = readDepth(pos) + bonusDepth
          val tDepth         = templateDepth(pos) + bonusDepth

          if (!ok(fStartCount, rStartCount, depth, params)) None else Some(CutSiteInfo(
            sample                = sample,
            guide                 = params.guide.getOrElse(""),
            enzyme                = params.enzyme.getOrElse(""),
            expected_overhang     = params.overhang,
            window_size           = params.maxOffset * 2 + 1,
            mapq_cutoff           = this.minMapQ,
            chrom                 = chrom,
            pos                   = pos - 1,
            strand                = "F",
            low_mapq_fraction     = lowMapqFraction(pos),
            forward_starts        = fStartCount,
            reverse_starts        = rStartCount,
            read_depth            = depth,
            template_depth        = tDepth,
            read_fraction_cut     = (fStartCount + rStartCount) / depth.toDouble,
            read_score            = score(pValue(fStartCount+rStartCount, depth, params.readStartProbability)),
            template_fraction_cut = (fStartCount + rStartCount) / tDepth.toDouble,
            template_score        = score(pValue(fStartCount+rStartCount, tDepth, params.templateStartProbability)),
            median_overhang       = medianOverhang,
            overhang_distribution = fs.mkString(";"),
            neighborhood_ratio    = nearby,
            aln_start             = alignment.map(_.start - 1).getOrElse(-1),
            aln_end               = alignment.map(_.end).getOrElse(-1),
            aln_strand            = alignment.map(_.strand.toString).getOrElse(""),
            aln_padded_guide      = alignment.map(_.paddedGuide).getOrElse(""),
            aln_alignment_string  = alignment.map(_.paddedAlignment).getOrElse(""),
            aln_padded_target     = alignment.map(_.paddedTarget).getOrElse(""),
            aln_mismatches        = alignment.map(_.mismatches).getOrElse(-1),
            aln_gap_bases         = alignment.map(_.gapBases).getOrElse(-1),
            aln_mm_and_gaps       = alignment.map(_.edits).getOrElse(-1)
          ))
        }

        // rev = Guide matches rev sequence, binds to fwd strand, yielding clean cut on fwd and variable cut on rev strand
        val rev = {
          val rPos           = pos - 1 + params.overhang
          val rRange         = Range.inclusive(rPos-params.maxOffset, rPos+params.maxOffset)
          val fStartCount    = fStarts(pos)
          val rs             = rRange.map(this.rStarts).reverse
          val rStartCount    = rs.sum
          val medianOverhang = median(rs.reverse, rStartCount) - params.maxOffset + params.overhang
          val bonusDepth     = rRange.iterator.filter(_ < pos).map(this.rStarts).sum
          val depth          = readDepth(pos) + bonusDepth
          val tDepth         = templateDepth(pos) + bonusDepth

          if (!ok(fStartCount, rStartCount, depth, params)) None else Some(CutSiteInfo(
            sample                = sample,
            guide                 = params.guide.getOrElse(""),
            enzyme                = params.enzyme.getOrElse(""),
            expected_overhang     = params.overhang,
            window_size           = params.maxOffset * 2 + 1,
            mapq_cutoff           = this.minMapQ,
            chrom                 = chrom,
            pos                   = pos - 1,
            strand                = "R",
            low_mapq_fraction     = lowMapqFraction(pos),
            forward_starts        = fStartCount,
            reverse_starts        = rStartCount,
            read_depth            = depth,
            template_depth        = tDepth,
            read_fraction_cut     = (fStartCount + rStartCount) / depth.toDouble,
            read_score            = score(pValue(fStartCount+rStartCount, depth, params.readStartProbability)),
            template_fraction_cut = (fStartCount + rStartCount) / tDepth.toDouble,
            template_score        = score(pValue(fStartCount+rStartCount, tDepth, params.templateStartProbability)),
            median_overhang       = medianOverhang,
            overhang_distribution = rs.mkString(";"),
            neighborhood_ratio    = nearby,
            aln_start             = alignment.map(_.start - 1).getOrElse(-1),
            aln_end               = alignment.map(_.end).getOrElse(-1),
            aln_strand            = alignment.map(_.strand.toString).getOrElse(""),
            aln_padded_guide      = alignment.map(_.paddedGuide).getOrElse(""),
            aln_alignment_string  = alignment.map(_.paddedAlignment).getOrElse(""),
            aln_padded_target     = alignment.map(_.paddedTarget).getOrElse(""),
            aln_mismatches        = alignment.map(_.mismatches).getOrElse(-1),
            aln_gap_bases         = alignment.map(_.gapBases).getOrElse(-1),
            aln_mm_and_gaps       = alignment.map(_.edits).getOrElse(-1)
          ))
        }

        (fwd, rev) match {
          case (None,    None   ) => None
          case (Some(f), None   ) => fwd
          case (None,    Some(r)) => rev
          case (Some(f), Some(r)) => if (r.template_score > f.template_score) rev else fwd
        }
      }
    }

    /** Turn a putative p-value into a score, accounting for the fact that the p-value may have underflowed. */
    private def score(pvalue: Double): Double = {
      require(pvalue >= 0 && pvalue <= 1, s"pvalue $pvalue is not in range 0-1.")
      if (pvalue == 0) - math.floor(-log10(Double.MinPositiveValue))
      else -log10(pvalue)
    }

    /** Given an array of counts, gives the index of the median. */
    private[digenome] def median(ns: Seq[Int], total: Int) = {
      // If there aren't any cut reads, arbitrarily return the middle index
      if (total == 0) {
        ns.length / 2
      }
      else {
        var sum      = 0
        var index    = -1
        var midpoint = total / 2
        while (sum < midpoint) {
          index += 1
          sum   += ns(index)
        }

        index
      }
    }

    /** Checks some basic filters and returns true if filters are passed, false if one or more filters fail. */
    private[digenome] def ok(fStarts:Int, rStarts:Int, depth: Int, filters: CutSiteParams): Boolean =
      fStarts >= filters.minForwardReads &&
        rStarts >= filters.minReverseReads &&
        (fStarts + rStarts) >= filters.minSupportingReads &&
        (fStarts + rStarts) / depth.toDouble >= filters.minSupportingFraction

    /** Calculates a p-value using a poisson approximation of a binomial process.*/
    private[digenome] def pValue(starts: Int, depth: Int, p: Double): Double = {
      if (depth == 0) 1 else {
        val dist = new PoissonDistribution(p*depth)
        dist.cumulativeProbability(starts-1) match {
          case 1 =>
            var sum = 0.0
            forloop(from=starts, until=depth+1) { i => sum += dist.probability(i) }
            sum
          case x =>
            1 - x
        }
      }
    }
  }


  /** Generate an iterator which filters out reads that are not part of a primary FR read pair. */
  private[digenome] def filtered(iter: Iterator[SamRecord], maxInsertSize: Int = 1200): Iterator[SamRecord] = {
    iter.filter(r => r.mapped && r.paired)
      .filterNot(r => r.duplicate || r.secondary || r.supplementary)
      .filter(r => r.isFrPair && math.abs(r.insertSize) <= maxInsertSize)
      .filter(r => math.abs(r.insertSize) >= CoordMath.getLength(r.start, r.end))
  }

  /** Samples out two million reads and calculates the mean read length and median insert size.
    * @param iterators iterators of SamRecords that can be consumed from and closed when used
    * @return a tuple of (read_length, insert_size)
    */
  private[digenome] def determineReadLengthAndInsertSize(iterators: Seq[Iterator[SamRecord] with Closeable]): (Int, Int) = {
    val sampleSize  = 2e6.toInt
    val isizes      = new NumericCounter[Int]
    val readLengths = new NumericCounter[Int]

    for (iter <- iterators; rec <- filtered(iter).take(sampleSize)) {
      isizes.count(math.abs(rec.insertSize))
      readLengths.count(rec.length)
    }

    iterators.foreach(_.safelyClose())
    (readLengths.mean().toInt, isizes.median().toInt)
  }

  /** Builds up an interval list to use for the whole genome by splitting on large blocks of Ns. */
  private[digenome] def wholeGenomeIntervalList(ref: PathToFasta, nCount: Int): IntervalList = {
    val walker = new ReferenceSequenceFileWalker(ref)
    val dict = walker.getSequenceDictionary
    val ns = new IntervalList(dict)

    val N= 'N'.toByte
    val n = 'n'.toByte

    dict.getSequences.foreach { chrom =>
      val name = chrom.getSequenceName
      val bases = walker.get(chrom.getSequenceIndex).getBases
      val length = bases.length

      var offset = 0
      while (offset < length) {
        if (bases(offset) == N || bases(offset) == n) {
          val start = offset + 1
          while (offset+1 < length && (bases(offset+1) == N || bases(offset+1) == n)) offset += 1
          val end = offset + 1
          if (end - start + 1 >= nCount) ns.add(new Interval(name, start, end))
        }

        offset += 1
      }
    }

    walker.safelyClose()
    IntervalList.invert(ns)
  }
}


@clp(group=ClpGroups.Digenome, description=
  """
    |Identifies cut sites from aligned whole digenome sequencing data.
    |
    |Ignores:
    |* Unmapped reads
    |* Reads flagged as duplicates
    |* Reads with mapping quality less than `min-map-q`
    |* Records flagged as secondary or supplementary
    |
  """)
class IdentifyCutSites
( @arg(flag='i', doc="One or more input BAM files of aligned digenome-sequencing data.") val input: Seq[PathToBam],
  @arg(flag='o', doc="Output file of predicted cut sites.") val output: FilePath,
  @arg(flag='g', doc="Guide sequence.") val guide: Option[String] = None,
  @arg(flag='p', doc="Optional PAM at the 5' end of the guide.") val pamFivePrime: Option[String] = None,
  @arg(flag='P', doc="Optional PAM at the 3' end of the guide.") val pamThreePrime: Option[String] = None,
  @arg(flag='e', doc="Name of the specific CAS9 enzyme.") val enzyme: Option[String] = None,
  @arg(flag='r', doc="Reference genome fasta file.") val ref: PathToFasta,
  @arg(flag='l', doc="Optional set of intervals to restrict analysis to.") val intervals: Option[PathToIntervals] = None,
  @arg(flag='m', doc="Discard reads below the given mapping quality.") val minMapQ: Int = 30,
  @arg(          doc="The expected overhang at the cut site. Positive for 5' overhang, negative for 3'.") val overhang: Int = 0,
  @arg(flag='x', doc="Counts reads starting this many bases from expected start sites.") val maxOffset: Int = 2,
  @arg(flag='d', doc="Don't check that reference genome and BAM file have the same set of sequences.") val dontCheckRef: Boolean = false,
  // Filtering arguments
  @arg(flag='N', doc="Minimum depth of coverage to output site.") val minDepth: Int = 10,
  @arg(flag='X', doc="Maximum depth of coverage to output site.") val maxDepth: Int = 300,
  @arg(flag='Q', doc="Maximum fraction of reads at site with mapq <= min-map-q.") val maxLowMapqFraction: Double = 0.3,
  @arg(flag='F', doc="Minimum F reads starting at the cut site.") val minForwardReads: Int = 4,
  @arg(flag='R', doc="Minimum R reads starting at the cut site.") val minReverseReads: Int = 4,
  @arg(flag='M', doc="Minimum total reads starting at the cut site.") val minSupportingReads: Int = 8,
  @arg(flag='S', doc="Minimum fraction of covering reads supporting the cut site") val minSupportingFraction: Double = 0.2
) extends EditasTool {
  import IdentifyCutSites._

  // Check inputs and output are readable/writable
  Io.assertReadable(input)
  Io.assertReadable(ref)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertReadable)

  // Some input validations
  validate(minForwardReads >=1 || minReverseReads >=1 || minSupportingReads >= 1,
    "One of --min-forward-reads, --min-reverse-reads or --min-supporting-reads must be set to >= 1.")

  validate(maxDepth >= minDepth, "max-depth must be >= min-depth")
  validate(minForwardReads    <= maxDepth, "min-forward-reads must be <= max-depth")
  validate(minReverseReads    <= maxDepth, "min-reverse-reads must be <= max-depth")
  validate(minSupportingReads <= maxDepth, "min-supporting-reads must be <= max-depth")

  validate(guide.isEmpty == enzyme.isEmpty, "Must provide both guide and enzyme or neither.")

  val validBases = "ATUGCSWRYKMBVHDN"
  guide.foreach         { s => validate(s.forall(ch => validBases.indexOf(ch) >= 0), s"Guide $s contained invalid bases. May only contain $validBases.") }
  pamFivePrime.foreach  { s => validate(s.forall(ch => validBases.indexOf(ch) >= 0), s"PAM $s contained invalid bases. May only contain $validBases.") }
  pamThreePrime.foreach { s => validate(s.forall(ch => validBases.indexOf(ch) >= 0), s"PAM $s contained invalid bases. May only contain $validBases.") }

  override def execute(): Unit = {
    val sources = input.map(bam => SamSource(bam))
    val sample  = sources.flatMap(bam => bam.header.getReadGroups.map(_.getSample)).distinct.sorted.mkString(",")

    val ref  = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.ref)
    val dict = ref.getSequenceDictionary

    sources.foreach { source =>
      val name = source.toSamReader.getResourceDescription
      require(source.header.getSortOrder == SortOrder.coordinate, s"Input BAM is not coordinate sorted: $name")
      require(source.indexed, s"Input BAM is not indexed: $name")
      if (!dontCheckRef) require(source.dict.isSameDictionary(dict), s"Input BAMs doesn't match provided reference: $name.")
    }

    logger.info(f"Sampling reads to determine read length and insert size.")
    val (readLength, insertSize) = determineReadLengthAndInsertSize(sources.map(_.iterator))
    logger.info(f"Determined read_length=$readLength and median_insert_size=$insertSize")

    // Either query over the provided regions or go over the whole file
    val regions = intervals match {
      case Some(path) => IntervalList.fromPath(path).uniqued().getIntervals.toIndexedSeq
      case None       => IdentifyCutSites.wholeGenomeIntervalList(this.ref, 2*insertSize).getIntervals.toIndexedSeq
    }

    val params = CutSiteParams(
      readLength            = readLength,
      insertSize            = insertSize,
      guide                 = this.guide,
      pamFivePrime          = this.pamFivePrime,
      pamThreePrime         = this.pamThreePrime,
      enzyme                = this.enzyme,
      overhang              = this.overhang,
      maxOffset             = this.maxOffset,
      minDepth              = this.minDepth,
      maxDepth              = this.maxDepth,
      maxLowMapqFraction    = this.maxLowMapqFraction,
      minForwardReads       = this.minForwardReads,
      minReverseReads       = this.minReverseReads,
      minSupportingReads    = this.minSupportingReads,
      minSupportingFraction = this.minSupportingFraction
    )

    val progress = ProgressLogger(logger, unit=2e6.toInt)
    val out = Metric.writer[CutSiteInfo](output)

    regions.foreach { region =>
      // Setup some tracking objects
      val chrom = region.getContig
      val start = region.getStart
      val end   = region.getEnd
      val accumulator = new CutSiteAccumulator(sample=sample, chrom=chrom, start=start, end=end, minMapQ=this.minMapQ, ref=ref)

      // Now iterate over the reads in the region
      for (bam <- sources; rec <- filtered(bam.query(chrom, start, end, QueryType.Overlapping))) {
        accumulator.accumulate(rec)
        progress.record(rec)
      }

      // Window size is intended to capture the maximum distance between two "calls" that are likely
      // generated by the same actual cut site.  E.g with an overhang of -3 we might make one call
      // on the F strand and another call on the R strand 4 bases away.
      val windowSize = math.abs(this.overhang) + 1

      // Buffer the cut sites as we generate them and then only flush them once we have all the ones
      // that are nearby each other
      val buffer = new mutable.ArrayBuffer[CutSiteInfo]

      def flushBuffer(): Unit = {
        if (buffer.size <= 1 || buffer.size > 10) buffer.foreach(out.write(_))
        else out.write(buffer.maxBy(x => (x.template_score, x.template_fraction_cut)))
        buffer.clear()
      }

      forloop(from=start, until=end+1) { pos =>
        accumulator.build(pos, params) match {
            case None => Unit
            case Some(best) =>
              if (buffer.nonEmpty && buffer.last.pos < best.pos-windowSize) flushBuffer()
              buffer += best
        }
      }

      flushBuffer()
      out.flush()
    }

    out.close()
  }
}
