package com.editasmedicine.digenome

import com.fulcrumgenomics.util.Metric

object CutSiteInfo {
  val ScientificFormatString = "0.0000E00"
  private val scientificFormat = new java.text.DecimalFormat(ScientificFormatString)
  private val binomial = new com.fulcrumgenomics.math.BinomialDistribution
}

case class CutSiteInfo(digenomitas_version: String = getClass.getPackage.getImplementationVersion,
                       sample: String,
                       guide: String,
                       enzyme: String,
                       expected_overhang: Int,
                       window_size: Int,
                       mapq_cutoff: Int,
                       chrom: String,
                       pos: Int,
                       strand: String,
                       low_mapq_fraction: Double,
                       forward_starts: Int,
                       reverse_starts: Int,
                       read_depth: Int,
                       read_fraction_cut: Double,
                       read_score: Double,
                       template_depth: Int,
                       template_fraction_cut: Double,
                       template_score: Double,
                       median_overhang: Int,
                       overhang_distribution: String,
                       neighborhood_ratio: Double,
                       aln_start: Int,
                       aln_end: Int,
                       aln_strand: String,
                       aln_padded_guide: String,
                       aln_alignment_string: String,
                       aln_padded_target: String,
                       aln_mismatches: Int,
                       aln_gap_bases: Int,
                       aln_mm_and_gaps: Int
                      ) extends Metric {

  override protected def formatValue(value: Any): String = value match {
    case d: Double if d == 0.0 => "0.0"
    case x   => super.formatValue(x)
  }
}

