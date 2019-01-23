package com.editasmedicine.aligner

import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.testing.ReferenceSetBuilder
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.SequenceUtil

class AlignerTest extends UnitSpec {
  private val ref = {
    val builder = new ReferenceSetBuilder()
    builder.add("chr1") // 25 lines x 100 bases each = 2500b
      .add("AAATAGACCTTTCCCATTTATAACTTATTTGTAAAATGATTTCTATTATAAACATAACATATACATTGTATAACAATTAGAAAACCTGTCTGTTTTGATG") // 1-100
      .add("GATCTCAAGATTTAAGAAGGCTTAGACTTCAGCTATAAGATGCACATGCCACTGTGGGAGGCCGAGGCGGGCAGATCACGAGGTCAGGAGTTCTAGACCA") // 101-200
      .add("GCCTGACCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCATGGCAGCAGACACCTGTAATCCCAGTTATTCGGGAGGCTGA") // 201-300
      .add("GGCAGGAGAATTGCTTGAATGCAGGAGGCAGAGGTTGCAGTGAGCCGAGACGGCGCCACTGCACTCCAGCCTGGGCAACAGAGCAGATGGAGACCATCCT") // 301-400
      .add("GACCAACATGATGAAACTCTGTCTCTACTAAAAATACAAAAATTAGCTGGGCATGGTGGCGTGCACCTACTAGTCCCAGCTACTCGGGAGGCTGAGGCAG") // 401-500
      .add("GAGAATTGCTTGAACCCAGGAGGCGGAGGTTTCAGTGAGCCGATACCGCGCCATTGCACTCCAGCCTGGGCAACAGAGCGAGACTGTGTCTCAAAAAAAA") // 501-600
      .add("AAAAAAAAAGGAGATGCACATGTTTAAGTCTATTTCAGGCGGTTAGCTGGTGGATTGCTACAATTCCTCTGTAAGTTTAAAAAATCATGTAAGTGCTGTT") // 601-700
      .add("TTGGAGTACTGTAATAACTCTTGAGATGTAGAACACATCTGCAAAATGAGGGTAGTATAAAAGAGACGAGGGGATGAGGGTAATACATAAGAAATAGGGG") // 701-800
      .add("AAAGGACAAGAACAGGTAAATTAAACTTCAAGTACTATTTTTGCTATTGCTGTCTACACTCAACTAGCAAGGAAAAAGCCTTGCTTCTGCTCTGCGGGTT") // 801-900
      .add("TTCTTCGGGTTTAACTTGACCAAGCAAAACAGACCATCTGGGATTAACTTTTTCCTTTTCACTGTAGGTCACAGGCTCTACGTGTAGGGTGTTGGCCACC") // 901-1000
      .add("TGTTCTTCCACCATCTCTACCTCCACCTCCTCCTTTGTGGCCACAGCAATGTCACAGCCCATACATGGGGGAGGGGAGCATTCAGGAACTCGGAGGCAGA") // 1001-1100
      .add("TGCATTTTTTTCCAAACACAATAACCTCAAACAGTGGTCTCTAAGCACTTTCCTATGCTCTTCCAAAACGTGACCTCCCCTCTTACTCACACATCCCCTA") // 1101-1200
      .add("CACACGGAAAAGGACCACTATCCGTCCAGCCTGCGCTCGAGGGAGAAGTTTATACCTTCGTCCTAGAGATGCCAAATGCAGCAGGGAAGGCTGGACCGAG") // 1201-1300
      .add("GCAGCCGAGTGCTGGAAAGGGAGGCAAGAGGTGCGGGAGCGGGGAGAGGGGGAGGGGAGGCCGGGGCGCCGCGGGAGTAACCTCCACCGCACCCCACCGC") // 1301-1400
      .add("TCCGAGGGGCAGCCGGCCCGGCCCGAGTTTCTCCCCAGAAGCCTCCAGCCGCGGCTCTCGGGGAGGAGGAAGGAAGGGGTTCCCCGTCCAGGAAGCAGCA") // 1401-1500
      .add("CCAGCGGCGACCGCCTCCAGCCTCACCCTCCTCAGCCCCGCACCGCCCATTCCTCACTCCCCGCGCCGCCGCGTCCGCGCGCCTCCCCCCTGCAGACCCC") // 1501-1600
      .add("TCTCACCCAGCCCGCCCCGACCCCGCGCCCGCGCCCCCCACCCGCCCCTCCGGGGACCCCTAATTCATTCACTCGCCGCCGGCCCCGCCCGGCGCCGGCA") // 1601-1700
      .add("AAGAGGGTCGGGACCCGGGCAGGGGCCCAGGAGGGGTGGTCCGCTCCGTACCTCTCTCCCGCACCTGGGAGCCGCTGAGCCTCTGGCCCCGCCGCCGCCT") // 1701-1800
      .add("TCAGTGCCTGCGCCGCGCTCGCTCCCAGTCCGAAATGGCGGGGGCCGGGAGTACTGGCCGAGCCGCCGCCACCTTCGCCGCCGCCACTGCCGCCGCCGCT") // 1801-1900
      .add("GCTGCCTCCGCCGCCGCGGCCGCCGCCTAGGAAAATCGAGCTCCGAGCACACCGATGAGTTCGGGGCCGGGCGGCCGCAGAGGGCAGAGCTATCGATGCG") // 1901-2000
      .add("TTCCGCGCTCGATTCTTCTTCAGACGGGCGTACGAGAGGGAGCGGCTGAGGGCGGTGTGGGAAGAGGGAAGAGGGGGAGGCAGCGAGCGCCGGCGGGGAG") // 2001-2100
      .add("AAGGAGGGGGCCGGGCCGGGCCGGCGGGGGAGGAGCGGGGGCCGGGCCGGCGGAGGAAGGGGTGGCTGGGGCGGTCTAGGGTGGCGAGCCGGGCCGGCTG") // 2101-2200
      .add("GAGAGCGGGTCTGGGCGGCGCCTTGGCGGGAGGAGGGACTGCCGGACCCACGCGGCGGCCCGCCCCCTGCCTAGCCGCAAGGCTGTCCCCGCAGCCGCCA") // 2201-2300
      .add("ATTCTGACCCGGAGCGGGACCGGACCGCGGCGGGCTGTGCGGATGCCACCAGGGAGACGCCGCGAGCGGCCACGCCGCCCCGCTGACCGGTCTCCACAGA") // 2301-2400
    builder.add("chr2")
      .add("GATACaaCTCGTACTGTCAGT")
      .add("GATACGTCTCGTACTGTCAtT")
    val path = builder.toTempFile()
    ReferenceSequenceFileFactory.getReferenceSequenceFile(path)
  }

  val aligner = new Aligner(ref)

  "Aligner" should "align a perfect sequence without gaps or mismatches in the right spot on the F strand" in {
    val query  = ref.getSubsequenceAt("chr1", 50, 69).getBaseString
    val result = aligner.align(query, "chr1", 65)

    result.chrom  shouldBe "chr1"
    result.start  shouldBe 50
    result.end    shouldBe 69
    result.strand shouldBe 'F'
    result.paddedGuide shouldBe result.paddedTarget
    result.paddedAlignment.forall(_ == '|') shouldBe true
    result.score should be >= 0
  }

  it should "handle Us the same as Ts" in {
    val tQuery  = ref.getSubsequenceAt("chr1", 50, 69).getBaseString
    val uQuery  = tQuery.replace('T', 'U')
    uQuery should not be tQuery

    val tResult = aligner.align(tQuery, "chr1", 65)
    val uResult = aligner.align(uQuery, "chr1", 65)
    uResult.score shouldBe tResult.score
    uResult.paddedAlignment shouldBe tResult.paddedAlignment
  }

  it should "align a perfect sequence without gaps or mismatches in the right spot on the R strand" in {
    val aligner = new Aligner(ref)
    val query   = SequenceUtil.reverseComplement(ref.getSubsequenceAt("chr1", 50, 69).getBaseString)
    val result  = aligner.align(query, "chr1", 65)

    result.chrom  shouldBe "chr1"
    result.start  shouldBe 50
    result.end    shouldBe 69
    result.strand shouldBe 'R'
    result.paddedAlignment.forall(_ == '|') shouldBe true
    result.score should be >= 0
  }

  it should "align with a mismatch on the F strand" in {
    val query  = "GAGAATTGtTTGAACCCAGGnGG" // start of 5th line == 501-523
    val aligns = "||||||||.||||||||||||||"
    val result = aligner.align(query.toUpperCase, "chr1", 515)

    result.chrom shouldBe "chr1"
    result.start shouldBe 501
    result.end   shouldBe 523
    result.strand shouldBe 'F'
    result.paddedAlignment shouldBe aligns
    result.mismatches shouldBe 1
  }

  it should "handle ambiuity codes correctly in the PAM" in {
    // ref       "TCAGTGCCTGCGCCGCGCTCGCTCCCAGTCCGAAA" // start of 18th line == 1801-1835523
    val query  = "TCAGTGCCTGCGCCGCGCTCGCTCCCNRVCWSHDB" // start of 5th line == 501-523
    val aligns = "||||||||||||||||||||||||||||.|.|||."
    val result = aligner.align(query, "chr1", 1820)

    result.chrom shouldBe "chr1"
    result.start shouldBe 1801
    result.end   shouldBe 1835
    result.strand shouldBe 'F'
    result.paddedAlignment shouldBe aligns
    result.mismatches shouldBe 3
  }

  it should "align with a two bulges on the R strand" in {
    val query  = "AGGCTGG-GGCGGTCGCtCGCNGG" // revcomp of start of 16th line == 1501-1523
    val aligns = "||||||| ||||||||| ||||||"
    val result = aligner.align(query.filter(_.isLetter).toUpperCase, "chr1", 1510)

    result.chrom shouldBe "chr1"
    result.start shouldBe 1501
    result.end   shouldBe 1523
    result.strand shouldBe 'R'
    result.paddedAlignment shouldBe aligns
  }

  it should "prefer an alignment with two mismatches in the guide vs. one mismatch in the PAM" in {
    val query = "GATACGTCTCGTACTGTnrg"
    val result = aligner.align(query, "chr2", 22)
    result.chrom shouldBe "chr2"
    result.start shouldBe 1
    result.end shouldBe 20
    result.gapBases shouldBe 0
    result.mismatches shouldBe 2
  }
}
