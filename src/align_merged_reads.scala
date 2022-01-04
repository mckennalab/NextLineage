import scala.io._
import java.io._
import scala.sys.process._


// inputs
val edaFile = args(0)
val alignerOpt = args(1)
val mergedReads = args(2)
val ref = args(3)

// outputs
val outputMergedReads = args(4)

// gap open and extend
val gapOpen = args(5).toDouble
val gapExt = args(6).toDouble
val paired = args(7).toBoolean

// concat two fasta files into a temp file
def concatTwoTemp(in1: String, readName: String, readString: String): String = {
  val output = File.createTempFile("mafft","alignment")
  //println(output)
  output.deleteOnExit()
  val outputFl = new PrintWriter(output.getAbsolutePath)
  Source.fromFile(in1).getLines.foreach{line => outputFl.write(line + "\n")}
  outputFl.write(">" + readName + "\n")
  outputFl.write(readString + "\n")
  outputFl.close
  output.getAbsolutePath
}

// do the alignments
def doAlignment(fq: String, outFasta: String): Boolean = {
  // check what aligner they'd like to use
  val aligner = alignerOpt.toLowerCase match {

    // NEEDLEALL
    case "needle" => {
      val aligner = "/dartfs/rc/lab/M/McKennaLab/resources/tools/bin/needleall -datafile " + edaFile + " -snucleotide1 -snucleotide2 -aformat3 fasta -gapextend " + gapExt + " -gapopen " + gapOpen + " -asequence " + ref + " -bsequence " + fq + " -outfile " + outFasta
      println(aligner)
      val result = (aligner).!
      if (result == 0)
        return true
      return false
    }

      // Align with our homebrew convex / general gap penality code
    case "convex" => {
      var aligner = "java -Xmx24g -jar /dartfs/rc/lab/M/McKennaLab/projects/Zon_lab_GESTALT/2021_12_05_GESTALT/2021_01_12_10X_data/HMMAlign-assembly-1.2.jar " + 
      "--inputReads " + fq + " --reference " + ref + " --outputReads " + outFasta + " --threads 20 --matchScore 3 --mismatchCost 2 --gapOpenCost 10"
      if (paired)
      aligner += " --paired TRUE"
      println(aligner)
      val result = (aligner).!
      if (result == 0)
        return true
      return false
    }
    case _ => throw new IllegalStateException("Unable to run aligner: " + alignerOpt + " as it's unknown")
  }
  return false
}

// touch the output fasta file
def touchFile(fq: String, outFasta: String): Boolean = {
  println("touch " + outFasta)
  val result = ("touch " + outFasta).!
  if (result == 0)
    return true
  return false
}

// check that there are reads in each of the input read files
println(mergedReads + " " + (new File(mergedReads)).exists)
val mergedLength = if ((new File(mergedReads)).exists) Source.fromFile(mergedReads).getLines().size else 0
println("Merged length " + mergedLength)

if (mergedLength > 0) {
  println("Aligning merged file to " + outputMergedReads)
  val mergedAlignmentSuccess = doAlignment(mergedReads, outputMergedReads)
  if (!mergedAlignmentSuccess) {
    throw new IllegalStateException("Unable to align " + mergedReads + " to " + outputMergedReads)
  }
} else {
  println("Touching merged file to " + outputMergedReads)
  val mergedAlignmentSuccess = touchFile(mergedReads, outputMergedReads)
  if (!mergedAlignmentSuccess) {
    throw new IllegalStateException("Unable to touch " + mergedReads + " to " + outputMergedReads)
  }
}
