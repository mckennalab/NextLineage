import scala.io._
import java.io._
import scala.sys.process._


// inputs
val edaFile = args(0)
val alignerOpt = args(1)
val mergedReads = args(2)
val alignedPairs = args(3)
val ref = args(4)

// outputs
val outputMergedReads = args(5)
val outputAlignedPairs = args(6)

// gap open and extend
val gapOpen = args(7).toDouble
val gapExt = args(8).toDouble

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
      val aligner = "java -Xmx7g -jar /dartfs-hpc/rc/lab/M/McKennaLab/resources/tools/bin/HMMAlign-assembly-1.2.jar " + 
      "--inputReads " + fq + " --reference " + ref + " --outputReads " + outFasta + " --threads 24"
      println(aligner)
      val result = (aligner).!
      if (result == 0)
        return true
      return false
    }
      // MAFFT
    case "mafft" => {
      // this is a little more involved -- we have to split out each read from the input, align it, and add it to the output file

      val aligner = "/net/gs/vol1/home/aaronmck/tools/bin/mafft --maxiterate 1000 --op " + gapOpen + " --genafpair "
      val outputWriter = new PrintWriter(outFasta)


      Source.fromFile(fq).getLines().grouped(4).zipWithIndex.foreach{case(group,index) => {
        val merged = concatTwoTemp(ref,group(0).stripPrefix("@"),group(1))

        // make an array of SequenceReads to store the result in
        /**
          * Run MAFFT and capture the output
          */
        val out = new StringBuilder
        val err = new StringBuilder
        val logger = ProcessLogger(
          (o: String) => out.append(o + "\n"),
          (e: String) => err.append(e + "\n"))

        (aligner + merged) ! logger

        out.toString().split("\n") foreach { line => {
          //println(line)
          if (line.startsWith(">"))
            outputWriter.write(line + "\n")
          else
            outputWriter.write(line.toUpperCase + "\n")
        }}

        if (index % 1000 == 0)
          println("Processed " + index + " reads")
      }}
      outputWriter.close()

      return true
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
println(alignedPairs + " " + (new File(alignedPairs)).exists)
val pairedLength = Source.fromFile(alignedPairs).getLines().size
println("Merged length " + mergedLength)
println("Paired length " + pairedLength)

if (mergedLength > 0) {
  println("Aligning merged file to " + outputMergedReads)
  val mergedAlignmentSuccess = doAlignment(mergedReads, outputMergedReads)
  if (!mergedAlignmentSuccess) {
    throw new IllegalStateException("Unable to align " + mergedReads + " to " + outputMergedReads)
  }
} else {
  println("Touching merged file to " + outputAlignedPairs)
  val mergedAlignmentSuccess = touchFile(mergedReads, outputMergedReads)
  if (!mergedAlignmentSuccess) {
    throw new IllegalStateException("Unable to touch " + mergedReads + " to " + outputMergedReads)
  }
}

if (pairedLength > 0) {
  println("Aligning pair file to " + outputAlignedPairs)
  val pairedAlignmentSuccess = doAlignment(alignedPairs, outputAlignedPairs)
  if (!pairedAlignmentSuccess) {
    throw new IllegalStateException("Unable to align " + alignedPairs + " to " + outputAlignedPairs)
  }
} else {
  println("Touching pair file to " + outputAlignedPairs)
  val pairedAlignmentSuccess = touchFile(alignedPairs, outputAlignedPairs)
  if (!pairedAlignmentSuccess) {
    throw new IllegalStateException("Unable to touch " + alignedPairs + " to " + outputAlignedPairs)
  }
}
