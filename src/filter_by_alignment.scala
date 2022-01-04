import scala.io._
import java.io._
import scala.sys.process._
import scala.collection.mutable._

// convert a run where the UMI is read out as an indexed run to a more traditional run where
// we prepend it to the front of the first read
import java.io._
import scala.io._

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

def gos(s: String) = new DataOutputStream(new GZIPOutputStream(new FileOutputStream(s)))
def gis(s: String) = new DataInputStream(new GZIPInputStream(new FileInputStream(s)))


// inputs
val reference = args(0)
val targetContigName = args(1)
val input_read1 = args(2)
val input_read2 = args(3)
val input_read3 = args(4)

val outputReads1 = args(5)
val outputReads2 = args(6)
val outputReads3 = args(7)

val tempReads = "temp.sam"

val aligner = Seq("bwa","mem",reference,input_read1)
println(aligner)

// align it
val result = (aligner #>> new File(tempReads)).!
if (result != 0)
  throw new IllegalStateException("Unable to align " + tempReads)

// process the reads, making a list of reads that align to the reference
val alignedReads = new HashMap[String, Boolean]()

Source.fromFile(tempReads).getLines().foreach{line => {
  if (!(line startsWith "@")) {
    val sp = line.split("\t")

    if (sp(2) == targetContigName) {
      //println(sp(0))
      alignedReads(sp(0)) = true
    }
  }
}}

val inputReads1 = Source.fromInputStream(gis(input_read1)).getLines().grouped(4)
val inputReads2 = Source.fromInputStream(gis(input_read2)).getLines().grouped(4)
val inputReads3 = Source.fromInputStream(gis(input_read3)).getLines().grouped(4)

val outReads1 = new PrintWriter(gos(outputReads1))
val outReads2 = new PrintWriter(gos(outputReads2))
val outReads3 = new PrintWriter(gos(outputReads3))

inputReads1.foreach{read1 => {
  val read2 = inputReads2.next()
  val read3 = inputReads3.next()

  if (alignedReads contains read1(0).stripPrefix("@").split(" ")(0)) {
    outReads1.write(read1(0) + "\n" + read1(1) + "\n" + read1(2) + "\n" + read1(3) + "\n")
    outReads2.write(read2(0) + "\n" + read2(1) + "\n" + read2(2) + "\n" + read2(3) + "\n")
    outReads3.write(read3(0) + "\n" + read3(1) + "\n" + read3(2) + "\n" + read3(3) + "\n")
  }
}}

outReads1.close()
outReads2.close()
outReads3.close()
