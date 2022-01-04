// convert a run where the UMI is read out as an indexed run to a more traditional run where
// we prepend it to the front of the first read
import java.io._
import scala.io._

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

def gos(s: String) = new DataOutputStream(new GZIPOutputStream(new FileOutputStream(s)))
def gis(s: String) = new DataInputStream(new GZIPInputStream(new FileInputStream(s)))

val inputReads = Source.fromInputStream(gis(args(0))).getLines().grouped(4)

val inputBarcode = Source.fromInputStream(gis(args(1))).getLines().grouped(4)

val inputUMI = Source.fromInputStream(gis(args(2))).getLines().grouped(4)

val outputReads = new PrintWriter(gos(args(3)))

inputReads.foreach{case(read) => {
    val barcode = inputBarcode.next()
    val umi = inputUMI.next()
    
    val barcodePlusUMI = barcode(1) + umi(1)
    val barcodePlusUMIQual = barcode(3) + umi(3)
    
    outputReads.write(read(0) + "\n" + barcodePlusUMI + read(1) + "\n" + read(2) + "\n" + barcodePlusUMIQual + read(3) + "\n")
}}
outputReads.flush()
outputReads.close()