// #############################################################################
// Validate reference and create cutsites file
//
// #############################################################################
import scala.io._
import java.io._


// load the command line parameters 
val inputRef = Source.fromFile(args(0)).getLines().drop(1).map{line => line}.mkString
val inputSites = Source.fromFile(args(1)).getLines().drop(1).map{line => line.split("\t")(0)}.toArray
val forwardPrimer = args(2)
val reversePrimer = args(3)
val primerExtend = args(4).toInt


val outputPrimersFile = new PrintWriter(args(0) + ".primers")
// validate that the foward and reverse primers are in the target reference
if (!(inputRef contains forwardPrimer))
  throw new IllegalStateException("Unable to find forward primer " + forwardPrimer + " in reference")
if (!(inputRef contains reversePrimer))
  throw new IllegalStateException("Unable to find reverse primer " + reversePrimer + " in reference")

outputPrimersFile.write(inputRef.slice(inputRef.indexOf(forwardPrimer), inputRef.indexOf(forwardPrimer) + forwardPrimer.length() + primerExtend) + "\n" + 
  inputRef.slice(inputRef.indexOf(reversePrimer) - primerExtend, inputRef.indexOf(reversePrimer) + reversePrimer.length()) + "\n")
outputPrimersFile.close()

// automatic extension files we generate
val outputFile = new PrintWriter(args(0) + ".cutSites")

def reverseComp(c: Char): Char = if (c == 'A') 'T' else if (c == 'C') 'G' else if (c == 'G') 'C' else if (c == 'T') 'A' else c
def reverseCompString(str: String): String = str.map{reverseComp(_)}.reverse.mkString

println("Running...")
outputFile.write("site\tposition\tcutPos\n")
inputSites.foreach{site => {
  val pos = inputRef.indexOf(site)
  if (pos >= 0) {
    println("FWD")
    outputFile.write(site + "\t" + (pos+1) + "\t" + (pos + 1 + (site.length - 7)) + "\n")
  } else {
    val posRev = reverseCompString(inputRef).indexOf(site)
    if (posRev >= 0) {
      println("REV")
      outputFile.write(site + "\t" + (inputRef.size - (posRev+site.size)) + "\t" + ((inputRef.size - (posRev + site.size)) + 7) + "\n")
    } else {
      println("No position found")
    }
  }

}}
outputFile.close()