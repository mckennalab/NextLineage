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

println("Finding target sites...")
outputFile.write("site\tposition\tcutPos\n")
var fwd = 0
var rev = 0
var outputArray = Array[Tuple3[String,Int,Int]]()

inputSites.foreach{site => {
  val sitePattern = ( site.slice(4,24) + "|" + reverseCompString(site.slice(4,24))).r
  val positions = sitePattern.findAllMatchIn(inputRef).toArray
  
  println(site + " " + positions.size)
  positions.foreach{pos => {
    if (inputRef.slice(pos.start,pos.start+20) ==  site.slice(4,24)) {
      outputArray :+= (inputRef.slice(pos.start-4,pos.start+20),(pos.start - 4),(pos.start + 17))
      fwd += 1
    }else{
      outputArray :+= (inputRef.slice(pos.start-4,pos.start+20),(pos.start - 4),(pos.start - 1))
      rev += 1
    }
  }}
}}

scala.util.Sorting.stableSort(outputArray,(e1: (String, Int, Int), e2: (String, Int, Int)) => e1._2 < e2._2)
outputArray.foreach{case(site,start,cutsite) => {
  outputFile.write(site + "\t" + start + "\t" + cutsite + "\n")
}}

outputFile.close()
println("forward " + fwd)
println("rev " + rev)

