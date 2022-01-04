import scala.io._
import java.io._
import scala.collection.mutable._

def extract_non_insertion_sequence(read: String,reference: String,target_position: Int,length: Int): String = {
    var bases = ""
    var ref_pos = 0
    (0 until reference.size).foreach{char_position => {
        if ((ref_pos >= target_position) && (ref_pos < (target_position + length))) {
            if (reference(char_position) != '-') {
                bases = bases + read(char_position)
            }
        } 
        if (reference(char_position) != '-') {
            ref_pos += 1
        }
    }}
    bases
}

def smith_waterman(seq1: String, seq2: String, matchScore: Int = 1, mismatchScore: Int = -1): Int = {
  val mat = Array.ofDim[Int]((seq1.size + 1),(seq2.size + 1))
  (0 until seq1.size + 1).foreach{i => mat(0)(i) = i * mismatchScore}
  (0 until seq2.size + 1).foreach{j => mat(j)(0) = j * mismatchScore}
  
  (1 until seq1.size + 1).foreach{i => {
      (1 until seq2.size + 1).foreach{j => {
        val matchMismatchScore = if (seq1(i-1) == seq2(j-1)) matchScore else mismatchScore
        mat(i)(j) = math.max(matchMismatchScore + mat(i-1)(j-1),math.max(mat(i-1)(j) + mismatchScore,mat(i)(j-1) + mismatchScore))    
      }}
  }}
  mat(seq1.size)(seq2.size)
}

val known_closest = HashMap[String,Tuple2[String,Int]]()
def closest(guide: String, known_list: Array[String]): Tuple2[String,Int] = {
  if (known_closest contains guide) {
    return(known_closest(guide))
  } 

  var maxScore = Int.MinValue
  var minGuide: String = "----------------"
  known_list.foreach{known => {
    val score = smith_waterman(known,guide)
    if (score > maxScore) {
      maxScore = score
      minGuide = known
    }
  }}
  known_closest(guide) = (minGuide,maxScore)
  known_closest(guide)
}


def base_to_index(base: Char): Int = {
    base match {
       case x if base == 'A' || base == 'a' => 0 
       case x if base == 'C' || base == 'c' => 1
       case x if base == 'G' || base == 'g' => 2
       case x if base == 'T' || base == 't' => 3
        case _ => 4
    }
    
}
base_to_index('c')

def comp(base: Char): Char = base match {
  case 'A' => 'T'
  case 'T' => 'A'
  case 'C' => 'G'
  case 'G' => 'C'
  case _ => 'N'
}
def revcomp(str: String): String = str.map{case(s) => comp(s)}.reverse
def miror_complete(guide: String): String = revcomp(guide.slice(guide.size - 10,guide.size)) + guide.slice(guide.size - 10,guide.size)

def edit_distance_right_align(str1: String, str2: String): Int = {
    if (str1.size > str2.size) 
        str1.slice(str1.size - str2.size,str1.size).zip(str2).map{case(a,b) => if (a == b) 0 else 1}.sum
    else
        str2.slice(str2.size - str1.size,str2.size).zip(str1).map{case(a,b) => if (a == b) 0 else 1}.sum
}
    
def left_right_both_neither(read: String,ref: String): (String, Array[Int]) = {
    // forward strand convert A to G, 
    // reverse strand convert T to C
    var edited_positions = Array[Int]()
    var pos = 0
    var left = false
    var right = false
    read.slice(0,10).zip(ref.slice(0,10)).foreach{case(read_base,ref_base) => {
        if (ref_base == 'A' && read_base == 'G') {
            edited_positions :+= pos
            left = true
        }
        pos += 1
    }}
    read.slice(10,20).zip(ref.slice(10,20)).foreach{case(read_base,ref_base) => {
        if (ref_base == 'T' && read_base == 'C') {
            edited_positions :+= pos
            right = true
        }
        pos += 1
    }}
    val outcome = (left,right) match {
        case (x,y) if left && right => "BOTH"
        case (x,y) if left => "LEFT"
        case (x,y) if right => "RIGHT"
        case _ => "NEITHER"
    }
    return((outcome,edited_positions))
}

import java.util.zip._

// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))


def extract_guide_target_pair(read: String, reference: String): Tuple2[String,String] = {
  val front_half_read = read.slice(50,190)
  val front_half_ref  = reference.slice(50,190)
  val back_half_read  = read.slice(190,read.size - 50)
  val back_half_ref   = reference.slice(190,reference.size - 50)

  val guide_start = front_half_ref.indexOf('N')
  val guide_stop  = front_half_ref.lastIndexOf('N') + 1
  val target_start = back_half_ref.indexOf('N')
  val target_stop  = back_half_ref.lastIndexOf('N') + 1

  val guide = front_half_read.slice(guide_start,guide_stop)
  val target = back_half_read.slice(target_start,target_stop)

  (guide,target)
}


def output_stats(input_file_gz: String,output_file: PrintWriter,sample_name: String,guide: String, knownTargets: Array[String]) {
  val stats_file = Source.fromFile(input_file_gz).getLines() // fromInputStream(gis(input_file_gz)).getLines()
  val header = stats_file.next()

  var cnt = 0
  var match_count = 0
  var known_count = 0
  var left = 0
  var right = 0
  var both = 0

  stats_file.foreach{line => {
    val sp = line.split("\t")
    if ((line contains "MERGED") && (line contains "PASS")) {
      val count = sp(0).split("_")(2)
      val guide_target = extract_guide_target_pair(sp(18),sp(21))
      val extend_guide = miror_complete(guide_target._1)
      val closest_guide = closest(extend_guide, knownTargets)
      val closest_target = closest(guide_target._2.slice(3,23), knownTargets)
      val edits = left_right_both_neither(guide_target._2.slice(3,26),extend_guide)

      if (edits._1 == "RIGHT")
        right += 1
      if (edits._1 == "LEFT")
        left += 1
      if (edits._1 == "BOTH")
        both += 1

      val guide_target_dist = smith_waterman(extend_guide,guide_target._2.slice(3,23))
      output_file.write(sample_name + "\t" + count + "\t" + guide_target._1 + "\t" + closest_guide._1 + "\t" + closest_guide._2 + "\t" + closest_target._1 + "\t" + 
        closest_target._2 + "\t" + (closest_guide._1 == closest_target._1) + "\t" + guide_target_dist + "\t" + extend_guide + "\t" + guide_target._2 + "\t" + sp(22) + "\t" + sp(1) + "\t" + sp(6) + "\t" + edits._1 + 
        "\t" + edits._2.size + "\t" + edits._2.mkString(",") + "\n")
    }
  }}

}
val output_file = new PrintWriter(args(2))
val inputs = Array[String](args(0))
val sample = args(1)
val cutSites = Source.fromFile(args(3)).getLines().drop(1).next().split("\t")

val knownTargets = Source.fromFile("/dartfs/rc/lab/M/McKennaLab/projects/base_editing_targets/2021_07_28_sequencing_twist_plus_Steven/known_targets.txt").getLines().map{l => l.slice(0,20)}.toArray



output_file.write("sample\tcount\tguide\tclosestGuide\tclosestGuideScore\tclosestTarget\tclosestTargetScore\tclosestMatch\tgtDist\tpalguide\ttarget\toutcome\tpass\ttype\tediting\tnumberEdits\tpositions\n")
inputs.foreach{input_fl => {
  val input_file = (new File(input_fl)).getName()
  println(input_file)
  val tag = sample
  val name = input_file.split("\\.")(0)
  println("Running tag: " + tag + " and name " + name)
  output_stats(input_fl,output_file,name,cutSites(0),knownTargets)
}}
output_file.close()



