#!/bin/sh
exec scala "$0" "$@"
!#


import scala.io._
import java.io._
import scala.collection.mutable._

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



def output_stats(input_file_gz: String,
            output_file: PrintWriter, 
                sample_name: String,
                guides: Array[String], 
                knownTargets: Array[String]) {

  val stats_file = Source.fromFile(input_file_gz).getLines() // fromInputStream(gis(input_file_gz)).getLines()
  val header = stats_file.next()

    val targetSites = header.split("\t").zipWithIndex.filter{case(k,index) => k contains "sequence" }.map{case(k,v) => v}.toArray
    val outcomeSites = header.split("\t").zipWithIndex.filter{case(k,index) => k contains "target" }.map{case(k,v) => v}.toArray

    assert(guides.size == targetSites.size, "guide length " + guides.size + " does not equal target size " + targetSites.size)
  var cnt = 0
  var match_count = 0
  var known_count = 0
  var left = 0
  var right = 0
  var both = 0

  stats_file.foreach{line => {
    val sp = line.split("\t")
    if ((line contains "MERGED") && (line contains "PASS")) {
      
      targetSites.map{k => sp(k).slice(0,23)}.toArray.zipWithIndex.foreach{case(target,index) => {
        val guide  = guides(index)
        val outcome = sp(outcomeSites(index))
        val edits = left_right_both_neither(target,guide)
        if (edits._1 == "RIGHT")
          right += 1
        if (edits._1 == "LEFT")
          left += 1
        if (edits._1 == "BOTH")
          both += 1

        val guide_target_dist = edit_distance_right_align(target,guide)
        output_file.write(sample_name + "\t" + guide + "\t" + target + "\t" + outcome + "\t" + edits._1 + "\tpos[" + edits._2.mkString(",") + "]\t" + edit_distance_right_align(guide,target) + "\n")
          
      }}
    }
  }}

}
val output_file = new PrintWriter(args(2))
val inputs = Array[String](args(0))
val sample = args(1)
val cutSites = Source.fromFile(args(3)).getLines().drop(1).map{case(ln) => ln.split("\t")(0)}.toArray

val knownTargets = Source.fromFile("/dartfs/rc/lab/M/McKennaLab/projects/base_editing_targets/2021_07_28_sequencing_twist_plus_Steven/known_targets.txt").getLines().map{l => l.slice(0,20)}.toArray



output_file.write("sample\tguide\ttarget\toutcome\teditType\teditPositions\tdistance\n")
inputs.foreach{input_fl => {
  val input_file = (new File(input_fl)).getName()
  println(input_file)
  val tag = sample
  val name = sample
  println("Running tag: " + tag + " and name " + name)
  output_stats(input_fl,output_file,name,cutSites,knownTargets)
}}
output_file.close()



