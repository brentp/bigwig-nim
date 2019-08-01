import unittest
import bigwig
import math

proc `~~`[T: float32|float64](a: T, b: T): bool =
  return abs(a - b) < 1e-4

suite "test writing":
  test "that setting header works":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")

    var hdr = bw.header
    bw.close

    var wtr:BigWig
    check true == open(wtr, "tests/writer.bw", fmWrite)
    wtr.setHeader(hdr)
    check wtr.header == hdr

    wtr.setHeader(hdr[0..<1])
    check wtr.header == hdr[0..<1]

  test "that adding intervals works":
    var wtr:BigWig
    check true == open(wtr, "tests/writer.bw", fmWrite)
    wtr.setHeader(@[(name:"chr1", length: 2000, tid: 0'u32)])
    wtr.writeHeader

    wtr.add("chr1", @[(start: 22, stop: 33, value: 0.01'f32), (start: 44, stop: 55, value: 155'f32)])
    wtr.close

    var rdr: BigWig
    check true == open(rdr, "tests/writer.bw")
    var i = 0
    for iv in rdr.intervals("chr1"):
      if i == 0: check iv == (start: 22, stop: 33, value: 0.01'f32)
      if i == 1: check iv == (start: 44, stop: 55, value: 155'f32)
      i.inc
    rdr.close

  test "that adding spans works":
    var wtr:BigWig
    check true == open(wtr, "tests/writer.bw", fmWrite)
    wtr.setHeader(@[(name:"chr1", length: 2000, tid: 0'u32)])
    wtr.writeHeader

    # add intervals with span of 100
    wtr.add("chr1", 100, @[(start: 22, value: 0.01'f32), (start: 44, value: 155'f32)])
    wtr.close

    var rdr: BigWig
    check true == open(rdr, "tests/writer.bw")
    var i = 0
    for iv in rdr.intervals("chr1"):
      if i == 0: check iv == (start: 22, stop: 122, value: 0.01'f32)
      if i == 1: check iv == (start: 44, stop: 144, value: 155'f32)
      i.inc
    rdr.close


  test "that add span step works":

    var values = @[1'f32, 2222.2'f32, 555.5'f32, 666.6'f32]

    var wtr:BigWig
    check true == open(wtr, "tests/writer.bw", fmWrite)
    wtr.setHeader(@[(name:"chr1", length: 2000, tid: 0'u32)])
    wtr.writeHeader

    # add 1-base intervals starting at 100
    wtr.add("chr1", 100, values)
    wtr.close

    var rdr: BigWig
    check true == open(rdr, "tests/writer.bw")
    var i = 0
    for iv in rdr.intervals("chr1"):
      if i == 0: check iv == (start: 100, stop: 101, value: 1'f32)
      if i == 1: check iv == (start: 101, stop: 102, value: 2222.2'f32)
      check iv.value == values[i]
      check iv.start == 100 + i
      check iv.stop == 100 + i + 1
      i.inc

  test "that add span step works with span and step":

    var values = @[1'f32, 2222.2'f32, 555.5'f32, 666.6'f32]

    var wtr:BigWig
    check true == open(wtr, "tests/writer.bw", fmWrite)
    wtr.setHeader(@[(name:"chr1", length: 2000, tid: 0'u32)])
    wtr.writeHeader
    let span = 200'u32
    let step = 33'u32

    # add intervals with start of 100
    wtr.add("chr1", 100, values, span=span, step=step)
    wtr.close

    var rdr: BigWig
    check true == open(rdr, "tests/writer.bw")
    var i = 0
    for iv in rdr.intervals("chr1"):
      check iv.start == 100 + i * step.int
      check iv.stop == 100 + i * step.int + span.int
      check iv.value == values[i]
      i.inc
