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
