import unittest
import bigwig
import math

proc `~~`[T: float32|float64](a: T, b: T): bool =
  return abs(a - b) < 1e-4

suite "test reading":
  test "that missing file returns false":
    var bw: BigWig
    check false == open(bw, "xxxxxx.bw")

  test "that reading values works":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")

    var values: seq[float32]
    bw.values(values, "1", 0, 50)
    check values[0] ~~ 0.1
    check values[1] ~~ 0.2
    check values[2] ~~ 0.3
    check $values[3] == "nan"

  test "that interval iteration works":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")
    var expect : seq[tuple[start:int, stop:int, value:float32]] = @[
        (start: 0, stop: 1, value: 0.1'f32),
        (start: 1, stop: 2, value: 0.2'f32),
        (start: 2, stop: 3, value: 0.3'f32),
        (start: 100, stop: 150, value: 1.4'f32),
        (start: 150, stop: 151, value: 1.5'f32)]

    var i = 0
    for iv in bw.intervals("1"):
      check iv == expect[i]
      i += 1


  test "that stats work":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")

    var mean = bw.stats("1", 0, 3)
    check mean.len == 1
    check 0.2 ~~ mean[0]

    var mins = bw.stats("1", 0, 4, Stat.min, 4)
    check mins[0] ~~ 0.1
    check mins[1] ~~ 0.2
    check mins[2] ~~ 0.3
    check $mins[3] == "nan"

  test "header":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")
    check bw.header == @[(name: "1", length: 195471971, tid: 0'u32), (name: "10", length: 130694993, tid: 1'u32)]
