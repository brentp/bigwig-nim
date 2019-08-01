import unittest
import bigwig
import math

proc `~~`(a: float32, b: float32): bool =
  return abs(a - b) < 1e-9

suite "test reading":
  test "that missing file returns false":
    var bw: BigWig
    check false == open(bw, "xxxxxx.bw")

  test "that reading values works":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")

    var values: seq[float32]
    bw.values(values, "1")
    check values[0] ~~ 0.1
    check values[1] ~~ 0.2
    check values[2] ~~ 0.3

  test "that interval iteration works":
    var bw: BigWig
    check true == open(bw, "tests/test.bw")
    var expect : seq[tuple[start:int, stop:int, value:float32]] = @[
        (start: 0, stop: 1, value: 0.1000000014901161'f32),
        (start: 1, stop: 2, value: 0.2000000029802322'f32),
        (start: 2, stop: 3, value: 0.300000011920929'f32),
        (start: 100, stop: 150, value: 1.399999976158142'f32),
        (start: 150, stop: 151, value: 1.5'f32)]

    var i = 0
    for iv in bw.intervals("1"):
      check iv == expect[i]
      i += 1

