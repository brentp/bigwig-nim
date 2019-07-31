import unittest
import bigwig

suite "test reading":
  test "that missing file returns false":
    var bw: BigWig
    check false == open(bw, "xxxxxx.bw")
