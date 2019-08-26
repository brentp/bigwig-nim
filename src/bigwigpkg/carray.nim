from strformat import `&`

type CArray*[T] = object
  c: ptr UncheckedArray[T]
  L: int

proc `[]`*[T](A: CArray[T], i:SomeOrdinal): T {.inline.} =
  ## get the value at the i'th index of a CArray
  let i = int(if i < 0: A.L + i else: i)
  when compileOption("boundChecks"):
    if i >= A.len:
     raise newException(IndexError, &"index: {i} out of bounds in CArray of length {A.L}")
  return A.c[i]

proc `[]`*[T](A: CArray[T], b:BackwardsIndex): T {.inline.} =
  ## get the value at the b'th index from the end of a CArray
  return A[A.len - b.int]

proc `[]=`*[T](A: CArray[T], i:SomeOrdinal, value: T) {.inline.} =
  ## set the value at the i'th index of a CArray
  let i = int(if i < 0: A.L + i else: i)
  when compileOption("boundChecks"):
    if i >= A.len:
     raise newException(IndexError, &"index: {i} out of bounds in CArray of length {A.L}")
  A.c[i] = value

proc fix_slice*[U: SomeOrdinal, V: SomeOrdinal, T](A: CArray[T], sl: HSlice[U, V]): Slice[int] {.inline.} =
  result.a = int(if sl.a < 0: A.L + sl.a else: sl.a)
  result.b = int(if sl.b < 0: A.L + sl.b else: sl.b)

proc `[]`*[U: SomeOrdinal, V: SomeOrdinal, T](A: CArray[T], sl:HSlice[U, V]): CArray[T] {.inline.} =
  ## extract a new CArray for the given slice
  let sl = A.fix_slice(sl)
  when compileOption("boundChecks"):
    if sl.b >= A.len or sl.a >= A.len:
      raise newException(IndexError, &"index: {sl.b} out of bounds in CArray of length {A.L}")
  result.c = cast[ptr UncheckedArray[T]](A.c[sl.a].addr)
  result.L = sl.b - sl.a + 1

proc toSeq*[T](A: CArray[T]): seq[T] =
  result = newSeq[T](A.len)
  copyMem(result[0].addr, A.c, sizeof(T) * A.L)

proc len*[T](A: CArray[T]): int =
  ## get the length of a CArray
  return A.L

proc carray*[T](c: ptr UncheckedArray[T], L:int): CArray[T] {.inline.} =
  ## create a new CArray from an UncheckedArray.
  result.c = c
  result.L = L

proc carray*[T](c: ptr T, L:int): CArray[T] {.inline.} =
  ## create a new CArray from a pointer.
  result.c = cast[ptr UncheckedArray[T]](c)
  result.L = L

proc c_free(p: pointer) {.
  importc: "free", header: "<stdlib.h>".}

proc free*[T](A: var CArray[T]) =
  ## helper function to free the memory associated with the CArray
  ## this should only be called when the user is responsible for freeing the
  ## underlying memory.
  if A.c != nil:
    c_free(A.c)
    A.c = nil
  A.L = 0

proc `$`*[T](A: CArray[T]): string =
  return &"CArray[{$T}](len: {A.len})"

when isMainModule:
  import unittest
  import math
  var x = newSeq[float32](10)
  x[2] = 43.2

  proc `~~`[T: float32|float64](a: T, b: T): bool =
    return abs(a - b) < 1e-4

  suite "CArray":
    test "creatin and length work":
      let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)
      check v.len == x.len

    test "that carray indexing works":
      let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)
      check v[2] == x[2]

    test "that carray slicing":
      let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)
      check v[(2..<4)][0] == v[2]


    when compileOption("boundChecks"):
      test "that indexing out of bounds raises error":

        expect IndexError:
          let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)
          discard v[10]

    test "toSeq":
          let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)
          let s = v.toSeq
          for i, sv in s:
            check sv == x[i]

    test "setting and negative indexing":
        let v = carray(cast[ptr UncheckedArray[float32]](x[0].addr), x.len)

        v[3] = 44.2
        check x[3] ~~ 44.2
        check v[3] ~~ 44.2

        v[-2] = 98.1
        check x[^2] ~~ 98.1
        check v[^2] ~~ 98.1
        check v[^1] ~~ 0



