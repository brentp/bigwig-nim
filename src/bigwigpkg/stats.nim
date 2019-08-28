import math

type StatFun* = proc (a: float32): float64

#https://stackoverflow.com/a/19045659
template med3[T](a, b, c:T): T = max(min(a,b), min(max(a,b),c))

proc qsort[T](a: var openarray[T], inl = 0, inr = -1) =
  # https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Nim
  # this is 3X faster than the algorithm.sort distributed with nim.
  var r = if likely(inr >= 0): inr else: a.high
  var l = inl
  let n = r - l + 1
  if n < 2: return
  let p = med3(a[l], a[r], a[l + n div 2])
  while l <= r:
    if a[l] < p:
      inc l
      continue
    if a[r] > p:
      dec r
      continue
    if l <= r:
      swap a[l], a[r]
      inc l
      dec r
  qsort(a, inl, r)
  qsort(a, l, inr)

proc mean*(a: openarray[float32]): float64 =
  for v in a:
    result += v
  result /= a.len.float64

proc fmin*(a: openarray[float32]): float64 =
  return min(a).float64

proc fmax*(a: openarray[float32]): float64 =
  return max(a).float64

proc ipercentiles*(a: var seq[float32], pcts: seq[float32]): seq[float64] =
  ## calculate percentiles on values that are guaranteed to be integers.
  var counts = newSeq[uint32](65536)
  let L = a.high.float64
  for v in a:
    counts[min(v.int, counts.high)].inc

  result = newSeqUninitialized[float64](pcts.len)
  for i, p in pcts:
    doAssert 0 <= p and p <= 1, "[error] pcts values to `percentiles` should be between 0 and 1."

    var S = 0'f64
    var j = 0
    var stop = p * L
    while S <= stop:
      S += counts[j].float64
      j.inc

    result[i] = (j - 1).float64

proc percentiles*(a: var seq[float32], pcts: seq[float32]): seq[float64] =
  #var allints = true
  #for v in a:
  #  if float(int(v)) != v:
  #    allints = false
  #    break
  #if allints:
  #  return ipercentiles(a, pcts)
  qsort(a)
  let L = a.high.float64
  result = newSeqUninitialized[float64](pcts.len)
  for i, p in pcts:
    doAssert 0 <= p and p <= 1, "[error] pcts values to `percentiles` should be between 0 and 1."
    result[i] = a[int(0.5 + p * L)]

proc std*(a: openarray[float32]): tuple[mean:float64, std: float64] =
  # s0 = sum(1 for x in samples)
  # s1 = sum(x for x in samples)
  # s2 = sum(x*x for x in samples)
  let N: float64 = a.len.float64
  if N == 0: return (0'f64, 0'f64)
  var S: float64
  var V: float64

  for v in a:
    S += v
    V += (v * v)
  result = (S/N, sqrt((N * V - S * S)/(N * (N - 1))))
