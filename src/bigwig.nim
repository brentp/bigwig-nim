import nimbigwig/bigWig
import os

import ./bigwigpkg/version

type BigWig* = ref object
  bw : ptr bigWigFile_t
  path: string


type Stat* {.pure.} = enum
    #doesNotExist = -1 #!< This does nothing */
    mean = 0
    stdev = 1
    max = 2
    min = 3
    # The number of bases covered
    coverage = 4
    # The sum of per-base values */
    sum = 5

proc close*(bw: BigWig) =
  ## close the file and free up resources
  if bw.bw != nil:
    bwClose(bw.bw)
    bw.bw = nil

proc destroy_bigwig(bw: BigWig) =
  bw.close

proc open*(bw: var BigWig, path: string, mode: FileMode=fmRead): bool =
  ## open the bigwig file
  var fmode: string
  new(bw, destroy_bigwig)
  if mode == fmRead: fmode = "r" elif mode == fmWrite: fmode = "w" elif mode == fmAppend: fmode = "a"
  bw = BigWig(bw: bwOpen(path, nil, fmode), path: path)
  return bw.bw != nil

type CPtr[T] = ptr UncheckedArray[T]

proc get_stop(bw: var BigWig, chrom: string, stop:int): int {.inline.} =
  if stop >= 0 : return stop
  let tid = bw.bw.bwGetTid(chrom)
  if tid == uint32.high:
    raise newException(KeyError, "[bigwig] unknown chromosome:" & chrom)
  result = cast[CPtr[uint32]](bw.bw.cl.len)[tid].int

proc values*(bw: var BigWig, values: var seq[float32], chrom: string, start:int=0, stop:int= -1, includeNA:bool=false) =
  ## exctract values for the given range into `values`
  var stop = bw.get_stop(chrom, stop)
  var ivs = bwGetValues(bw.bw, chrom, start.uint32, stop.uint32, includeNA.cint)
  values.setLen(ivs.l)
  copyMem(values[0].addr, ivs.value, sizeof(values[0]) * ivs.l.int)
  ivs.bwDestroyOverlappingIntervals()

iterator intervals*(bw: var BigWig, chrom: string, start:int=0, stop:int= -1): tuple[start: int, stop: int, value: float32] =
  ## iterate over the values in the given region
  var stop = bw.get_stop(chrom, stop)
  var it = bwOverlappingIntervalsIterator(bw.bw, chrom, start.uint32, stop.uint32, 1000)

  while it.data != nil:
    if it.entries != nil:
      raise newException(ValueError, "iterator not implemented for bigbed")

    elif it.intervals != nil:
      let starts = cast[CPtr[uint32]](it.intervals.start)
      let stops = cast[CPtr[uint32]](it.intervals.end)
      let values = cast[CPtr[cfloat]](it.intervals.value)
      for i in 0..<it.intervals.l:
        yield (starts[i].int, stops[i].int, values[i])

    it = bwIteratorNext(it)

  it.bwIteratorDestroy
