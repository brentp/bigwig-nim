import nimbigwig/bigWig
import os

import ./bigwigpkg/version

type BigWig* = ref object
  bw : ptr bigWigFile_t
  path: string

type BigWigHeader* = seq[tuple[name: string, length: int, tid: uint32]]

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

proc open*(bw: var BigWig, path: string, mode: FileMode=fmRead, maxZooms:int=8): bool =
  ## open the bigwig file. maxZooms is only used when opening in write mode.
  var fmode: string
  new(bw, destroy_bigwig)
  if mode == fmRead: fmode = "r" elif mode == fmWrite: fmode = "w" elif mode == fmAppend: fmode = "a"
  bw = BigWig(bw: bwOpen(path, nil, fmode), path: path)
  if bw.bw == nil: return false
  result = true
  if mode  == fmWrite:
    result = 0 == bw.bw.bwCreateHdr(maxZooms.int32)

type CPtr[T] = ptr UncheckedArray[T]

proc get_stop(bw: var BigWig, chrom: string, stop:int): int {.inline.} =
  if stop >= 0 : return stop
  let tid = bw.bw.bwGetTid(chrom)
  if tid == uint32.high:
    raise newException(KeyError, "[bigwig] unknown chromosome:" & chrom)
  result = cast[CPtr[uint32]](bw.bw.cl.len)[tid].int

proc header*(bw: var BigWig): BigWigHeader =
  result = newSeq[tuple[name: string, length: int, tid:uint32]](bw.bw.cl.nKeys)
  var lens = cast[CPtr[uint32]](bw.bw.cl.len)
  var names = cast[cstringArray](bw.bw.cl.chrom)
  for i in 0..<bw.bw.cl.nKeys:
    result[i] = ($names[i], lens[i].int, i.uint32)

proc values*(bw: var BigWig, values: var seq[float32], chrom: string, start:int=0, stop:int= -1, includeNA:bool=true) =
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

proc c_free(p: pointer) {.
  importc: "free", header: "<stdlib.h>".}

proc stats*(bw: var BigWig, chrom: string, start:int=0, stop:int= -1, stat:Stat=Stat.mean, nBins=1): seq[float64] =
  var stop = bw.get_stop(chrom, stop)
  var st = bw.bw.bwStats(chrom, start.uint32, stop.uint32, nBins.uint32, stat.int.bwStatsType)
  result = newSeqUninitialized[float64](nBins)
  copyMem(result[0].addr, st, nBins * sizeof(float64))
  c_free(st)

proc setHeader*(bw:BigWig, header:BigWigHeader) =
  ## set the header of a bigwig file opened for writing
  var nimChroms = newSeq[string](header.len)
  var lens = newSeq[uint32](header.len)
  for i, c in header:
    nimChroms[i] = c.name
    lens[i] = c.length.uint32
  var chroms = allocCStringArray(nimChroms)
  bw.bw.cl = bwCreateChromList(chroms, lens[0].addr, header.len.int64)
  doAssert bw.bw.cl != nil, "[bigwig] error creating header"

proc writeHeader*(bw:BigWig, header:BigWigHeader) =
  ## write the header (which must have been added in `setHeader` to file.
  doAssert bw.bw.cl != nil, "[bigwig] attempted to call writeHeader on empty header; use setHeader first"
  doAssert 0 == bw.bw.bwWriteHdr


