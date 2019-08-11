import nimbigwig/bigWig
export bbIsBigBed, bwIsBigWig

import ./bigwigpkg/version
export version

type BigWig* = ref object
  bw : ptr bigWigFile_t
  path: string
  isBigBed: bool

  # these are internal, re-used cache before sending
  # to get data in format for bw from more common format
  starts: seq[uint32]
  stops: seq[uint32]
  values: seq[float32]
  cs: cstringArray

proc c_free(p: pointer) {.
  importc: "free", header: "<stdlib.h>".}

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

template isBigBed*(b:BigWig): bool =
  ## indicate wether file is BigBed (true) or BigWig (false)
  b.isBigBed

proc close*(bw: BigWig) =
  ## close the file and free up resources
  if bw.bw != nil:
    bwClose(bw.bw)
    bw.bw = nil
  if bw.cs != nil:
    deallocCStringArray(bw.cs)
    bw.cs = nil

proc destroy_bigwig(bw: BigWig) =
  bw.close

proc open*(bw: var BigWig, path: string, mode: FileMode=fmRead, maxZooms:int=8): bool =
  ## open the bigwig file. maxZooms is only used when opening in write mode.
  var fmode: string
  new(bw, destroy_bigwig)
  if mode == fmRead: fmode = "r" elif mode == fmWrite: fmode = "w" elif mode == fmAppend: fmode = "a"
  if mode == fmRead and bbIsBigBed(path, nil) == 1:
    bw = BigWig(bw: bbOpen(path, nil), path: path, isBigBed: true)
  else:
    bw = BigWig(bw: bwOpen(path, nil, fmode), path: path)
  if bw.bw == nil: return false
  result = true
  if mode  == fmWrite:
    result = 0 == bw.bw.bwCreateHdr(maxZooms.int32)
    bw.cs = allocCStringArray(@[""])

type CPtr[T] = ptr UncheckedArray[T]

proc SQL*(bw: BigWig): string =
  # return any SQL associated with a bigbed file; this can be used to parse the
  # extra columns in bigbed
  var cs = bbGetSQL(bw.bw)
  result = $cs
  c_free(cs)

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


iterator entries*(bw: var BigWig, chrom: string, start:int=0, stop:int= -1): tuple[start: int, stop: int, value: cstring] =
  ## yield bigbed entries. any values is returned as a string
  var stop = bw.get_stop(chrom, stop)
  var it = bbOverlappingEntriesIterator(bw.bw, chrom, start.uint32, stop.uint32, 0.cint, 1000)
  while it.data != nil:
    let starts = cast[CPtr[uint32]](it.entries.start)
    let stops = cast[CPtr[uint32]](it.entries.end)
    for i in 0..<it.entries.l:
      yield (starts[i].int, stops[i].int, if it.entries.str != nil: cast[cstringArray](it.entries.str)[i] else: nil)

    it = bwIteratorNext(it)
  it.bwIteratorDestroy


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
  deallocCStringArray(chroms)

proc writeHeader*(bw:BigWig) =
  ## write the header (which must have been added in `setHeader` to file.
  doAssert bw.bw.cl != nil, "[bigwig] attempted to call writeHeader on empty header; use setHeader first"
  doAssert 0 == bw.bw.bwWriteHdr

template setLens(bw:BigWig, L:int) =
  bw.starts.setLen(L)
  bw.stops.setLen(L)
  bw.values.setLen(L)

proc add*[T: int|uint32|uint64|int32|int64](bw:BigWig, chrom: string, intervals: seq[tuple[start:T, stop: T, value: float32]]) =
  ## add intervals to the bigwig.
  if intervals.len == 0: return
  bw.setLens(intervals.len)

  for i, iv in intervals:
    bw.starts[i] = iv.start.uint32
    bw.stops[i] = iv.stop.uint32
    bw.values[i] = iv.value

  bw.cs[0] = chrom.cstring
  doAssert 0 == bw.bw.bwAddIntervals(bw.cs, bw.starts[0].addr, bw.stops[0].addr, bw.values[0].addr, 1'u32), "[bigwig] error adding intervals"

  if intervals.len > 1:
    doAssert 0 == bw.bw.bwAppendIntervals(bw.starts[1].addr, bw.stops[1].addr, bw.values[1].addr, intervals.high.uint32), "[bigwig] error appending intervals"


proc add*[T: int|uint32|uint64|int32|int64, U: int|uint32|uint64|int32|int64](bw:BigWig, chrom: string, span: U, intervals: seq[tuple[start:T, value: float32]]) =
  ## add spans to the bigwig. this adds fixed-length (span) intervals starting at the given positions.
  if intervals.len == 0: return
  bw.setLens(intervals.len)

  for i, iv in intervals:
    bw.starts[i] = iv.start.uint32
    bw.values[i] = iv.value

  doAssert 0 == bw.bw.bwAddIntervalSpans(chrom.cstring, bw.starts[0].addr, span.uint32, bw.values[0].addr, 1'u32), "[bigwig] error adding interval spans"
  if intervals.len > 1:
    doAssert 0 == bw.bw.bwAppendIntervalSpans(bw.starts[1].addr, bw.values[1].addr, intervals.high.uint32), "[bigwig] error appending interval spans"

proc add*(bw:BigWig, chrom:string, start: uint32, values: var seq[float32], step:uint32=1, span:uint32=1) =
  ## add values to the bigwig starting at start and stepping by step.
  ## this is the most efficient way (space and performance) to add to a bigwig file if your intervals match this format.
  if values.len == 0: return
  doAssert 0 == bw.bw.bwAddIntervalSpanSteps(chrom, start, span, step, values[0].addr, 1)
  if values.len > 1:
    doAssert 0 == bw.bw.bwAppendIntervalSpanSteps(values[1].addr, values.high.uint32)
