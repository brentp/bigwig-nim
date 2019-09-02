import strutils
import tables
import strformat
import ./lib
import hts/files
import argparse

type region = tuple[chrom: string, start: int, stop: int]

proc isdigit2(s:string): bool {.inline.} =
  for c in s:
    if c < '0' or c > '9': return false
  return true

proc looks_like_region_file(f:string): bool =
  if ':' in f and '-' in f: return false
  if not f.fileExists: return false
  var fh:HTSFile
  if not open(fh, f):
    stderr.write_line &"[bigwig] tried '{f}' as a region file but couldn't open. Trying as an actual region"
    return false
  defer:
    fh.close()
  for l in fh.lines:
    if l[0] == '#' or l.strip().len == 0: continue
    var toks = l.strip().split("\t")
    if toks.len >= 3 and toks[1].isdigit2 and toks[2].isdigit2: return true
    stderr.write_line &"[bigwig] tried '{f}' as a region file but it did not have proper format. Trying as an actual region"
    return false

proc parse_colon_region(reg: string): region {.inline.} =
  let chrom_rest = reg.split(':', maxsplit=1)
  if chrom_rest.len == 1:
    return (chrom_rest[0], 0, -1)
  doAssert chrom_rest.len == 2, ("[bigwig] invalid region:" & reg)
  var ss = chrom_rest[1].split('-')
  result.chrom = chrom_rest[0]
  result.start = max(0, parseInt(ss[0]) - 1)
  result.stop = parseInt(ss[1])
  if result.stop < result.start:
    quit ("[bigwig] ERROR. invalid region:" & reg)

proc parse_one_region(reg:string): region {.inline.} =
  if reg == "": return ("", 0, -1)
  let chrom_rest = reg.rsplit('\t', maxsplit=4)
  if chrom_rest.len == 1:
    return parse_colon_region(reg)
  result.chrom = chrom_rest[0]
  result.start = max(0, parseInt(chrom_rest[1]))
  result.stop = parseInt(chrom_rest[2])
  if result.stop < result.start:
    quit ("[bigwig] ERROR. invalid region:" & reg)

iterator parse_region(reg_or_bed:string): region {.inline.} =
  if reg_or_bed.looks_like_region_file:
    for l in reg_or_bed.hts_lines:
      yield parse_one_region(l.strip(leading=false, chars={'\n', '\r'}))
  else:
    yield parse_one_region(reg_or_bed)


proc from_fai(path: string): BigWigHeader =
  ## create a bigwig header from an fai (fasta index) or a genome file
  for l in path.lines:
    let vals = l.strip().split('\t')
    result.add((name: vals[0], length: parseInt(vals[1]), tid: result.len.uint32))

proc ffloat(f:float, precision:int=5): string {.inline.} =
  result = format_float(f, ffDecimal, precision=precision)
  result = result.strip(leading=false, chars={'0'})
  if result[result.high] == '.': result.setLen(result.high)

proc write_region_from(ofh:File, bw:var BigWig, reg:region) =
  for iv in bw.intervals(reg.chrom, reg.start, reg.stop):
    var v = ffloat(iv.value, precision=5)
    ofh.write_line(&"{reg.chrom}\t{iv.start}\t{iv.stop}\t{v}")

type chunk = seq[tuple[start: int, stop:int, value:float32]]

iterator chunks(bw: var BigWig, reg: region, n:int): chunk =
  var cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)
  for iv in bw.intervals(reg.chrom, reg.start, reg.stop):
    cache.add(iv)
    if cache.len == n:
      yield cache
      cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)

  if cache.len != 0:
    yield cache

proc make_interval(toks: seq[string], col: int): tuple[start: int, stop: int, value: float32] =
  return (parseInt(toks[1]), parseInt(toks[2]), parseFloat(toks[col]).float32)

iterator chunks(bed_path: string, chrom: var string, n:int, value_column: int= 4): chunk =
  let col = value_column - 1

  var cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)
  for l in bed_path.hts_lines:
    let toks = l.strip.split('\t')
    if toks[0] != chrom and cache.len > 0:
      yield cache
      cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)

    chrom = toks[0]
    var iv = make_interval(toks, col)
    # split on large chunks of 0 bases.
    if iv.value == 0 and iv.stop - iv.start > 100 and (cache.len == 0 or iv.stop - iv.start != cache[cache.high].stop - cache[cache.high].start):
      if cache.len > 0:
        yield cache
        cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)
      yield @[iv]
      continue

    cache.add(iv)
    if cache.len == n:
      yield cache
      cache = newSeqOfCap[tuple[start: int, stop:int, value:float32]](n)

  if cache.len != 0:
    yield cache

proc looks_like_single_base(chunk: chunk): bool =
  var n = chunk.len.float32
  if chunk.len < 2: return false
  var nsmall = 0
  var nskip = 0
  var last_stop = chunk[0].start
  var total_bases = 0
  for c in chunk:
    nsmall += int(c.stop - c.start < 8)
    if last_stop > c.start: return false
    nskip += c.start - last_stop
    last_stop = c.stop
    total_bases += c.stop - c.start

  #echo "nskip:", nskip , "p:", nsmall.float32 / n
  return nsmall.float32 / n > 0.80 and nskip == 0

proc looks_like_fixed_span(chunk: chunk): bool =
  if chunk.len < 2: return false
  var sp = chunk[0].stop - chunk[0].start
  result = true
  for i, c in chunk:
    if likely(i < chunk.high) and c.stop - c.start != sp: return false

proc write_fixed_span(ofh: var BigWig, chunk:chunk, chrom: string, span:int) =
  var values = newSeqOfCap[float32](chunk.len)

  # check for end of chrom
  let end_of_chrom = chunk.len > 1 and chunk[chunk.high].stop - chunk[chunk.high].start < span

  for c in chunk:
    for s in countup(c.start, c.stop - 1, span):
      values.add(c.value)

  let eoc_val = values[values.high]
  let eoc_start = chunk[chunk.high].start
  let eoc_stop = chunk[chunk.high].stop

  if end_of_chrom:
    values.setLen(values.high)
  ofh.add(chrom, chunk[0].start.uint32, values, span=span.uint32, step=span.uint32)
  if end_of_chrom:
    # intervals: seq[tuple[start:T, stop: T, value: float32]]
    ofh.add(chrom, @[(start: eoc_start.uint32, stop: eoc_stop.uint32, value: eoc_val)])

proc write_single_base(ofh: var BigWig, chunk:chunk, chrom: string) =
  ofh.write_fixed_span(chunk, chrom, 1)

proc write_region_from(ofh:var BigWig, bw:var BigWig, reg:region, chunksize:int) =
  ## read from bw and write to ofh. try to do this efficiently
  ## read a chunk of a particular size and guess what the best bigwig
  ## representation might be
  for chunk in bw.chunks(reg, chunksize):
    if chunk.looks_like_single_base:
      ofh.write_single_base(chunk, reg.chrom)
    elif chunk.looks_like_fixed_span:
      ofh.write_fixed_span(chunk, reg.chrom, chunk[0].stop - chunk[0].start)
    else:
      ofh.add(reg.chrom, chunk)

proc write_from(ofh:var BigWig, bed_path: string, value_column: int, chunksize:int) =
  ## read from bw and write to ofh. try to do this efficiently
  ## read a chunk of a particular size and guess what the best bigwig
  ## representation might be
  var chrom: string
  for chunk in bed_path.chunks(chrom, n=chunksize, value_column=value_column):
    #echo chunk[0..<min(chunk.len, 4)]
    #echo ""
    if chunk.looks_like_single_base:
      ofh.write_single_base(chunk, chrom)
    elif chunk.looks_like_fixed_span:
      ofh.write_fixed_span(chunk, chrom, chunk[0].stop - chunk[0].start)
    else:
      ofh.add(chrom, chunk)

proc isBig(path: string): bool =
  return bwIsBigWig(path, nil) == 1 or bbIsBigBed(path, nil) == 1

proc stats_main*() =
  var p = newParser("bigwig stats"):
    option("-s", "--stat", choices= @["mean", "coverage", "min", "max", "sum", "header"], default="mean", help="statistic to output. 'header' will show the lengths, mean and coverage for each chromosome in the bigwig.")
    option("--bins", default="1", help="integer number of bins")
    arg("input", nargs=1)
    arg("region", nargs=1, help="BED file or regions or chromosome, or chrom:start-stop region to extract stats")

  var args = commandLineParams()
  if len(args) > 0 and args[0] == "stats":
    args = args[1..args.high]
  if len(args) == 0: args = @["--help"]

  try:
    discard p.parse(args)
  except UsageError:
    echo p.help
    echo "error:", getCurrentExceptionMsg()
    echo "specify a dummy region, even for --stat header"
    quit 1
  let opts = p.parse(args)
  if opts.help:
    quit 0

  var bw: BigWig
  if not bw.open(opts.input):
    quit "[bigwig] unable to open input file"
  defer:
    bw.close

  if opts.stat == "header":
    echo "#chrom\tlength\tmean_depth\tcoverage"
    for h in bw.header:
      var m = bw.stats(h.name, 0, h.length, stat=Stat.mean, nBins=1)
      var c = bw.stats(h.name, 0, h.length, stat=Stat.coverage, nBins=1)
      echo &"{h.name}\t{h.length}\t{ffloat(m[0])}\t{ffloat(c[0], 8)}"
    quit 0

  if opts.region == "":
    echo p.help
    quit "error: region is required"

  var L = {"mean": Stat.mean, "coverage": Stat.coverage, "min": Stat.min, "max": Stat.max, "sum": Stat.sum}.toTable
  var stat = L[opts.stat]
  var bins = parseInt(opts.bins)

  try:
    for region in opts.region.parse_region:
      var st = bw.stats(region.chrom, region.start, region.stop, stat=stat, nBins=bins)
      for v in st:
        echo &"{region.chrom}\t{region.start}\t{region.stop}\t{ffloat(v)}"
  except:
    echo "error:", getCurrentExceptionMsg()
    quit 1

proc view_main*() =

  var p = newParser("bigwig view"):
    option("-r", "--region", help="optional chromosome, or chrom:start-stop region to view")
    option("-c", "--chrom-sizes", help="file indicating chromosome sizes (can be .fai), only used for converting BED->BigWig")
    option("-i", "--value-column", help="column-number (1-based) of the value to encode in to BigWig, only used for encoding BED->BigWig", default="4")
    option("-O", "--output-fmt", choices= @["bed", "bigwig"], default="bed", help="output format")
    option("-o", "--output-file", default="/dev/stdout", help="output bed or bigwig file")
    arg("input", nargs=1)

  var args = commandLineParams()
  if len(args) > 0 and args[0] == "view":
    args = args[1..args.high]
  if len(args) == 0: args = @["--help"]

  let opts = p.parse(args)
  if opts.help:
    quit 0
  if opts.input == "":
    # TODO: check for stdin (can't get libbigwig to open stdin)
    echo p.help
    echo "[bigwig] input file is required"
    quit 2

  let chunksize = 4096
  if opts.input.isBig:

    var bw:BigWig
    if not bw.open(opts.input):
      quit "[bigwig] couldn't open file:" & opts.input

    if opts.output_fmt == "bed":
      #####################
      ### BigWig To BED ###
      #####################
      var ofh: File
      if not ofh.open(opts.output_file, fmWrite):
        quit "[bigwig] couldn't open output file:" & opts.output_file

      if opts.region == "":
        for chrom in bw.header:
          var reg: region = (chrom.name, 0, chrom.length)
          ofh.write_region_from(bw, reg)

      else:
        for region in opts.region.parse_region:
          ofh.write_region_from(bw, region)

      ofh.close

    elif opts.output_fmt == "bigwig":
      ########################
      ### BigWig To BigWig ###
      ########################
      var ofh: BigWig
      if  not ofh.open(opts.output_file, fmWrite):
        quit "[bigwig] couldn't open output bigwig file:" & opts.output_file
      ofh.setHeader(bw.header)
      ofh.writeHeader

      if opts.region == "":
        for chrom in bw.header:
          var reg: region = (chrom.name, 0, chrom.length)
          ofh.write_region_from(bw, reg, chunksize)

      else:
        for region in opts.region.parse_region:
          ofh.write_region_from(bw, region, chunksize)

      ofh.close
    bw.close
  else:
    if opts.chrom_sizes == "":
      quit "[bigwig] --chrom-sizes is required when input is not bigwig."
    if opts.region != "":
      quit "[bigwig] --region is not supported for BED input"
    var h = opts.chrom_sizes.from_fai
    var ofh: BigWig
    if  not ofh.open(opts.output_file, fmWrite):
      quit "[bigwig] couldn't open output bigwig file:" & opts.output_file
    ofh.setHeader(h)
    ofh.writeHeader
    ofh.write_from(opts.input, parseInt(opts.value_column), chunksize)
    ofh.close
