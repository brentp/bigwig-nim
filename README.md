# bigwig for nim

[![badge](https://img.shields.io/badge/docs-latest-blue.svg)](https://brentp.github.io/libbigwig-nim/bigwig.html)

## Command Line

bigwig-nim includes a command-line tool distributed as a static binary. It supports converting
bed to bigwig and bigwig to bed and extracting stats (mean, coverage, etc) for regions in a bigwig.

There are other tools to do this, including [kentTools](https://hgwdev.gi.ucsc.edu/~kent/src/) which has a more restrictive license and does not supported (b)gzipped input and [bwtools](https://github.com/CRG-Barcelona/bwtool) which seems to provide similar functionality (but I am not able to build it).


### view 

To convert a bed with the value in the 4th column to bigwig, use:

```Shell
bigwig view $bed_in --value-column 4 --chrom-sizes $fai -O bigwig -o $bigwig_out
```
`bigwig` will automatically determine the best data format for each block (fixed span and step or per-base) most of the
CPU time is spent parsing the input bed file.

### stats

To get the mean value for a given region (in this case on chromosome 22)

```Shell
bigwig stats --stat mean $bigwig 22:145000-155000
```

The supported stats are `mean`, `min`, `max`, `coverage`, `sum` with a special-case for the stat of `header` which
shows the chromosomes, lengths and mean coverages for each chromosome in the bigwig file.


## Reading

```Nim
var bw: BigWig
bw.open(path, fmRead)

# avoid allocating when possible
var values: seq[float32]
bw.values(values, "chr1", 0, 2222)

for iv in bw.intervals("chr2", 999, 88888): # iterator.
  # tuple[start: int, stop: int, value: float32]

# for bigbed
for iv in bw.entries("chr2", 999, 88888): # iterator.
  # tuple[start: int, stop: int, value: cstring]
  # value contains "SQL" for bigbed entry.

# single value
var m: seq[float32] = bw.stats("chr2", 999, 9999, stat=Stat.mean)

# multiple bins:
var L: seq[float32] = bw.stats("chr2", 999, 9999, stat=Stat.min, nBins=10)

echo bw.header # @[(name: "1", length: 195471971, tid: 0'u32), (name: "10", length: 130694993, tid: 1'u32)]

bw.close
```

## Writing

```Nim
var wtr:BigWig
doAssert wtr.open("tests/writer.bw", fmWrite)
wtr.setHeader(@[(name:"chr1", length: 2000, tid: 0'u32)])
wtr.writeHeader

# add intervals with tuples
wtr.add("chr1", @[(start: 22, stop: 33, value: 0.01'f32), (start: 44, stop: 55, value: 155'f32)])

# or with, for example a span of 15 bases:
wtr.add("chr1", 15, @[(start: 20, value: 0.01'f32), (start: 30, value: 155'f32)])

# or an array of values with a given span and step:
var values = @[0.1'f32, 0.2, 0.3, 0.4]
wtr.add("chr1", 100, values, span=100, step=200) # 100-200 is 0.1, 300-400 is 0.2 ...
wtr.close()

```
