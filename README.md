# bigwig for nim

[![badge](https://img.shields.io/badge/docs-latest-blue.svg)](https://brentp.github.io/libbigwig-nim/bigwig.html)

## Reading

```Nim
var bw: BigWig
bw.open(path, fmRead)

# avoid allocating when possible
var values: seq[float32]
bw.values(values, "chr1", 0, 2222)

for iv in bw.intervals("chr2", 999, 88888): # iterator.
    ## tuple[start: int, stop: int, value: float32]

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
wtr.add("chr1", 100, @[(start: 22, stop: 33, value: 0.01'f32), (start: 44, stop: 55, value: 155'f32)])

# or with, for example a span of 10 bases:
wtr.add("chr1", 10, @[(start: 20, value: 0.01'f32), (start: 30, value: 155'f32)])

# or an array of values with a given span and step:
var values = @[0.1'f32, 0.2, 0.3, 0.4]
wtr.add("chr1", 100, values, span=100, step=200) # 100-200 is 0.1, 300-400 is 0.2 ...
wtr.close()

```
