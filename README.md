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

bw.close
```
