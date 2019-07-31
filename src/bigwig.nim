import nimbigwig/bigWig
import os

import ./bigwigpkg/version

type BigWig* = ref object
  bw : ptr bigWigFile_t
  path: string

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

