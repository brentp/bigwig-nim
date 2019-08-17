import ospaths
template thisModuleFile: string = instantiationInfo(fullPaths = true).filename

when fileExists(thisModuleFile.parentDir / "src/bigwig.nim"):
  # In the git repository the Nimble sources are in a ``src`` directory.
  import src/bigwigpkg/version as _
else:
  # When the package is installed, the ``src`` directory disappears.
  import bigwigpkg/version as _

# Package

version       = bigwigVersion
author        = "Brent Pedersen"
description   = "ergonomic wrapper for libbigwig"
license       = "MIT"


# Dependencies

requires "nimbigwig", "argparse", "hts >= 0.2.20"
srcDir = "src"
installExt = @["nim"]

bin = @["bigwig"]

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r --threads:on tests/all"

task docs, "Builds documentation":
  mkDir("docs"/"bigwig")
  for file in @["src/bigwig.nim", "src/bigwigpkg/lib.nim"]:
    var f = file.changefileext("html").split("/")
    var fn = f[f.high]
    exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ fn.changefileext("html") & " " & file

