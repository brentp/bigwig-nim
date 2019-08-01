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

requires "nimbigwig"
srcDir = "src"
installExt = @["nim"]

bin = @["bigwig"]

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r --threads:on tests/test_read"

task docs, "Builds documentation":
  mkDir("docs"/"bigwig")
  var file = "src/bigwig.nim"
  exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file
  for file in listfiles("src/bigwig"):
    if file.endswith("value.nim"): continue
    if splitfile(file).ext == ".nim":
      exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file

