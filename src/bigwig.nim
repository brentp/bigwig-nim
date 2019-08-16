import ./bigwigpkg/lib
import ./bigwigpkg/cli
import tables
import strformat
import os
export lib

proc main*() =
  type pair = object
    f: proc()
    description: string

  var dispatcher = {
    "view": pair(f:view_main, description:"view and convert bigwig"),
    "stats": pair(f:stats_main, description:"extract stats (coverage a gnotate zip file for a given VCF"),
    }.toOrderedTable

  var args = commandLineParams()

  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""

  dispatcher[args[0]].f()

when isMainModule:
  main()

