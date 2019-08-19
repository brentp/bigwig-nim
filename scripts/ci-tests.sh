#!/bin/sh
set -euo pipefail
nimble install -y
nimble test
