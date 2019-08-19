#!/bin/bash
set -euo pipefail
nimble install -y
nimble test
