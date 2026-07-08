#!/usr/bin/env bash
set -euo pipefail
exec "$(dirname "$0")/run_interpolate_field.sh" "$@"
