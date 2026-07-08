#!/usr/bin/env bash
set -euo pipefail

unset LD_LIBRARY_PATH

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec julia --project="${repo_root}/julia" "${repo_root}/scripts/interpolate_field.jl" "$@"
