#!/bin/sh
set -eu
script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$script_directory/../run-example.sh" "$(basename -- "$script_directory")" "$@"
