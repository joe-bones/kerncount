#!/bin/sh
# Double-click to open kerncount in your browser. Keep this window open while you work.
cd "$(dirname "$0")"
exec .venv/bin/kerncount gui "$@"
