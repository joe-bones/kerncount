#!/bin/sh
# One-time install (Mac/Linux). Needs Python 3.12.
cd "$(dirname "$0")"
python3.12 -m venv .venv && .venv/bin/python -m pip install --upgrade pip && .venv/bin/python -m pip install -e . \
  && echo "Done. Run:  $(pwd)/.venv/bin/kerncount run <project folder>"
