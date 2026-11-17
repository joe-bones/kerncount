"""Small, standard-library-only file helpers shared by the browser app and the pipeline: atomic writes, per-file
locks, and building folders next to their final place before swapping them in."""
from __future__ import annotations

import os
import tempfile
import threading
from pathlib import Path

_FILE_LOCKS: dict[str, threading.RLock] = {}
_GUARD = threading.Lock()


def file_lock(path: Path) -> threading.RLock:
    """One re-entrant lock per file, for read-modify-write from several threads of one process (the browser app).
    Readers in other processes are safe because every write is atomic (atomic_write)."""
    key = str(Path(path).absolute())
    with _GUARD:
        return _FILE_LOCKS.setdefault(key, threading.RLock())


def atomic_write(path: Path, text: str, encoding: str = "utf-8", newline: str | None = None):
    """Write via a temporary file in the same folder, then os.replace: a reader never sees a half-written file and
    a failed write leaves the old file untouched."""
    path = Path(path)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline=newline) as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
