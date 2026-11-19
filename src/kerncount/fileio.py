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
        os.chmod(tmp, _mode(path))  # mkstemp creates 0600; keep the file's permissions (or the umask default)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _mode(path: Path) -> int:
    try:
        return Path(path).stat().st_mode & 0o777
    except OSError:
        um = os.umask(0)
        os.umask(um)
        return 0o666 & ~um


class FolderInUse(RuntimeError):
    """The old folder could not be replaced (on Windows: a file in it is open, typically in Excel)."""


KEPT = "_kept"  # files found in an output folder that kerncount did not write: copied here, never deleted


def _files(root: Path) -> set[str]:
    return {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()} if root.is_dir() else set()


def swap_in(tmp: Path, final: Path, known: set[str] | None) -> list[str]:
    """Replace the folder `final` by the freshly built `tmp` (same parent) with renames only: a build that fails
    or a file held open elsewhere never leaves a half-written folder behind.

    known: relative paths kerncount itself wrote into the old `final` (from its manifest); None = no record, so
    every file the new build does not rewrite counts as a user file. User files are copied into
    <final>/_kept/<relative path> (earlier _kept files stay), never deleted. Returns the kept paths.
    Raises FolderInUse when the old folder cannot be moved (Windows: a file in it is open, e.g. in Excel); the new
    results then stay in `tmp` and the old folder is untouched."""
    import shutil

    tmp, final = Path(tmp), Path(final)
    kept: list[str] = []
    if not final.exists():
        os.replace(tmp, final)
        return kept
    for rel in sorted(_files(final) - _files(tmp)):
        if rel.startswith(KEPT + "/"):
            dst = tmp / rel
        elif known is not None and rel in known:
            continue  # an output of an earlier build that this build no longer writes
        else:
            dst = tmp / KEPT / rel
            kept.append(rel)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(final / rel, dst)
    old = final.with_name(f".{final.name}.old")
    if old.exists():
        shutil.rmtree(old, ignore_errors=True)
    try:
        os.replace(final, old)
    except PermissionError as e:
        raise FolderInUse(f"cannot replace {final.name}: is a file in {final.name} open in Excel (or another "
                          f"program)? Close it and run again; the new results wait in {tmp.name}.") from e
    try:
        os.replace(tmp, final)
    except OSError:
        os.replace(old, final)  # put the old folder back
        raise
    shutil.rmtree(old, ignore_errors=True)
    return kept
