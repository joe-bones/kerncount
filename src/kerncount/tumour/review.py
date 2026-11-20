"""Blinded review sheets (ticket 06, ADR 0002): one image per tumour + overview sheets, labelled by a random code only.

Blinding: the reviewer gets the `review/` folder, which holds nothing that names a group —
- every tumour gets a random code (`review_key.csv` maps code -> tumour ID; it lives in the project folder, NOT in
  `review/`, and codes stay fixed across runs so corrections keep pointing at the same tumour);
- file names, image titles and `review/corrections.csv` carry only the code (no factor level, session folder, date,
  section label or file name); tumours are listed / tiled in code order, which is random, not grouped;
- each image is contrast-stretched on its own (percentiles), so brightness differences between tumour lines are not
  shown as such; JPEGs are written without metadata.

Per tumour (`review/<code>.jpg`, downsampled to `[tumour] review_px`): left the maximum projection over all layers
with every object's footprint outline, right the mid-section centre layer with the outlines in that layer. Main mass
(object 1) orange, satellites cyan; the largest objects are numbered — the numbers are what a correction re-picks.
"""
from __future__ import annotations

import math
import secrets
from pathlib import Path

import numpy as np
import pandas as pd

from .config import tumour_settings

REVIEW = "review"
KEY_FILE = "review_key.csv"
CORRECTIONS = "corrections.csv"
ALPHABET = "BCDFGHJKLMNPQRSTVWXZ23456789"  # no vowels: codes never spell a word or a level name
CODE_LEN = 5
MAIN_COLOUR, SAT_COLOUR = "#ff9f1c", "#2ec4f0"
MAX_NUMBERED = 40
PER_SHEET = 12


def review_dir(project_dir: Path) -> Path:
    return Path(project_dir) / REVIEW


def read_key(project_dir: Path) -> pd.DataFrame:
    p = Path(project_dir) / KEY_FILE
    if not p.exists():
        return pd.DataFrame(columns=["code", "tumour_id"])
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def codes(project_dir: Path, tumour_ids: list[str]) -> dict[str, str]:
    """{tumour ID: code}; new tumours get a fresh random code, existing codes never change. Writes the key."""
    from ..fileio import atomic_write

    key = read_key(project_dir)
    have = dict(zip(key["tumour_id"], key["code"]))
    used = set(key["code"])
    new = []
    for tid in tumour_ids:
        if tid in have:
            continue
        while True:
            c = "".join(secrets.choice(ALPHABET) for _ in range(CODE_LEN))
            if c not in used:
                break
        used.add(c)
        have[tid] = c
        new.append({"code": c, "tumour_id": tid})
    if new or not (Path(project_dir) / KEY_FILE).exists():
        out = pd.concat([key, pd.DataFrame(new, columns=["code", "tumour_id"])], ignore_index=True)
        atomic_write(Path(project_dir) / KEY_FILE, out.sort_values("code").to_csv(index=False))
    return {t: have[t] for t in tumour_ids}


def _block(img: np.ndarray, f: int, how: str) -> np.ndarray:
    h, w = (img.shape[0] // f) * f, (img.shape[1] // f) * f
    if h == 0 or w == 0:
        return img
    b = img[:h, :w].reshape(h // f, f, w // f, f)
    return b.mean(axis=(1, 3)) if how == "mean" else b.max(axis=(1, 3))


def _label_projection(labels: np.ndarray) -> np.ndarray:
    """2D label image: at each pixel the smallest object number present in any layer (main mass wins)."""
    proj = np.zeros(labels.shape[1:], labels.dtype)
    for lay in labels:
        m = (lay > 0) & ((proj == 0) | (lay < proj))
        proj[m] = lay[m]
    return proj


def _stretch(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, (1.0, 99.7))
    return np.clip((img - lo) / (hi - lo if hi > lo else 1.0), 0, 1)


def render(stack: np.ndarray, labels: np.ndarray, centre: int, code: str, out: Path,
           px_um: tuple[float, float], max_px: int = 640, quality: int = 80) -> Path:
    """Write one blinded review image (JPEG) for a tumour; `labels` already carry any re-picked main mass."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import ndimage as ndi

    centre = min(max(centre, 0), stack.shape[0] - 1)
    f = max(1, math.ceil(max(stack.shape[1:]) / max_px))
    panels = [("maximum projection, all layers", stack.max(axis=0), _label_projection(labels)),
              (f"mid-section centre layer ({centre + 1} of {stack.shape[0]})", stack[centre], labels[centre])]
    n_obj = int(labels.max())
    h, w = _block(panels[0][1], f, "mean").shape
    scale = max(1.0, 700 / (2 * w))  # small stacks: draw at least ~700 px wide so the titles fit
    fig, axes = plt.subplots(1, 2, figsize=((2 * w / 100 + 0.2) * scale, h / 100 * scale + 0.75), dpi=100)
    for ax, (title, img, lab) in zip(axes, panels):
        ax.imshow(_stretch(_block(img.astype(np.float32), f, "mean")), cmap="gray", vmin=0, vmax=1,
                  interpolation="nearest")
        for mask, col, lw in ((lab == 1, MAIN_COLOUR, 1.3), (lab > 1, SAT_COLOUR, 0.9)):
            m = _block(mask, f, "max").astype(np.float32)
            if m.any() and not m.all():
                ax.contour(m, [0.5], colors=col, linewidths=lw)
        ids = [k for k in range(1, min(n_obj, MAX_NUMBERED) + 1)]
        present = np.bincount(lab.ravel(), minlength=n_obj + 1)
        ids = [k for k in ids if present[k]]
        if ids:
            for k, (cy, cx) in zip(ids, ndi.center_of_mass(lab > 0, lab, ids)):
                ax.text(cx / f, cy / f, str(k), color=MAIN_COLOUR if k == 1 else SAT_COLOUR, fontsize=7,
                        ha="center", va="center", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.1", fc="black", ec="none", alpha=0.55))
        ax.set_title(title, fontsize=7, color="#444444")
        ax.set_axis_off()
    bar_px = 100.0 / px_um[1] / f  # 100 µm scale bar
    axes[0].plot([w * 0.04, w * 0.04 + bar_px], [h * 0.95] * 2, color="white", linewidth=2)
    axes[0].text(w * 0.04, h * 0.93, "100 µm", color="white", fontsize=6)
    fig.text(0.01, 1 - 0.08 / fig.get_figheight(), code, fontsize=12, fontweight="bold", va="top")
    fig.text(0.99, 1 - 0.1 / fig.get_figheight(),
             "orange = main mass (object 1), cyan = satellites\nnumber = object to type in corrections.csv",
             fontsize=6.5, va="top", ha="right", color="#444444")
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0.005, top=1 - 0.6 / fig.get_figheight(), wspace=0.02)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, format="jpeg", dpi=100, pil_kwargs={"quality": int(quality), "optimize": True})
    plt.close(fig)
    return out


def overviews(project_dir: Path, codes_in_order: list[str], quality: int = 75, width: int = 600) -> list[Path]:
    """Contact sheets `review/overview_<k>.jpg`: PER_SHEET tumour images each, in code order (3 columns)."""
    from PIL import Image

    rd = review_dir(project_dir)
    imgs = [rd / f"{c}.jpg" for c in sorted(codes_in_order) if (rd / f"{c}.jpg").exists()]
    out = []
    for k in range(0, len(imgs), PER_SHEET):
        tiles = []
        for p in imgs[k:k + PER_SHEET]:
            with Image.open(p) as im:
                rgb = im.convert("RGB")
            rgb.thumbnail((width, width))
            tiles.append(rgb)
        cols = min(3, len(tiles))
        rows = math.ceil(len(tiles) / cols)
        tw, th = max(t.width for t in tiles) + 8, max(t.height for t in tiles) + 8  # white gutter between tiles
        sheet = Image.new("RGB", (cols * tw, rows * th), "white")
        for i, t in enumerate(tiles):
            sheet.paste(t, ((i % cols) * tw, (i // cols) * th))
        path = rd / f"overview_{k // PER_SHEET + 1}.jpg"
        sheet.save(path, quality=quality, optimize=True)
        out.append(path)
    for stale in rd.glob("overview_*.jpg"):
        if stale not in out:
            stale.unlink()
    return out


def write_corrections_template(project_dir: Path, settings: dict, rows: pd.DataFrame, rounds: int = 3) -> Path:
    """`review/corrections.csv`: one row per coded tumour in code order (code, objects, qc_flags, manual_1..,
    note). What the reviewer typed (manual_*, note) is kept; rows of tumours no longer measured stay as they are."""
    from ..fileio import atomic_write

    path = review_dir(project_dir) / CORRECTIONS
    path.parent.mkdir(parents=True, exist_ok=True)
    delim = settings.get("project", {}).get("csv_delimiter", ",")
    old = (pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False) if path.exists()
           else pd.DataFrame(columns=["code"]))
    user = [c for c in old.columns if c.startswith("manual_") or c == "note"]
    manual = sorted({*(f"manual_{i}" for i in range(1, rounds + 1)), *(c for c in user if c != "note")},
                    key=lambda c: int(c.split("_")[1]) if c.split("_")[1].isdigit() else 10 ** 6)
    fresh = rows[["code", "objects", "qc_flags"]].astype(str)
    keep_old = old[~old["code"].isin(fresh["code"])]
    out = fresh.merge(old[["code", *user]], on="code", how="left")
    out = pd.concat([out, keep_old], ignore_index=True)
    for c in [*manual, "note"]:
        if c not in out.columns:
            out[c] = ""
    out = out[["code", "objects", "qc_flags", *manual, "note"]].fillna("").sort_values("code")
    atomic_write(path, out.to_csv(index=False, sep=delim))
    return path


def options(settings: dict) -> tuple[int, int]:
    """(review_px, review_jpeg_quality) from [tumour]."""
    ts = tumour_settings(settings)
    return int(ts.get("review_px", 640)), int(ts.get("review_jpeg_quality", 80))
