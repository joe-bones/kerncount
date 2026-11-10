"""Review images: original | original + semi-transparent segmentation, with a header banner."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from skimage.segmentation import find_boundaries

COLOURS = {  # RGB 0-1
    "nucleus": (0.10, 0.85, 0.35),
    "merged": (1.00, 0.55, 0.00),
    "debris": (0.85, 0.20, 0.85),
    "out_of_focus": (0.30, 0.55, 1.00),
    "edge": (1.00, 0.90, 0.10),
}
LEGEND = [("counted", "nucleus"), ("possible merge (counted 1)", "merged"), ("debris", "debris"),
          ("out of focus", "out_of_focus"), ("edge (not counted)", "edge")]


def _font(size: int):
    for name in ("DejaVuSans.ttf", "Arial.ttf", "arial.ttf", "Helvetica.ttc"):
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def _stretch(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, [0.5, 99.8])
    return np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1)


def render(img: np.ndarray, labels: np.ndarray, objects: pd.DataFrame, title: str, subtitle: str,
           alpha: float = 0.35) -> Image.Image:
    g = _stretch(img)
    base = np.stack([g] * 3, -1)
    over = base.copy()
    if len(objects):
        lut = np.zeros((int(labels.max()) + 1, 3))
        cat_of = np.full(int(labels.max()) + 1, "", dtype=object)
        for lab, cat in zip(objects["label"].astype(int), objects["category"]):
            lut[lab] = COLOURS.get(cat, (1, 1, 1))
            cat_of[lab] = cat
        m = labels > 0
        over[m] = (1 - alpha) * over[m] + alpha * lut[labels[m]]
        bnd = find_boundaries(labels, mode="inner")
        over[bnd] = lut[labels[bnd]]
    w = img.shape[1]
    pair = np.concatenate([base, np.ones((img.shape[0], 8, 3)), over], 1)
    pim = Image.fromarray((pair * 255).astype(np.uint8))

    fs = max(14, w // 55)
    head = fs * 3 + 10
    canvas = Image.new("RGB", (pim.width, pim.height + head), (20, 20, 20))
    canvas.paste(pim, (0, head))
    d = ImageDraw.Draw(canvas)
    d.text((10, 5), title, fill=(255, 255, 255), font=_font(fs))
    d.text((10, 8 + fs), subtitle, fill=(200, 200, 200), font=_font(int(fs * 0.85)))
    x = w + 18
    small = _font(int(fs * 0.75))
    for text, cat in LEGEND:
        col = tuple(int(255 * v) for v in COLOURS[cat])
        d.rectangle([x, 8 + fs * 2, x + fs // 2, 8 + fs * 2 + fs // 2], fill=col)
        d.text((x + fs // 2 + 4, 4 + fs * 2), text, fill=(220, 220, 220), font=small)
        x += int(d.textlength(text, font=small)) + fs + 10
    return canvas


QC_COLOURS = {"high": (40, 170, 70), "medium": (230, 160, 0), "low": (220, 40, 40)}


def contact_sheet(items: list[tuple[Path, str, str]], out: Path, title: str = "", cols: int = 3,
                  thumb_w: int = 900, quality: int = 85) -> None:
    """items: (review image, label, qc level). Big labels and a QC-coloured border per tile."""
    if not items:
        return
    border, lab_h = 8, 46
    tiles = []
    for path, label, level in items:
        im = Image.open(path)
        im = im.resize((thumb_w, int(im.height * thumb_w / im.width)))
        col = QC_COLOURS.get(level, (120, 120, 120))
        tile = Image.new("RGB", (thumb_w + 2 * border, im.height + lab_h + 2 * border), col)
        tile.paste(im, (border, border + lab_h))
        d = ImageDraw.Draw(tile)
        d.rectangle([border, border, thumb_w + border, border + lab_h], fill=(15, 15, 15))
        d.text((border + 10, border + 6), label, fill=col, font=_font(30))
        tiles.append(tile)
    tw, th = tiles[0].width, max(t.height for t in tiles)
    head = 60 if title else 0
    rows = (len(tiles) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * tw, rows * th + head), (0, 0, 0))
    if title:
        ImageDraw.Draw(sheet).text((12, 10), title, fill=(255, 255, 255), font=_font(36))
    for i, t in enumerate(tiles):
        sheet.paste(t, ((i % cols) * tw, head + (i // cols) * th))
    sheet.save(out, quality=quality)
