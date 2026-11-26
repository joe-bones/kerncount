"""Synthetic mCherry tumour z-stacks with known analytic volumes (tests only).

A tumour is a solid ellipsoid (the main mass) plus optional small spherical satellites, rendered on an anisotropic
voxel grid (sub-sampled for partial-volume accuracy), blurred with an anisotropic Gaussian point-spread function,
put on a background with noise, scaled by a brightness factor and clipped to 8 bit (so bright tumours saturate,
as in the real stacks). Written as OME-TIFF with physical voxel sizes via tifffile.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi


@dataclass
class Satellite:
    centre_um: tuple[float, float, float]  # (z, y, x)
    radius_um: float

    @property
    def volume_um3(self) -> float:
        return 4 / 3 * math.pi * self.radius_um ** 3


@dataclass
class SynthTumour:
    semi_axes_um: tuple[float, float, float] = (40.0, 120.0, 90.0)  # (z, y, x) of the main-mass ellipsoid
    centre_um: tuple[float, float, float] | None = None  # default: middle of the stack
    voxel_um: tuple[float, float, float] = (2.36, 0.6239, 0.6239)  # (z, y, x)
    shape: tuple[int, int, int] = (48, 512, 512)  # layers, rows, columns
    satellites: list[Satellite] = field(default_factory=list)
    plateau: float = 150.0  # tumour signal above background (grey levels) before brightness scaling
    background: float = 6.0
    noise_sd: float = 4.0
    psf_sigma_um: tuple[float, float, float] = (3.0, 1.2, 1.2)
    brightness: float = 1.0  # uniform scaling of the whole recorded signal (detector gain / expression)
    hollow_frac: float = 0.0  # dark (necrotic) core: inner ellipsoid of this relative size without signal
    seed: int = 0

    @property
    def volume_um3(self) -> float:
        a, b, c = self.semi_axes_um
        return 4 / 3 * math.pi * a * b * c

    def _centre(self) -> np.ndarray:
        if self.centre_um is not None:
            return np.asarray(self.centre_um, float)
        return np.asarray(self.shape, float) * np.asarray(self.voxel_um) / 2

    def occupancy(self, sub: int = 3) -> np.ndarray:
        """Fraction of each voxel inside tumour tissue (main mass + satellites), float32 in [0, 1]."""
        vz, vy, vx = self.voxel_um
        nz, ny, nx = self.shape
        occ = np.zeros(self.shape, np.float32)
        offs = (np.arange(sub) + 0.5) / sub
        zc = (np.arange(nz)[:, None] + offs[None, :]).ravel() * vz
        yc = (np.arange(ny)[:, None] + offs[None, :]).ravel() * vy
        xc = (np.arange(nx)[:, None] + offs[None, :]).ravel() * vx
        cz, cy, cx = self._centre()
        a, b, c = self.semi_axes_um
        ry = ((yc - cy) / b) ** 2
        rx = ((xc - cx) / c) ** 2
        rxy = (ry[:, None] + rx[None, :]).astype(np.float32)
        objects = [((cz, cy, cx), (a, b, c))] + [(s.centre_um, (s.radius_um,) * 3) for s in self.satellites]
        for k, z in enumerate(zc):
            iz = k // sub
            sl = np.zeros((ny * sub, nx * sub), bool)
            for (oz, oy, ox), (ra, rb, rc) in objects:
                dz = ((z - oz) / ra) ** 2
                if dz > 1:
                    continue
                if (oz, oy, ox) == (cz, cy, cx) and (ra, rb, rc) == (a, b, c):
                    sl |= rxy <= 1 - dz
                else:
                    ys = ((yc - oy) / rb) ** 2
                    xs = ((xc - ox) / rc) ** 2
                    sl |= (ys[:, None] + xs[None, :]) <= 1 - dz
            occ[iz] += sl.reshape(ny, sub, nx, sub).mean(axis=(1, 3)) / sub
        return occ

    def stack(self) -> np.ndarray:
        """The recorded 8-bit z-stack (ZYX uint8)."""
        rng = np.random.default_rng(self.seed)
        occ = self.occupancy()
        if self.hollow_frac > 0:
            from dataclasses import replace

            core = replace(self, semi_axes_um=tuple(a * self.hollow_frac for a in self.semi_axes_um),
                           centre_um=tuple(self._centre()), satellites=[], hollow_frac=0.0)
            occ -= core.occupancy()
        img = occ * np.float32(self.plateau)
        sigma = [s / v for s, v in zip(self.psf_sigma_um, self.voxel_um)]
        img = ndi.gaussian_filter(img, sigma)
        img += np.float32(self.background)
        img += rng.normal(0, self.noise_sd, img.shape).astype(np.float32)
        img *= np.float32(self.brightness)
        return np.clip(np.rint(img), 0, 255).astype(np.uint8)


def write_ome(path: Path, stack: np.ndarray, voxel_um: tuple[float, float, float] | None) -> Path:
    """Write a ZYX uint8 stack as OME-TIFF; voxel_um=None writes no physical sizes (uncalibrated)."""
    import tifffile

    md: dict = {"axes": "ZYX"}
    if voxel_um is not None:
        vz, vy, vx = voxel_um
        md.update(PhysicalSizeZ=vz, PhysicalSizeY=vy, PhysicalSizeX=vx,
                  PhysicalSizeZUnit="µm", PhysicalSizeYUnit="µm", PhysicalSizeXUnit="µm")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, stack, ome=True, metadata=md, photometric="minisblack")
    return path
