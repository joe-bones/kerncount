"""Null simulation behind the table in docs/statistics.md (not collected by pytest).

Usage: python tests/sim_null.py <repo root> <data sets per row> [unit:classes,...] [scenario filter] [n list] [out.csv]
e.g.   python tests/sim_null.py . 400 well:auto,field:two "misspecified,uneven" 3,4,6
"""
import sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, sys.argv[1] + "/src")
from kerncount import models, metrics, settings as st_mod

TR = ["Control", "LyoVec", "CL075", "miR-1", "miR-2", "miR-3", "miR-4", "miR-5"]
CTRL = ("Control", "LyoVec")
COUNT = metrics.COUNT


def settings(unit="well", classes="auto"):
    a = dict(reference="LyoVec", normalise_to=["Control", "LyoVec"], treatment_order=TR, vehicle_check="",
             model_unit=unit, classes=classes)
    return st_mod._merge(st_mod.DEFAULTS, {"analysis": a})


def sim(rng, n_n, var_nt, res_sd=None, var_n=0.05, var_well=0.03, var_pos=0.1, base=np.log(70), drop_lyo=False):
    """Poisson-lognormal field counts. var_nt: {treatment: n x treatment variance}; res_sd: {treatment: field
    log-SD} (heterogeneous field dispersion; default 0.45)."""
    rows = []
    for i in range(n_n):
        n = f"n{i + 1}"
        bn = rng.normal(0, np.sqrt(var_n))
        pos = rng.normal(0, np.sqrt(var_pos), 6)
        for t in TR:
            if drop_lyo and i == 1 and t == "LyoVec":
                continue
            bnt = rng.normal(0, np.sqrt(var_nt.get(t, 0))) if var_nt.get(t, 0) > 0 else 0.0
            sd = (res_sd or {}).get(t, 0.45)
            for w in (1, 2):
                bw = rng.normal(0, np.sqrt(var_well))
                for k in range(6):
                    # field log-normal with mean-preserving correction so arithmetic means are equal under the null
                    mu = base + bn + bnt + bw + pos[k] + rng.normal(0, sd) - sd ** 2 / 2
                    rows.append(dict(n=n, treatment=t, well=w, position=k + 1, image=(w - 1) * 6 + k + 1,
                                     timepoint="", value=float(rng.poisson(np.exp(mu)))))
    return pd.DataFrame(rows)


MIR = [t for t in TR if t.startswith("miR")]
SCEN = {
    "treatments 0.1, controls 0": lambda: {t: 0.1 for t in TR if t not in CTRL},
    "no n x treatment variance": lambda: {},
    "treatments 0.3, controls 0": lambda: {t: 0.3 for t in TR if t not in CTRL},
    "misspecified: miR-1 0.3, other miRs 0.02": lambda: {**{t: 0.02 for t in MIR}, "miR-1": 0.3, "CL075": 0.02},
    "controls unequal: LyoVec 0.1, Control 0": lambda: {"LyoVec": 0.1, **{t: 0.1 for t in TR if t not in CTRL}},
    "uneven fields: miR field SD 1.0 vs 0.4": lambda: {t: 0.1 for t in TR if t not in CTRL},
    "vehicle varies: LyoVec 0.3, Control 0, treatments 0": lambda: {"LyoVec": 0.3},
}


def one(args):
    scen, n_n, seed, unit, classes = args
    rng = np.random.default_rng(seed)
    res_sd = {t: (1.0 if t in MIR else 0.4) for t in TR} if scen.startswith("uneven") else None
    d = sim(rng, n_n, SCEN[scen](), res_sd, drop_lyo=(n_n == 4))
    r = models.analyse(d, COUNT, settings(unit, classes), positions_ok=True, leave_one_out=False, density=False)
    if not r.ok:
        return None
    c = r.contrasts[r.contrasts["family"] == "vs LyoVec"]
    m1 = c[c["group"] == "miR-1"]
    return (c["p"].to_numpy(float), float((c["p_adj"] < 0.05).any()), float(m1["p"].iloc[0]) if len(m1) else np.nan)


if __name__ == "__main__":
    root, reps = sys.argv[1], int(sys.argv[2])
    configs = [("well", "auto")]
    if len(sys.argv) > 3:
        configs = [tuple(x.split(":")) for x in sys.argv[3].split(",")]
    out = []
    with Pool(8) as pool:
        for unit, classes in configs:
            cls = {} if classes == "two" else classes
            sfilt = sys.argv[4].split(",") if len(sys.argv) > 4 else None
            nfilt = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else (3, 4, 6)
            for scen in SCEN:
                if sfilt and not any(f in scen for f in sfilt):
                    continue
                for n_n in nfilt:
                    t0 = time.time()
                    jobs = [(scen, n_n, 1000 * n_n + s, unit, cls) for s in range(reps)]
                    res = [r for r in pool.map(one, jobs) if r is not None]
                    ps = np.concatenate([r[0] for r in res])
                    row = dict(unit=unit, classes=classes, scenario=scen, n=n_n, datasets=len(res),
                               per_comparison=float(np.mean(ps < 0.05)), family_wise=float(np.mean([r[1] for r in res])),
                               miR1=float(np.nanmean([r[2] < 0.05 for r in res])), secs=round(time.time() - t0, 1))
                    print(row, flush=True)
                    out.append(row)
    pd.DataFrame(out).to_csv(sys.argv[6] if len(sys.argv) > 6 else "sim_null_results.csv", index=False)
