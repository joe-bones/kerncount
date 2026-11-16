"""REML for a linear mixed model with independent variance components nested in groups (one group = one n).

    y_g = X_g beta + sum_k M_kg b_kg + e_g,   b_kg ~ N(0, theta_k I),   e_g ~ N(0, theta_e I)
    V_g = sum_k theta_k M_kg M_kg' + theta_e I

statsmodels' MixedLM builds the design (from_formula with vc_formula) and supplies one starting point, but its
optimiser works on square-root parameters and its inner solver breaks down when a variance is exactly zero - a
common, legitimate REML estimate here (very stable controls). This module maximises the same REML likelihood
directly over theta >= 0 (L-BFGS-B with bounds, analytic gradient, several starting points, best likelihood
wins), so boundary estimates are handled cleanly. It also gives the GLS covariance of beta, the expected REML
information of theta, Satterthwaite degrees of freedom for any L beta and the BLUPs of the random effects.
All work is done group by group with the low-rank M matrices, so memory grows with the largest group only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg, optimize


@dataclass
class Group:
    y: np.ndarray
    X: np.ndarray
    M: list[np.ndarray]  # one (N, q_k) matrix per random term (q_k may be 0)


class REML:
    def __init__(self, groups: list[Group], names: list[str]):
        self.groups = groups
        self.names = list(names) + ["residual"]  # theta order: random terms, then residual
        self.K = len(names)
        self.p = groups[0].X.shape[1]
        self.N = sum(len(g.y) for g in groups)
        ally = np.concatenate([g.y for g in groups])
        allX = np.vstack([g.X for g in groups])
        beta, *_ = np.linalg.lstsq(allX, ally, rcond=None)
        r = ally - allX @ beta
        self.var_y = float(max(r.var(), 1e-12))

    # ------------------------------------------------------------ core pieces
    def _V(self, g: Group, theta) -> np.ndarray:
        V = theta[-1] * np.eye(len(g.y))
        for t, M in zip(theta[:-1], g.M):
            if M.size and t > 0:
                V += t * (M @ M.T)
        return V

    def _solve(self, theta):
        """Per group: Vi, logdet; plus C = (X'V^-1X)^-1, beta and pieces for gradients."""
        parts = []
        XtViX = np.zeros((self.p, self.p))
        XtViy = np.zeros(self.p)
        logdet = 0.0
        for g in self.groups:
            cf = linalg.cho_factor(self._V(g, theta), lower=True, check_finite=False)
            Vi = linalg.cho_solve(cf, np.eye(len(g.y)), check_finite=False)
            logdet += 2 * np.log(np.diag(cf[0])).sum()
            W = Vi @ g.X
            XtViX += g.X.T @ W
            XtViy += W.T @ g.y
            parts.append((Vi, W))
        cfx = linalg.cho_factor(XtViX, lower=True, check_finite=False)
        C = linalg.cho_solve(cfx, np.eye(self.p), check_finite=False)
        beta = C @ XtViy
        logdet_x = 2 * np.log(np.diag(cfx[0])).sum()
        return parts, C, beta, logdet, logdet_x

    def loglike(self, theta, grad: bool = False):
        theta = np.asarray(theta, float)
        try:
            parts, C, beta, logdet, logdet_x = self._solve(theta)
        except (np.linalg.LinAlgError, linalg.LinAlgError, ValueError):
            return (-np.inf, np.zeros(len(theta))) if grad else -np.inf
        quad = 0.0
        g_out = np.zeros(self.K + 1)
        T = [np.zeros((self.p, self.p)) for _ in range(self.K + 1)]
        for g, (Vi, W) in zip(self.groups, parts):
            r = g.y - g.X @ beta
            s = Vi @ r
            quad += r @ s
            if grad:
                for k, M in enumerate(g.M):
                    if not M.size:
                        continue
                    ViM = Vi @ M
                    MW = M.T @ W
                    g_out[k] += np.sum(ViM * M) - np.sum((M.T @ s) ** 2)
                    T[k] += MW.T @ MW
                g_out[-1] += np.trace(Vi) - s @ s
                T[-1] += W.T @ W
        ll = -0.5 * (logdet + logdet_x + quad)
        if not grad:
            return ll
        for k in range(self.K + 1):
            g_out[k] -= np.sum(C * T[k])  # tr(C T_k)
        return ll, -0.5 * g_out

    # ------------------------------------------------------------ fitting
    def starts(self, extra: list | None = None) -> list[np.ndarray]:
        v = self.var_y
        K = self.K
        out = [np.r_[np.full(K, 0.1 * v), 0.5 * v], np.r_[np.full(K, v / (K + 1)), v / (K + 1)],
               np.r_[np.full(K, 0.01 * v), 0.9 * v]]
        for e in extra or []:
            e = np.asarray(e, float)
            if e.shape == (K + 1,) and np.all(np.isfinite(e)):
                out.insert(0, np.clip(e, 0, None) + np.r_[np.zeros(K), 1e-9])
        return out

    def fit(self, extra_starts: list | None = None, n_default: int | None = None) -> "REMLFit":
        lo = 1e-8 * self.var_y
        bounds = [(0.0, None)] * self.K + [(lo, None)]
        best = None
        tried = self.starts(extra_starts)
        if n_default is not None:
            tried = tried[: len(extra_starts or []) + n_default]
        for x0 in tried:
            x0 = np.maximum(x0, [b[0] for b in bounds])
            def f(th):
                ll, gr = self.loglike(th, grad=True)
                return (np.inf, np.zeros_like(gr)) if not np.isfinite(ll) else (-ll, -gr)

            try:
                res = optimize.minimize(f, x0, jac=True, method="L-BFGS-B", bounds=bounds,
                                        options=dict(maxiter=500, ftol=1e-12, gtol=1e-7))
            except Exception:  # noqa: BLE001
                continue
            if not np.isfinite(res.fun):
                continue
            if best is None or res.fun < best.fun - 1e-9:
                best = res
        if best is None:
            raise np.linalg.LinAlgError("REML optimisation failed from every starting point")
        theta = np.asarray(best.x, float)
        theta[: self.K][theta[: self.K] < 1e-7 * self.var_y] = 0.0  # numerically zero -> exactly on the boundary
        ll, gr = self.loglike(theta, grad=True)
        # projected gradient: zero on interior parameters, non-positive (pointing outward) at the boundary
        pg = np.where(theta > 0, gr, np.maximum(gr, 0))
        converged = bool(best.success or np.max(np.abs(pg) * np.maximum(theta, self.var_y)) < 1e-3)
        return REMLFit(self, theta, ll, converged, str(best.message))


class REMLFit:
    def __init__(self, model: REML, theta: np.ndarray, llf: float, converged: bool, message: str):
        self.model = model
        self.theta = theta
        self.llf = llf
        self.converged = converged
        self.message = message
        parts, C, beta, _, _ = model._solve(theta)
        self.C, self.beta = C, beta
        self._parts = parts
        self._info = None
        self._T = None

    def _prep(self):
        m, th = self.model, self.theta
        K1 = m.K + 1
        A = np.zeros((K1, K1))
        B = np.zeros((K1, K1))
        T = [np.zeros((m.p, m.p)) for _ in range(K1)]
        for g, (Vi, W) in zip(m.groups, self._parts):
            Ms = list(g.M) + [np.eye(len(g.y))]
            MW = [M.T @ W if M.size else np.zeros((0, m.p)) for M in Ms]
            ViM = [Vi @ M if M.size else np.zeros((len(g.y), 0)) for M in Ms]
            for k in range(K1):
                T[k] += MW[k].T @ MW[k]
            for k in range(K1):
                for l in range(k, K1):
                    if not Ms[k].size or not Ms[l].size:
                        continue
                    Q = Ms[k].T @ ViM[l]  # M_k' Vi M_l
                    A[k, l] += np.sum(Q * Q)
                    b = np.sum((MW[l].T @ Q.T @ MW[k]) * self.C.T)  # tr(C (M_l'W)' (M_l'Vi M_k) (M_k'W))
                    B[k, l] += b
        info = np.zeros((K1, K1))
        for k in range(K1):
            for l in range(k, K1):
                t3 = np.sum((self.C @ T[k]) * (self.C @ T[l]).T)
                info[k, l] = info[l, k] = 0.5 * (A[k, l] - 2 * B[k, l] + t3)
        self._info, self._T = info, T

    @property
    def info(self) -> np.ndarray:
        if self._info is None:
            self._prep()
        return self._info

    def theta_cov(self) -> tuple[np.ndarray, np.ndarray]:
        """(mask of parameters not on the boundary, inverse expected information for those)."""
        keep = self.theta > 0
        return keep, np.linalg.pinv(self.info[np.ix_(keep, keep)])

    def satterthwaite(self, L: np.ndarray) -> float:
        if self._T is None:
            self._prep()
        v = float(L @ self.C @ L)
        LC = L @ self.C
        g = np.array([-(LC @ T @ LC) for T in self._T])
        keep, Ii = self.theta_cov()
        if not keep.any():
            return np.inf
        denom = float(g[keep] @ Ii @ g[keep])
        return 2 * v * v / denom if denom > 0 else np.inf

    def blups(self) -> list[list[np.ndarray]]:
        """Per group, per random term: the predicted random effects theta_k M_k' V^-1 (y - X beta)."""
        out = []
        for g, (Vi, _) in zip(self.model.groups, self._parts):
            s = Vi @ (g.y - g.X @ self.beta)
            out.append([t * (M.T @ s) if M.size else np.zeros(0) for t, M in zip(self.theta[:-1], g.M)])
        return out

    def fitted_resid(self) -> tuple[np.ndarray, np.ndarray]:
        """Conditional fitted values (fixed + random effects) and residuals, in group order."""
        fit, res = [], []
        for g, b in zip(self.model.groups, self.blups()):
            f = g.X @ self.beta + sum((M @ bb for M, bb in zip(g.M, b) if M.size), np.zeros(len(g.y)))
            fit.append(f)
            res.append(g.y - f)
        return np.concatenate(fit), np.concatenate(res)
