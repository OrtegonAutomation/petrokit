"""
petrokit.nodal
==============

Análisis Nodal: intersección IPR–VLP.

Compatibilidad (API Fase 1, según README):
- nodal_analysis(p_res, q_max, well_depth, rho, mu, d, npts=50) -> (q_op, pwf_op)
- plot_nodal(...)

Extensión (PR2 Paso 2):
- Permite seleccionar modelos por nombre:
    ipr_model: "vogel" | "fetkovich" | "jones" | "standing"
    vlp_model: "darcy" | "beggs_brill" | "hagedorn_brown"

Notas:
- Unidades esperadas:
    * p_res, pwf: psi
    * q: STB/d
    * well_depth: ft
    * rho: lb/ft^3
    * mu: cP
    * d: in
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


def _norm_model(name: str) -> str:
    return (name or "").strip().lower().replace("-", "_").replace(" ", "_")


def _require_positive(name: str, value: float, allow_zero: bool = False) -> None:
    if allow_zero:
        if not (value >= 0):
            raise ValueError(f"{name} debe ser >= 0. Recibido: {value}")
    else:
        if not (value > 0):
            raise ValueError(f"{name} debe ser > 0. Recibido: {value}")


@dataclass(frozen=True)
class NodalResult:
    q_op: float
    pwf_op: float
    q_grid: np.ndarray
    pwf_ipr_on_q: np.ndarray
    pwf_vlp_on_q: np.ndarray
    ipr_model: str
    vlp_model: str


def _interp_pwf_from_ipr_curve(pwf_ipr: np.ndarray, q_ipr: np.ndarray, q_grid: np.ndarray) -> np.ndarray:
    pwf_ipr = np.asarray(pwf_ipr, dtype=float)
    q_ipr = np.asarray(q_ipr, dtype=float)

    # Para np.interp, x debe ir creciente: q_sorted crece de 0 -> qmax
    q_sorted = q_ipr[::-1]
    pwf_sorted = pwf_ipr[::-1]

    q_min, q_max = float(np.min(q_sorted)), float(np.max(q_sorted))
    q_clip = np.clip(q_grid, q_min, q_max)
    return np.interp(q_clip, q_sorted, pwf_sorted)


def _find_intersection(q_grid: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> Tuple[float, float]:
    """
    Intersección y1(q)=y2(q) sobre q_grid.
    Si no hay intersección:
      - si y1<y2 en todo q -> no-flow: q=0, pwf=y1(q=0)
      - si y1>y2 en todo q -> opera al máximo de la grilla: q=qmax, pwf=y1(qmax)
    """
    q_grid = np.asarray(q_grid, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    y2 = np.asarray(y2, dtype=float)

    diff = y1 - y2
    finite = np.isfinite(diff) & np.isfinite(q_grid) & np.isfinite(y1) & np.isfinite(y2)
    if not np.any(finite):
        raise ValueError("No hay datos finitos para calcular la intersección IPR–VLP.")

    qg = q_grid[finite]
    d = diff[finite]
    y1g = y1[finite]

    zeros = np.where(np.isclose(d, 0.0, atol=1e-12))[0]
    if zeros.size > 0:
        i0 = int(zeros[0])
        return float(qg[i0]), float(y1g[i0])

    sign = np.sign(d)
    changes = np.where(sign[:-1] * sign[1:] < 0)[0]
    if changes.size > 0:
        i = int(changes[0])
        q1, q2 = qg[i], qg[i + 1]
        d1, d2 = d[i], d[i + 1]
        t = d1 / (d1 - d2)

        q_op = q1 + t * (q2 - q1)
        y1_op = y1g[i] + t * (y1g[i + 1] - y1g[i])
        return float(q_op), float(y1_op)

    if np.all(d < 0):
        return float(qg[0]), float(y1g[0])

    if np.all(d > 0):
        return float(qg[-1]), float(y1g[-1])

    j = int(np.argmin(np.abs(d)))
    return float(qg[j]), float(y1g[j])


def nodal_analysis(
    p_res: float,
    q_max: float,
    well_depth: float,
    rho: float,
    mu: float,
    d: float,
    npts: int = 50,
    ipr_model: str = "vogel",
    vlp_model: str = "darcy",
    ipr_kwargs: Optional[Dict[str, Any]] = None,
    vlp_kwargs: Optional[Dict[str, Any]] = None,
) -> Tuple[float, float]:
    res = nodal_analysis_detail(
        p_res=p_res,
        q_max=q_max,
        well_depth=well_depth,
        rho=rho,
        mu=mu,
        d=d,
        npts=npts,
        ipr_model=ipr_model,
        vlp_model=vlp_model,
        ipr_kwargs=ipr_kwargs,
        vlp_kwargs=vlp_kwargs,
    )
    return res.q_op, res.pwf_op


def nodal_analysis_detail(
    p_res: float,
    q_max: float,
    well_depth: float,
    rho: float,
    mu: float,
    d: float,
    npts: int = 50,
    ipr_model: str = "vogel",
    vlp_model: str = "darcy",
    ipr_kwargs: Optional[Dict[str, Any]] = None,
    vlp_kwargs: Optional[Dict[str, Any]] = None,
) -> NodalResult:
    _require_positive("p_res", float(p_res))
    _require_positive("q_max", float(q_max), allow_zero=True)
    _require_positive("well_depth", float(well_depth))
    _require_positive("rho", float(rho))
    _require_positive("mu", float(mu))
    _require_positive("d", float(d))
    if int(npts) < 10:
        raise ValueError("npts debe ser >= 10 para nodal_analysis_detail.")

    ipr_kwargs = {} if ipr_kwargs is None else dict(ipr_kwargs)
    vlp_kwargs = {} if vlp_kwargs is None else dict(vlp_kwargs)

    ipr_m = _norm_model(ipr_model)
    vlp_m = _norm_model(vlp_model)

    # ---- IPR
    from . import ipr as ipr_mod

    if ipr_m == "vogel":
        pwf_ipr, q_ipr = ipr_mod.ipr_curve_vogel(p_res=float(p_res), q_max=float(q_max), npts=int(npts))
    elif ipr_m == "fetkovich":
        if "J" not in ipr_kwargs:
            raise ValueError("ipr_model='fetkovich' requiere ipr_kwargs={'J': ...}.")
        pwf_ipr, q_ipr = ipr_mod.ipr_curve_fetkovich(p_res=float(p_res), J=float(ipr_kwargs["J"]), npts=int(npts))
    elif ipr_m == "jones":
        if "C" not in ipr_kwargs or "D" not in ipr_kwargs:
            raise ValueError("ipr_model='jones' requiere ipr_kwargs={'C': ..., 'D': ...}.")
        pwf_ipr, q_ipr = ipr_mod.ipr_curve_jones(
            p_res=float(p_res), C=float(ipr_kwargs["C"]), D=float(ipr_kwargs["D"]), npts=int(npts)
        )
    elif ipr_m == "standing":
        if "p_b" not in ipr_kwargs or "J" not in ipr_kwargs:
            raise ValueError("ipr_model='standing' requiere ipr_kwargs={'p_b': ..., 'J': ...}.")
        pwf_ipr, q_ipr = ipr_mod.ipr_curve_standing(
            p_res=float(p_res), p_b=float(ipr_kwargs["p_b"]), J=float(ipr_kwargs["J"]), npts=int(npts)
        )
    else:
        raise ValueError("ipr_model inválido. Usa: vogel, fetkovich, jones, standing.")

    pwf_ipr = np.asarray(pwf_ipr, dtype=float)
    q_ipr = np.asarray(q_ipr, dtype=float)

    # ---- q grid
    q_hi = float(q_max) if float(q_max) > 0 else float(np.max(q_ipr))
    if q_hi <= 0:
        q_hi = 1.0
    q_grid = np.linspace(0.0, q_hi, int(npts))

    pwf_ipr_on_q = _interp_pwf_from_ipr_curve(pwf_ipr, q_ipr, q_grid)

    # ---- VLP
    from .vlp import vlp_curve

    try:
        from .vlp import vlp_curve_model  # ya lo tienes en vlp.py
        pwf_vlp_on_q = vlp_curve_model(
            vlp_m,
            q_grid,
            well_depth=float(well_depth),
            rho=float(rho),
            mu=float(mu),
            d=float(d),
            **vlp_kwargs,
        )
    except (ImportError, AttributeError):
        f = float(vlp_kwargs.get("f", 0.02))
        pwf_vlp_on_q = vlp_curve(q_grid, float(well_depth), float(rho), float(mu), float(d), f=f)

    pwf_vlp_on_q = np.asarray(pwf_vlp_on_q, dtype=float)

    q_op, pwf_op = _find_intersection(q_grid, pwf_ipr_on_q, pwf_vlp_on_q)

    return NodalResult(
        q_op=float(q_op),
        pwf_op=float(pwf_op),
        q_grid=q_grid,
        pwf_ipr_on_q=pwf_ipr_on_q,
        pwf_vlp_on_q=pwf_vlp_on_q,
        ipr_model=ipr_m,
        vlp_model=vlp_m,
    )


def plot_nodal(
    p_res: float,
    q_max: float,
    well_depth: float,
    rho: float,
    mu: float,
    d: float,
    npts: int = 50,
    ipr_model: str = "vogel",
    vlp_model: str = "darcy",
    ipr_kwargs: Optional[Dict[str, Any]] = None,
    vlp_kwargs: Optional[Dict[str, Any]] = None,
    show: bool = True,
) -> NodalResult:
    res = nodal_analysis_detail(
        p_res=p_res,
        q_max=q_max,
        well_depth=well_depth,
        rho=rho,
        mu=mu,
        d=d,
        npts=npts,
        ipr_model=ipr_model,
        vlp_model=vlp_model,
        ipr_kwargs=ipr_kwargs,
        vlp_kwargs=vlp_kwargs,
    )

    plt.figure()
    plt.plot(res.q_grid, res.pwf_ipr_on_q, label=f"IPR ({res.ipr_model})")
    plt.plot(res.q_grid, res.pwf_vlp_on_q, label=f"VLP ({res.vlp_model})")
    plt.scatter([res.q_op], [res.pwf_op], zorder=3, label="Punto de operación")
    plt.xlabel("q (STB/d)")
    plt.ylabel("pwf (psi)")
    plt.title("Análisis Nodal (IPR–VLP)")
    plt.grid(True)
    plt.legend()
    if show:
        plt.show()
    return res


__all__ = [
    "NodalResult",
    "nodal_analysis",
    "nodal_analysis_detail",
    "plot_nodal",
]
