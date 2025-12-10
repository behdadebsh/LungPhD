#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimal WSS pipeline (with order-element lists exposed):

What it does
------------
- Input: subject_id (for titles), list of scenario directories, labels
- Output: 6 figures (overlay across scenarios)
    1) Steady WSS (mean±SD) vs Strahler
    2) Steady WSS (mean±SD) vs Horsfield
    3) TAWSS (median±IQR/2) vs Strahler
    4) TAWSS (median±IQR/2) vs Horsfield
    5) OSI (median±IQR/2) vs Strahler
    6) OSI (median±IQR/2) vs Horsfield

- Additionally returns (in the run_analysis(...) return dict):
    results["orders"]["strahler"]  -> {order: [element_ids]}
    results["orders"]["horsfield"] -> {order: [element_ids]}

Assumptions
-----------
Each scenario directory contains:
  * WSS_arterial.txt                 (per-element time series: "<eid> v1 v2 ... vN")
  * lobe_imped.json                  (with "frequency": [f0, f1, ...], uses freqs[1])
  * <subject>_flow_perf.exelem       (steady flow per element)
  * <subject>_radius_perf.exelem     (steady radius per element)

Topology files live in: {subject_path}/Vessel/{subject}_Artery_Full.exelem|.exnode
"""

import os
import re
import math
import json
from collections import defaultdict
from typing import Any, Dict, List, Tuple, Optional, Union

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# ---------------------- CONFIG (edit) -----------------------
# ============================================================

data_path = '/hpc/bsha219/lung/Data'
study = 'CTEPH'
subject = 'Alfred4'
volume = 'Pre'
subject_path = os.path.join(data_path, study, subject, volume)

# Topology & steady scalar file names
element_file_name = subject + '_Artery_Full.exelem'
node_file_name    = subject + '_Artery_Full.exnode'
flow_file_name    = subject + '_flow_perf.exelem'
radius_file_name  = subject + '_radius_perf.exelem'

# Per-scenario filenames (inside each scenario dir)
frequency_json_filename   = 'lobe_imped.json'
wss_filename_in_scenario  = 'WSS_arterial.txt'

# Defaults (can be overridden via run_analysis args)
DEFAULT_NUM_TIMEPOINTS   = 50      # or None to infer from file
DEFAULT_STORE_WSS_SERIES = False   # we don't need to store series for plotting

# Physical constants
MU = 0.004        # Pa·s
MM3_TO_M3 = 1e-9  # mm^3/s → m^3/s
MM_TO_M  = 1e-3   # mm → m

# ============================================================
# ------------------------ Utilities -------------------------
# ============================================================

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def extract_fundamental_frequency(json_path) -> Optional[float]:
    if not os.path.exists(json_path):
        return None
    with open(json_path, "r") as f:
        data = json.load(f)
    freqs = data.get("frequency", [])
    if isinstance(freqs, list) and len(freqs) > 1 and freqs[1]:
        try:
            return float(freqs[1])
        except Exception:
            return None
    return None

def parse_wss_series_file(wss_file_path, num_timepoints: Optional[int] = None) -> Tuple[Dict[int, np.ndarray], int]:
    series: Dict[int, np.ndarray] = {}
    if not os.path.exists(wss_file_path):
        return series, 0
    with open(wss_file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                eid = int(parts[0])
            except ValueError:
                continue
            vals = np.array([safe_float(x) for x in parts[1:]], dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            if num_timepoints is not None and vals.size != num_timepoints:
                continue
            series[eid] = vals
    if num_timepoints is None and series:
        num_timepoints = next(iter(series.values())).size
    return series, (num_timepoints or 0)

def compute_wss_metrics_from_series(w: np.ndarray, dt: Optional[float] = None) -> Dict[str, Any]:
    w = np.asarray(w, dtype=float)
    n = w.size
    if n == 0:
        return {}
    abs_w = np.abs(w)
    abs_sum = abs_w.sum()
    signed_sum = w.sum()
    tawss = abs_w.mean()
    osi = 0.5 * (1.0 - (np.abs(signed_sum) / abs_sum)) if abs_sum > 0 else 0.0
    denom = (1.0 - 2.0 * osi) * tawss
    rrt = (1.0 / denom) if denom > 0 else float('inf')
    out = {"TAWSS": float(tawss), "OSI": float(osi), "RRT": float(rrt)}
    if dt is not None:
        out["TAWSS_time_integral"] = float(abs_sum * dt)
        out["net_shear_time_integral"] = float(signed_sum * dt)
        out["cycle_period"] = float(n * dt)
    return out

def _read_scalar_exelem(filepath):
    out = {}
    if not os.path.exists(filepath):
        return out
    with open(filepath, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Element:"):
            m = re.search(r"Element:\s+(\d+)", lines[i])
            if m and i + 2 < len(lines) and lines[i + 1].strip() == "Values:":
                eid = int(m.group(1))
                val = safe_float(lines[i + 2].strip().split()[0])
                if val is not None:
                    out[eid] = val
                i += 3
            else:
                i += 1
        else:
            i += 1
    return out

# ============================================================
# ---------- Build base topology + steady placeholders -------
# ============================================================

def build_master_elements(
    subject_path: str,
    element_file_name: str,
    node_file_name: str,
    flow_file_name: str,
    radius_file_name: str,
    flow_radius_dir: str,
) -> Dict[int, Dict[str, Any]]:
    # node coordinates (optional)
    node_coords = {}
    node_path = os.path.join(subject_path, 'Vessel', node_file_name)
    if os.path.exists(node_path):
        with open(node_path, "r") as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.strip().startswith("Node:"):
                m = re.search(r"Node:\s+(\d+)", line)
                if m:
                    nid = int(m.group(1))
                    if i + 3 < len(lines):
                        x = safe_float(lines[i + 1].strip())
                        y = safe_float(lines[i + 2].strip())
                        z = safe_float(lines[i + 3].strip())
                        if None not in (x, y, z):
                            node_coords[nid] = (x, y, z)
                        i += 4
                    else:
                        i += 1
                else:
                    i += 1
            else:
                i += 1

    # connectivity
    elements: Dict[int, Dict[str, Any]] = {}
    conn_path = os.path.join(subject_path, 'Vessel', element_file_name)
    with open(conn_path, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("Element:"):
            m = re.search(r"Element:\s+(\d+)", line)
            if m:
                eid = int(m.group(1))
                if i + 1 < len(lines) and lines[i + 1].strip() == "Nodes:":
                    node_line = lines[i + 2].strip() if i + 2 < len(lines) else ""
                    node_ids = [int(tok) for tok in node_line.split() if re.match(r"^-?\d+$", tok)]
                    if len(node_ids) >= 2:
                        start_node, end_node = node_ids[0], node_ids[-1]
                        elements[eid] = {
                            "start_node": start_node,
                            "end_node":   end_node,
                            "start_coord": node_coords.get(start_node),
                            "end_coord":   node_coords.get(end_node),
                            "generation": None,
                            "strahler_order": None,
                            "horsfield_order": None,
                            "flow_mm3_s": None,
                            "radius_mm":  None,
                            "steady_WSS_Poiseuille": None,
                        }
                    i += 3
                else:
                    i += 1
        else:
            i += 1

    # read steady scalars from provided directory
    def update_scalar_from_exelem(filepath, field_name):
        if not os.path.exists(filepath):
            return
        with open(filepath, "r") as f:
            lines = f.readlines()
        j = 0
        while j < len(lines):
            if lines[j].strip().startswith("Element:"):
                m = re.search(r"Element:\s+(\d+)", lines[j])
                if m and j + 2 < len(lines) and lines[j + 1].strip() == "Values:":
                    eid = int(m.group(1))
                    val = safe_float(lines[j + 2].strip().split()[0])
                    if eid in elements:
                        elements[eid][field_name] = val
                    j += 3
                else:
                    j += 1
            else:
                j += 1

    update_scalar_from_exelem(os.path.join(flow_radius_dir, flow_file_name), "flow_mm3_s")
    update_scalar_from_exelem(os.path.join(flow_radius_dir, radius_file_name), "radius_mm")

    # steady WSS (Poiseuille) baseline fill
    for eid, d in elements.items():
        Q_mm3_s = d["flow_mm3_s"]
        R_mm    = d["radius_mm"]
        if Q_mm3_s is None or R_mm is None or R_mm == 0:
            d["steady_WSS_Poiseuille"] = None
            continue
        Q = Q_mm3_s * MM3_TO_M3
        R = R_mm * MM_TO_M
        try:
            d["steady_WSS_Poiseuille"] = (4.0 * MU * Q) / (math.pi * (R**3))
        except ZeroDivisionError:
            d["steady_WSS_Poiseuille"] = None

    # children
    elem_by_start_node = defaultdict(list)
    for eid, d in elements.items():
        elem_by_start_node[d["start_node"]].append(eid)

    children = defaultdict(list)
    for eid, d in elements.items():
        children[eid] = elem_by_start_node.get(d["end_node"], [])

    # roots
    end_nodes_present = {d["end_node"] for d in elements.values()}
    root_start_nodes = {d["start_node"] for d in elements.values()} - end_nodes_present

    # generation
    visited = set()
    def assign_generation(eid, gen):
        if eid in visited:
            return
        visited.add(eid)
        elements[eid]["generation"] = gen
        kids = children[eid]
        next_gen = gen + 1 if len(kids) > 1 else gen
        for k in kids:
            assign_generation(k, next_gen)

    for node in root_start_nodes:
        for root_eid in elem_by_start_node.get(node, []):
            assign_generation(root_eid, 1)

    # Strahler/Horsfield
    str_cache = {}
    def compute_strahler(eid):
        if eid in str_cache:
            return str_cache[eid]
        kids = children[eid]
        if not kids:
            str_cache[eid] = 1
            return 1
        kid_orders = [compute_strahler(k) for k in kids]
        m = max(kid_orders)
        order = m + 1 if kid_orders.count(m) >= 2 else m
        str_cache[eid] = order
        return order

    hf_cache = {}
    def compute_horsfield(eid):
        if eid in hf_cache:
            return hf_cache[eid]
        kids = children[eid]
        if not kids:
            hf_cache[eid] = 1
            return 1
        order = max(compute_horsfield(k) for k in kids) + 1
        hf_cache[eid] = order
        return order

    for eid in elements.keys():
        elements[eid]["strahler_order"]  = compute_strahler(eid)
        elements[eid]["horsfield_order"] = compute_horsfield(eid)

    return elements

# ============================================================
# --------- Scenario loaders & simple summarizers ------------
# ============================================================

def compute_steady_elements_for_dir(
    elements_base: Dict[int, Dict[str, Any]],
    scenario_dir: str,
    flow_file_name: str,
    radius_file_name: str,
    field_name: str = "steady_WSS_Poiseuille",
) -> Dict[int, Dict[str, Any]]:
    flow_map   = _read_scalar_exelem(os.path.join(scenario_dir, flow_file_name))
    radius_map = _read_scalar_exelem(os.path.join(scenario_dir, radius_file_name))

    elems = {}
    for eid, d in elements_base.items():
        e = dict(d)
        if eid in flow_map:
            e["flow_mm3_s"] = flow_map[eid]
        if eid in radius_map:
            e["radius_mm"] = radius_map[eid]
        Q_mm3_s = e.get("flow_mm3_s")
        R_mm    = e.get("radius_mm")
        if Q_mm3_s is None or R_mm is None or R_mm == 0:
            e[field_name] = None
        else:
            Q = Q_mm3_s * MM3_TO_M3
            R = R_mm * MM_TO_M
            try:
                e[field_name] = (4.0 * MU * Q) / (math.pi * (R**3))
            except ZeroDivisionError:
                e[field_name] = None
        elems[eid] = e
    return elems

def summarize_mean_by_key(elements: Dict[int, Dict[str, Any]], key="strahler_order", field="steady_WSS_Poiseuille"):
    groups: Dict[Any, List[float]] = defaultdict(list)
    for _, d in elements.items():
        g = d.get(key, None)
        v = d.get(field, None)
        if g is None or v is None or not np.isfinite(v):
            continue
        groups[g].append(float(v))
    out = {}
    for g, vals in groups.items():
        arr = np.asarray(vals, dtype=float)
        n = arr.size
        if n == 0:
            continue
        mean = float(arr.mean())
        std  = float(arr.std(ddof=1)) if n > 1 else 0.0
        out[g] = {"mean": mean, "std": std, "n": int(n)}
    return out

def parse_dynamic_and_metrics_for_scenario(
    elements_base: Dict[int, Dict[str, Any]],
    scenario_dir: str,
    num_timepoints: Optional[int] = DEFAULT_NUM_TIMEPOINTS,
    store_series: bool = DEFAULT_STORE_WSS_SERIES
) -> Dict[int, Dict[str, Any]]:
    scen = {
        eid: {
            "strahler_order": d.get("strahler_order"),
            "horsfield_order": d.get("horsfield_order"),
        } for eid, d in elements_base.items()
    }
    wss_path = os.path.join(scenario_dir, wss_filename_in_scenario)
    series, N = parse_wss_series_file(wss_path, num_timepoints=num_timepoints)
    f0 = extract_fundamental_frequency(os.path.join(scenario_dir, frequency_json_filename))
    dt = (1.0 / f0) / N if (f0 and N) else None

    for eid, w in series.items():
        if N and w.size != N:
            continue
        if store_series:
            scen[eid]["wss_series"] = [float(x) for x in w.tolist()]
        metrics = compute_wss_metrics_from_series(w, dt=dt)
        scen[eid].update(metrics)
    return scen

def summarize_dynamic_by_order(scen: Dict[int, Dict[str, Any]], order_key: str, metric: str) -> Dict[Any, Dict[str, float]]:
    groups: Dict[Any, List[float]] = defaultdict(list)
    for _, d in scen.items():
        g = d.get(order_key)
        v = d.get(metric)
        if g is None or v is None or not np.isfinite(v):
            continue
        groups[g].append(float(v))
    out = {}
    for g, vals in groups.items():
        arr = np.asarray(vals, float)
        if arr.size == 0:
            continue
        q1, med, q3 = np.percentile(arr, [25, 50, 75])
        out[g] = {"median": float(med), "IQR": float(q3 - q1), "n": int(arr.size)}
    return out

# ============================================================
# ----------------- Order → element lists --------------------
# ============================================================

def elements_by_order(elements_like: Dict[int, Dict[str, Any]], key: str) -> Dict[int, List[int]]:
    """
    Return {order: [element_ids]} for the given order key ('strahler_order' or 'horsfield_order').
    Orders and element IDs are sorted for determinism.
    """
    buckets: Dict[int, List[int]] = defaultdict(list)
    for eid, d in elements_like.items():
        o = d.get(key, None)
        if o is None:
            continue
        try:
            o = int(o)
        except Exception:
            pass
        buckets[o].append(int(eid))
    out = {}
    for o in sorted(buckets.keys(), key=lambda x: int(x) if str(x).isdigit() else x):
        out[o] = sorted(buckets[o])
    return out

def get_orders_maps(elements_base: Dict[int, Dict[str, Any]]) -> Dict[str, Dict[int, List[int]]]:
    """
    Convenience: returns {"strahler": {...}, "horsfield": {...}}
    """
    return {
        "strahler":  elements_by_order(elements_base, "strahler_order"),
        "horsfield": elements_by_order(elements_base, "horsfield_order"),
    }

# ============================================================
# ----------------------- Plotting ---------------------------
# ============================================================

def _sorted_orders_union(summary_dicts: List[Dict[Any, Any]]) -> List[Any]:
    orders = set()
    for sd in summary_dicts:
        orders |= set(sd.keys())
    def _coerce(k):
        try: return int(k)
        except: return k
    try:
        return sorted(orders, key=lambda k: _coerce(k))
    except TypeError:
        return list(orders)

def _overlay_errorbar(
    summary_dicts: List[Dict[Any, Any]],
    labels: List[str],
    y_key_mean: str,
    y_key_err: str,
    xlabel: str,
    ylabel: str,
    title: str
):
    if not summary_dicts:
        return None, None
    orders = _sorted_orders_union(summary_dicts)
    if not orders:
        return None, None
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for sd, lab in zip(summary_dicts, labels):
        xs, ys, es = [], [], []
        for o in orders:
            if o in sd and y_key_mean in sd[o]:
                m = sd[o].get(y_key_mean)
                e = sd[o].get(y_key_err, 0.0)
                if m is None:
                    continue
                xval = int(o) if str(o).isdigit() else o
                xs.append(xval); ys.append(m); es.append(e)
        if xs:
            ax.errorbar(xs, ys, yerr=es, fmt='-o', label=lab)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    return fig, ax

# ============================================================
# ---------------------- Main runner -------------------------
# ============================================================

def run_analysis(
    subject_id: str,
    subject_path: str,
    scenario_dirs: List[str],
    labels: List[str],
    element_file_name: str,
    node_file_name: str,
    flow_file_name: str,
    radius_file_name: str,
    num_timepoints: Optional[int] = DEFAULT_NUM_TIMEPOINTS,
    store_series: bool = DEFAULT_STORE_WSS_SERIES,
    save_dir: Optional[str] = None,
    show: bool = True,
):
    """
    Returns a dict with:
      - steady_by_strahler, steady_by_horsfield
      - tawss_by_strahler, tawss_by_horsfield
      - osi_by_strahler,   osi_by_horsfield
      - orders: {"strahler": {order:[eids]}, "horsfield": {order:[eids]}}
    """
    assert len(scenario_dirs) == len(labels), "scenario_dirs and labels must have same length"

    # Build base topology from the first scenario's steady scalars (only used as initial fill)
    base_dir_for_scalars = scenario_dirs[0]
    elements_base = build_master_elements(
        subject_path=subject_path,
        element_file_name=element_file_name,
        node_file_name=node_file_name,
        flow_file_name=flow_file_name,
        radius_file_name=radius_file_name,
        flow_radius_dir=base_dir_for_scalars,
    )

    # Expose element lists by order
    orders_map = get_orders_maps(elements_base)
    # Access later in IDE:
    #   results = run_analysis(...)
    #   results["orders"]["strahler"][3]  -> list of element IDs with Strahler=3
    #   results["orders"]["horsfield"][2] -> list of element IDs with Horsfield=2

    # -------- Steady-state summaries (mean±SD) per scenario
    steady_by_st_list, steady_by_hf_list = [], []
    for d in scenario_dirs:
        elems = compute_steady_elements_for_dir(
            elements_base, d, flow_file_name, radius_file_name, field_name="steady_WSS_Poiseuille"
        )
        steady_by_st_list.append(summarize_mean_by_key(elems, key="strahler_order",  field="steady_WSS_Poiseuille"))
        steady_by_hf_list.append(summarize_mean_by_key(elems, key="horsfield_order", field="steady_WSS_Poiseuille"))

    # -------- Dynamic (TAWSS/OSI) summaries per scenario (median±IQR/2)
    dyn_by_st_T, dyn_by_hf_T = [], []
    dyn_by_st_O, dyn_by_hf_O = [], []
    for d in scenario_dirs:
        scen = parse_dynamic_and_metrics_for_scenario(
            elements_base, d, num_timepoints=num_timepoints, store_series=store_series
        )
        dyn_by_st_T.append(summarize_dynamic_by_order(scen, "strahler_order",  "TAWSS"))
        dyn_by_hf_T.append(summarize_dynamic_by_order(scen, "horsfield_order", "TAWSS"))
        dyn_by_st_O.append(summarize_dynamic_by_order(scen, "strahler_order",  "OSI"))
        dyn_by_hf_O.append(summarize_dynamic_by_order(scen, "horsfield_order", "OSI"))

    # ===== Figures =====
    titles_suffix = f" — scenarios for {subject_id}"

    # 1) Steady WSS vs Strahler (mean±SD)
    fig1, _ = _overlay_errorbar(
        steady_by_st_list, labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Strahler order", ylabel="Steady WSS (Pa)",
        title=f"Steady WSS by Strahler (mean ± SD){titles_suffix}"
    )
    if save_dir and fig1:
        os.makedirs(save_dir, exist_ok=True)
        fig1.savefig(os.path.join(save_dir, "steady_wss_by_strahler_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig1: plt.show()
    elif fig1: plt.close(fig1)

    # 2) Steady WSS vs Horsfield (mean±SD)
    fig2, _ = _overlay_errorbar(
        steady_by_hf_list, labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Horsfield order", ylabel="Steady WSS (Pa)",
        title=f"Steady WSS by Horsfield (mean ± SD){titles_suffix}"
    )
    if save_dir and fig2:
        fig2.savefig(os.path.join(save_dir, "steady_wss_by_horsfield_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig2: plt.show()
    elif fig2: plt.close(fig2)

    # Helper to convert median/IQR -> error=IQR/2
    def convert_med_iqr(sd: Dict[Any, Dict[str, float]]) -> Dict[Any, Dict[str, float]]:
        out = {}
        for k, v in sd.items():
            med = v.get("median")
            iqr = v.get("IQR", 0.0)
            if med is not None:
                out[k] = {"mean": float(med), "std": float(iqr) / 2.0}
        return out

    # 3) TAWSS vs Strahler (median±IQR/2)
    fig3, _ = _overlay_errorbar(
        [convert_med_iqr(sd) for sd in dyn_by_st_T], labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Strahler order", ylabel="TAWSS (Pa)",
        title=f"TAWSS by Strahler (median ± IQR/2){titles_suffix}"
    )
    if save_dir and fig3:
        fig3.savefig(os.path.join(save_dir, "tawss_by_strahler_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig3: plt.show()
    elif fig3: plt.close(fig3)

    # 4) TAWSS vs Horsfield (median±IQR/2)
    fig4, _ = _overlay_errorbar(
        [convert_med_iqr(sd) for sd in dyn_by_hf_T], labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Horsfield order", ylabel="TAWSS (Pa)",
        title=f"TAWSS by Horsfield (median ± IQR/2){titles_suffix}"
    )
    if save_dir and fig4:
        fig4.savefig(os.path.join(save_dir, "tawss_by_horsfield_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig4: plt.show()
    elif fig4: plt.close(fig4)

    # 5) OSI vs Strahler (median±IQR/2)
    fig5, _ = _overlay_errorbar(
        [convert_med_iqr(sd) for sd in dyn_by_st_O], labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Strahler order", ylabel="OSI",
        title=f"OSI by Strahler (median ± IQR/2){titles_suffix}"
    )
    if save_dir and fig5:
        fig5.savefig(os.path.join(save_dir, "osi_by_strahler_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig5: plt.show()
    elif fig5: plt.close(fig5)

    # 6) OSI vs Horsfield (median±IQR/2)
    fig6, _ = _overlay_errorbar(
        [convert_med_iqr(sd) for sd in dyn_by_hf_O], labels,
        y_key_mean="mean", y_key_err="std",
        xlabel="Horsfield order", ylabel="OSI",
        title=f"OSI by Horsfield (median ± IQR/2){titles_suffix}"
    )
    if save_dir and fig6:
        fig6.savefig(os.path.join(save_dir, "osi_by_horsfield_multi.png"), dpi=150, bbox_inches="tight")
    if show and fig6: plt.show()
    elif fig6: plt.close(fig6)

    return {
        "steady_by_strahler": steady_by_st_list,
        "steady_by_horsfield": steady_by_hf_list,
        "tawss_by_strahler": dyn_by_st_T,
        "tawss_by_horsfield": dyn_by_hf_T,
        "osi_by_strahler": dyn_by_st_O,
        "osi_by_horsfield": dyn_by_hf_O,
        "orders": orders_map,  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # Access example:
        #   results["orders"]["strahler"][3]  -> list of element IDs with Strahler=3
        #   results["orders"]["horsfield"][2] -> list of element IDs with Horsfield=2
    }

# ============================================================
# ------------------- Example invocation ---------------------
# ============================================================

if __name__ == "__main__":
    perfusion_folder = os.path.join(subject_path, 'Perfusion')
    scenario_dirs = [
        os.path.join(perfusion_folder, "WSS_test0.0"),
        os.path.join(perfusion_folder, "WSS_test_PE0.0"),
        os.path.join(perfusion_folder, "WSS_test_37.0"),
        os.path.join(perfusion_folder, "WSS_test_adj37.0"),
    ]
    labels = ["Baseline", "PE", "CTEPH", "Adjusted"]

    results = run_analysis(
        subject_id=subject,
        subject_path=subject_path,
        scenario_dirs=scenario_dirs,
        labels=labels,
        element_file_name=element_file_name,
        node_file_name=node_file_name,
        flow_file_name=flow_file_name,
        radius_file_name=radius_file_name,
        num_timepoints=DEFAULT_NUM_TIMEPOINTS,
        store_series=False,
        save_dir=None,   # or a path to save PNGs
        show=True
    )

    # Example: access orders in console:
    # print(results["orders"]["strahler"])
    # print(results["orders"]["horsfield"])

# ============================================================
# --------- (Optional) downstream/seed analysis (OFF) --------
# ============================================================

# The following functions are kept for future use but intentionally
# commented-out so they won't run or plot right now.

# def _build_children(elements_base, allowed=None):
#     if allowed is None:
#         allowed = set(elements_base.keys())
#     else:
#         allowed = set(allowed)
#     start_map = defaultdict(list)
#     for eid, d in elements_base.items():
#         if eid not in allowed:
#             continue
#         sn = d.get("start_node"); en = d.get("end_node")
#         if sn is None or en is None:
#             continue
#         start_map[sn].append(eid)
#     children = {eid: [] for eid in allowed}
#     for eid in allowed:
#         en = elements_base[eid].get("end_node")
#         if en is None:
#             continue
#         kids = [k for k in start_map.get(en, []) if k in allowed and k != eid]
#         children[eid] = kids
#     return children
#
# def _downstream_set(seed, children, include_self=False):
#     seen = set()
#     stack = list(children.get(seed, []))
#     while stack:
#         cur = stack.pop()
#         if cur in seen:
#             continue
#         seen.add(cur)
#         stack.extend(children.get(cur, []))
#     if include_self:
#         seen.add(seed)
#     return seen
#
# def _mean_std(values):
#     arr = np.array([v for v in values if v is not None and np.isfinite(v)], float)
#     n = arr.size
#     if n == 0:
#         return None, None, 0
#     if n == 1:
#         return float(arr.mean()), 0.0, 1
#     return float(arr.mean()), float(arr.std(ddof=1)), int(n)
#
# def compare_downstream_means(elements_base, scenario_dicts, labels, seeds, allowed=None, metrics=("TAWSS","OSI"), include_self=False, decimals=4):
#     children = _build_children(elements_base, allowed=allowed)
#     allowed_set = set(allowed) if allowed is not None else set(elements_base.keys())
#     print("\n=== Downstream averages (mean ± SD) per seed, across scenarios ===\n")
#     for seed in seeds:
#         if seed not in elements_base:
#             print(f"Seed {seed}: not in elements_base → skip\n")
#             continue
#         if seed not in allowed_set:
#             print(f"Seed {seed}: not in allowed set → skip\n")
#             continue
#         ds = _downstream_set(seed, children, include_self=include_self)
#         if not ds:
#             print(f"Seed {seed}: no downstream elements within the current graph.\n")
#             continue
#         print(f"-- Seed {seed} -- downstream elements: {len(ds)}")
#         header = "Scenario".ljust(18)
#         for m in metrics:
#             header += f" | {m} mean±SD (n)".ljust(28)
#         print(header)
#         print("-" * len(header))
#         for scen, lab in zip(scenario_dicts, labels):
#             line = lab.ljust(18)
#             for m in metrics:
#                 vals = [scen.get(e, {}).get(m) for e in ds]
#                 mu, sd, n = _mean_std(vals)
#                 if n == 0:
#                     cell = "—".ljust(26)
#                 else:
#                     cell = f"{mu:.{decimals}f}±{sd:.{decimals}f} (n={n})".ljust(26)
#                 line += " | " + cell
#             print(line)
#         print()
#     print("=== End ===\n")
#
# def plot_downstream_by_strahler_for_seed(elements_base, scenario_dicts, labels, seed, allowed=None, include_self=False, metrics=("TAWSS","OSI"), save_dir=None, show=True):
#     children = _build_children(elements_base, allowed=allowed)
#     ds = _downstream_set(seed, children, include_self=include_self)
#     if not ds:
#         print(f"[plot] Seed {seed}: no downstream elements under current settings.")
#         return
#     def _group_stats_by_order(elements_base, scen, eids, metric="TAWSS", order_key="strahler_order"):
#         buckets = {}
#         for eid in eids:
#             d = elements_base.get(eid, {})
#             order = d.get(order_key)
#             if order is None:
#                 continue
#             val = scen.get(eid, {}).get(metric)
#             if val is None or not np.isfinite(val):
#                 continue
#             buckets.setdefault(order, []).append(float(val))
#         if not buckets:
#             return [], [], [], []
#         orders = sorted(buckets.keys(), key=lambda x: int(x) if str(x).isdigit() else x)
#         means, sds, ns = [], [], []
#         for o in orders:
#             arr = np.asarray(buckets[o], dtype=float)
#             n = arr.size
#             ns.append(int(n))
#             means.append(float(arr.mean()) if n > 0 else np.nan)
#             sds.append(float(arr.std(ddof=1)) if n > 1 else 0.0)
#         return orders, means, sds, ns
#     for metric in metrics:
#         fig, ax = plt.subplots(figsize=(7, 4.5))
#         all_orders = set()
#         per_scen_data = []
#         for scen in scenario_dicts:
#             orders, means, sds, ns = _group_stats_by_order(elements_base, scen, ds, metric=metric)
#             per_scen_data.append((orders, means, sds, ns))
#             all_orders |= set(orders)
#         if not all_orders:
#             plt.close(fig)
#             continue
#         orders_sorted = sorted(all_orders, key=lambda x: int(x) if str(x).isdigit() else x)
#         for (orders, means, sds, ns), lab in zip(per_scen_data, labels):
#             mean_map = {o: m for o, m in zip(orders, means)}
#             sd_map   = {o: s for o, s in zip(orders, sds)}
#             y  = [mean_map.get(o, np.nan) for o in orders_sorted]
#             ey = [sd_map.get(o,   0.0)    for o in orders_sorted]
#             xs = [int(o) if str(o).isdigit() else o for o in orders_sorted]
#             xs_plot, y_plot, ey_plot = [], [], []
#             for xi, yi, eyi in zip(xs, y, ey):
#                 if np.isfinite(yi):
#                     xs_plot.append(xi); y_plot.append(yi); ey_plot.append(eyi)
#             if xs_plot:
#                 ax.errorbar(xs_plot, y_plot, yerr=ey_plot, fmt='-o', label=lab)
#         ax.set_xlabel("Strahler order")
#         ax.set_ylabel(metric)
#         ax.set_title(f"{metric} downstream of seed {seed} (mean ± SD by order)")
#         ax.grid(True, linestyle='--', alpha=0.5)
#         ax.legend()
#         fig.tight_layout()
#         if save_dir:
#             import os
#             os.makedirs(save_dir, exist_ok=True)
#             out = os.path.join(save_dir, f"seed_{seed}_{metric}_by_strahler.png")
#             fig.savefig(out, dpi=150, bbox_inches="tight")
#         if show:
#             plt.show()
#         else:
#             plt.close(fig)
