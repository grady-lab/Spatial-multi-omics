"""Centralized path resolution.

Read `config.yaml` (and optionally `config.local.yaml`) and expose paths as a
nested namespace so scripts and notebooks never hardcode locations.

Usage:
    from paths import P
    sc.read_h5ad(P.processed.adata.all_cells / P.fn.combined_adata.format(key="A1"))
    raw_dir = P.slides["A1"]   # absolute, from config.local.yaml

Override the project root by setting GDF15_ROOT (useful on HPC).
"""
from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import yaml

def _find_root() -> Path:
    """Locate the project root.

    Order of preference:
      1. `GDF15_ROOT` env var (useful on HPC).
      2. The directory containing this file's parent — i.e., the project root
         when paths.py is imported as a regular module from a script.
      3. Walk up from the current working directory looking for `config.yaml`
         — used when imported from a notebook (no `__file__`-relative anchor
         is reliable across Jupyter kernels).
    """
    env = os.environ.get("GDF15_ROOT")
    if env:
        return Path(env)
    try:
        return Path(__file__).resolve().parents[1]
    except NameError:
        cur = Path.cwd().resolve()
        while not (cur / "config.yaml").exists() and cur != cur.parent:
            cur = cur.parent
        return cur


ROOT = _find_root()

with open(ROOT / "config.yaml") as f:
    _cfg = yaml.safe_load(f)

_local_path = ROOT / "config.local.yaml"
_local = yaml.safe_load(_local_path.read_text()) if _local_path.exists() else {}


def _to_paths(node):
    if isinstance(node, dict):
        return SimpleNamespace(**{k: _to_paths(v) for k, v in node.items()})
    return ROOT / node


P = _to_paths(_cfg["paths"])
P.fn = SimpleNamespace(**_cfg["filenames"])
P.slides = _local.get("slides", {})


def ensure_dirs() -> None:
    """Create the full output tree. Call once at the top of pipeline scripts."""
    for leaf in [
        P.metadata, P.annotations,
        P.interim.adata, P.interim.sdata,
        P.processed.adata.all_cells,
        P.processed.adata.nonepithelial,
        P.processed.adata.immune, P.processed.adata.stromal,
        P.processed.adata.epithelial, P.processed.adata.tcells,
        P.results.figures, P.results.tables, P.results.xenium_explorer,
    ]:
        leaf.mkdir(parents=True, exist_ok=True)
