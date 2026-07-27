#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs
from cluster_utils import cluster_xenium

import numpy as np
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from scipy.spatial import cKDTree
from scipy import stats as sp_stats
import scanpy as sc
import anndata as ad
import warnings; warnings.simplefilter('ignore')
import spatialdata as sd
from scipy import stats
import pathlib
import squidpy as sq
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from shapely.geometry import Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as MPLPolygon

ensure_dirs()

immune = sc.read_h5ad(str(P.processed.adata.immune / "immune_cells.h5ad"))
stromal = sc.read_h5ad(str(P.processed.adata.stromal / "stromal_cells.h5ad"))


clustered_immune_harmony = cluster_xenium(immune, 16, 50, output_path=str(P.processed.adata.immune / "immune_clustered_harmonized.h5ad"), key="immune", scale="No", hvgs="No", harmony=["batch","patient_id"],thetas=[0.5,0.5])
clustered_stromal_harmony = cluster_xenium(stromal, 16, 50, output_path=str(P.processed.adata.stromal / "stromal_clustered_harmonized.h5ad"), key="stromal", scale="No", hvgs="No", harmony=["batch","patient_id"],thetas=[0.5,0.5])
