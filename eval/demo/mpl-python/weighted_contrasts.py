"""Python port of contr.poly.weighted() from src/weightedcontrasts.R."""

import numpy as np
import pandas as pd


def poly_weighted(levels, counts, width=1):
    """Polynomial contrasts weighted by group size.

    levels: ordered level names; counts: number of samples per level.
    Returns a DataFrame indexed by level with columns L, Q, ... (n - 1 columns).
    """
    n = len(levels)
    if n < 2:
        raise ValueError("contrasts need at least two levels")
    w = np.asarray(counts, dtype=float)
    w = w / w.sum()
    y = np.arange(1, n + 1) - np.sum(np.arange(1, n + 1) * w)
    X = np.sqrt(w)[:, None] * y[:, None] ** np.arange(n)
    q, r = np.linalg.qr(X)
    raw = q * np.diag(r) / np.sqrt(w)[:, None]
    contr = raw / np.sqrt((raw**2).sum(axis=0))
    scores = np.arange(1, width * n + 1, width, dtype=float)
    scores = scores - np.sum(scores * w)
    contr[:, 1] *= np.sqrt(np.sum(scores**2))
    names = ["L", "Q", "C"][: n - 1] + [f"^{k}" for k in range(4, n)]
    return pd.DataFrame(contr[:, 1:], index=levels, columns=names)
