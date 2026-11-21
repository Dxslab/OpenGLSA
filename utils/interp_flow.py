#!/usr/bin/env python3
# Copyright (c) 2026 Ruiqing Yang and Baofeng Ma
# SPDX-License-Identifier: Apache-2.0

"""Interpolate a 2D OpenGLSA flow field from an old mesh onto a new mesh.

File formats (Fortran sequential unformatted, see utils/x2op.f90 and 2D/IO.f90):

  grid.dat :  record 1 : Axx(nx_global, ny_global)   real(8)
              record 2 : Ayy(nx_global, ny_global)   real(8)

  flow.dat :  record 1 : Istep (int4), tt (real8)
              record 2 : d(nx_global, ny_global)     real(8)
              record 3 : u(...)
              record 4 : v(...)
              record 5 : T(...)

Arrays are stored in Fortran order (i fastest).  The grid files carry no
dimension header, so the shapes are taken from control.nml / the command line,
or guessed from the record length when neither is available.

Interpolation is piecewise linear over the old mesh, triangulated cell by cell
(each quadrilateral split along its shorter diagonal).  Using the mesh's own
cells rather than a Delaunay triangulation of the point cloud matters here:
wall-normal stretching is extreme, and Delaunay would happily connect points
across several cells.  New points that fall outside the old mesh -- boundary
nodes that miss by a rounding error, or a genuinely larger new domain -- are
projected onto the closest triangle and evaluated there; the largest projection
distance is reported so the extrapolation stays visible.

Usage (from the directory holding the files):

    python3 interp_flow.py                      # uses control.nml for the new shape
    python3 interp_flow.py --old-shape 1121x241 --new-shape 1641x321
    python3 interp_flow.py --method cubic
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
from matplotlib.tri import (
    CubicTriInterpolator,
    LinearTriInterpolator,
    Triangulation,
)
from scipy.io import FortranFile
from scipy.spatial import cKDTree

VARS = ("d", "u", "v", "T")


# --------------------------------------------------------------------------
# shapes
# --------------------------------------------------------------------------
def parse_shape(text):
    """Parse '1641x321' / '1641,321' / '1641 321' into (nx, ny)."""
    parts = re.split(r"[x,X*\s]+", text.strip())
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"cannot parse shape {text!r}, expected NXxNY")
    return int(parts[0]), int(parts[1])


def shape_from_control(path):
    """Read nx_global / ny_global from a control.nml, or return None."""
    if not os.path.isfile(path):
        return None
    with open(path, "r", errors="ignore") as fh:
        text = fh.read()
    got = {}
    for key in ("nx_global", "ny_global"):
        m = re.search(rf"^\s*{key}\s*=\s*(\d+)", text, re.MULTILINE | re.IGNORECASE)
        if m:
            got[key] = int(m.group(1))
    if len(got) == 2:
        return got["nx_global"], got["ny_global"]
    return None


def _roughness(a):
    """Relative curvature of a 2D array along both index directions.

    A structured grid reshaped with the correct nx is smooth; a wrong nx wraps
    the rows and produces large second differences.
    """
    score = 0.0
    for arr in (a, a.T):
        if arr.shape[0] < 3:
            continue
        d1 = np.abs(np.diff(arr, axis=0)).mean()
        d2 = np.abs(np.diff(arr, n=2, axis=0)).mean()
        score += d2 / (d1 + 1e-300)
    return score


def guess_shape(x, y):
    """Pick the (nx, ny) factorisation of len(x) that yields the smoothest mesh."""
    n = x.size
    best, best_score = None, np.inf
    for nx in range(2, int(n**0.5) + 1):
        if n % nx:
            continue
        for cand in ((nx, n // nx), (n // nx, nx)):
            xx = x.reshape(cand, order="F")
            yy = y.reshape(cand, order="F")
            score = _roughness(xx) + _roughness(yy)
            if score < best_score:
                best, best_score = cand, score
    if best is None:
        raise SystemExit(f"cannot factor {n} points into a 2D shape; pass --old-shape")
    return best


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def read_grid(path, shape=None, label="grid"):
    with FortranFile(path, "r") as f:
        x = f.read_record("<f8")
        y = f.read_record("<f8")
    if x.size != y.size:
        raise SystemExit(f"{path}: x and y records differ in length")
    if shape is None:
        shape = guess_shape(x, y)
        print(f"  {label}: shape guessed as {shape[0]} x {shape[1]}")
    if shape[0] * shape[1] != x.size:
        raise SystemExit(
            f"{path}: holds {x.size} points, but shape {shape[0]}x{shape[1]} "
            f"needs {shape[0] * shape[1]}"
        )
    return (x.reshape(shape, order="F"), y.reshape(shape, order="F"))


def read_flow(path, shape):
    npts = shape[0] * shape[1]
    with FortranFile(path, "r") as f:
        head = f.read_record(np.dtype([("Istep", "<i4"), ("tt", "<f8")]))
        istep, tt = int(head["Istep"][0]), float(head["tt"][0])
        fields = {}
        for name in VARS:
            rec = f.read_record("<f8")
            if rec.size != npts:
                raise SystemExit(
                    f"{path}: record '{name}' holds {rec.size} values, "
                    f"expected {npts} for a {shape[0]}x{shape[1]} mesh"
                )
            fields[name] = rec.reshape(shape, order="F")
    return istep, tt, fields


def write_flow(path, fields, istep=0, tt=0.0):
    with FortranFile(path, "w") as f:
        f.write_record(np.int32(istep), np.float64(tt))
        for name in VARS:
            # FortranFile serialises multidimensional arrays row-major, so flatten
            # in Fortran order ourselves -- otherwise the field comes out transposed.
            f.write_record(fields[name].ravel(order="F").astype("<f8", copy=False))


# --------------------------------------------------------------------------
# interpolation
# --------------------------------------------------------------------------
def build_triangulation(x, y):
    """Triangulate a structured mesh cell by cell, splitting on the shorter diagonal."""
    nx, ny = x.shape
    idx = np.arange(nx * ny).reshape((nx, ny), order="F")
    p00 = idx[:-1, :-1].ravel()
    p10 = idx[1:, :-1].ravel()
    p11 = idx[1:, 1:].ravel()
    p01 = idx[:-1, 1:].ravel()

    xf, yf = x.ravel(order="F"), y.ravel(order="F")
    d_a = (xf[p11] - xf[p00]) ** 2 + (yf[p11] - yf[p00]) ** 2  # 00-11 diagonal
    d_b = (xf[p01] - xf[p10]) ** 2 + (yf[p01] - yf[p10]) ** 2  # 10-01 diagonal
    use_a = d_a <= d_b

    tri = np.empty((2 * p00.size, 3), dtype=np.int32)
    tri[0::2] = np.where(
        use_a[:, None],
        np.column_stack((p00, p10, p11)),
        np.column_stack((p00, p10, p01)),
    )
    tri[1::2] = np.where(
        use_a[:, None],
        np.column_stack((p00, p11, p01)),
        np.column_stack((p10, p11, p01)),
    )

    # matplotlib wants counter-clockwise triangles
    ax, ay = xf[tri[:, 0]], yf[tri[:, 0]]
    bx, by = xf[tri[:, 1]], yf[tri[:, 1]]
    cx, cy = xf[tri[:, 2]], yf[tri[:, 2]]
    area2 = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
    flip = area2 < 0.0
    tri[flip] = tri[flip][:, ::-1]

    degenerate = int(np.count_nonzero(area2 == 0.0))
    if degenerate:
        print(f"  warning: {degenerate} degenerate triangles in the old mesh")

    return Triangulation(xf, yf, tri)


def _closest_point_weights(px, py, ax, ay, bx, by, cx, cy):
    """Barycentric weights of the point of triangle ABC closest to P.

    Vectorised over the triangles.  Returns (weights, squared distance); the
    weights are the true barycentric coordinates when P is inside, and those of
    the closest edge point (or vertex) when it is not.
    """
    v0x, v0y = bx - ax, by - ay
    v1x, v1y = cx - ax, cy - ay
    v2x, v2y = px - ax, py - ay

    den = v0x * v1y - v1x * v0y
    den = np.where(den == 0.0, np.finfo(float).tiny, den)
    l1 = (v2x * v1y - v1x * v2y) / den
    l2 = (v0x * v2y - v2x * v0y) / den
    l0 = 1.0 - l1 - l2
    w = np.stack((l0, l1, l2), axis=-1)

    inside = (l0 >= 0.0) & (l1 >= 0.0) & (l2 >= 0.0)
    if inside.all():
        return w, np.zeros_like(l0)

    # closest point on each edge, expressed in barycentric coordinates
    verts = np.stack((np.stack((ax, ay), -1), np.stack((bx, by), -1), np.stack((cx, cy), -1)), axis=-2)
    p = np.stack((px, py), axis=-1)
    best_d2 = np.full(l0.shape, np.inf)
    best_w = w.copy()
    for i, j in ((0, 1), (1, 2), (2, 0)):
        e = verts[..., j, :] - verts[..., i, :]
        t = np.einsum("...k,...k->...", p - verts[..., i, :], e)
        t = np.divide(t, np.einsum("...k,...k->...", e, e), out=np.zeros_like(t),
                      where=np.einsum("...k,...k->...", e, e) > 0.0)
        t = np.clip(t, 0.0, 1.0)
        q = verts[..., i, :] + t[..., None] * e
        d2 = np.einsum("...k,...k->...", p - q, p - q)
        cand = np.zeros_like(w)
        cand[..., i] = 1.0 - t
        cand[..., j] = t
        take = d2 < best_d2
        best_d2 = np.where(take, d2, best_d2)
        best_w = np.where(take[..., None], cand, best_w)

    w = np.where(inside[..., None], w, best_w)
    d2 = np.where(inside, 0.0, best_d2)
    return w, d2


def project_outside(trig, pts, k=8):
    """For points outside the mesh, find the closest triangle and its weights.

    The closest triangle is searched among those incident to the `k` nearest
    vertices, which is exact for the boundary-hugging points this is meant for.
    """
    tris = trig.triangles
    xs, ys = trig.x, trig.y

    # vertex -> incident triangles, as a CSR-style lookup
    order = np.argsort(tris.ravel(), kind="stable")
    tri_of = order // 3
    counts = np.bincount(tris.ravel(), minlength=xs.size)
    start = np.concatenate(([0], np.cumsum(counts)))

    tree = cKDTree(np.column_stack((xs, ys)))
    _, near = tree.query(pts, k=min(k, xs.size))
    near = np.atleast_2d(near)

    idx_out = np.empty(len(pts), dtype=np.int64)
    w_out = np.empty((len(pts), 3))
    d_out = np.empty(len(pts))

    for n, (p, vs) in enumerate(zip(pts, near)):
        cand = np.unique(np.concatenate([tri_of[start[v]:start[v + 1]] for v in vs]))
        t = tris[cand]
        w, d2 = _closest_point_weights(
            p[0], p[1],
            xs[t[:, 0]], ys[t[:, 0]], xs[t[:, 1]], ys[t[:, 1]], xs[t[:, 2]], ys[t[:, 2]],
        )
        b = int(np.argmin(d2))
        idx_out[n] = cand[b]
        w_out[n] = w[b]
        d_out[n] = np.sqrt(d2[b])

    return idx_out, w_out, d_out


def interpolate(old_xy, old_fields, new_xy, method="linear"):
    ox, oy = old_xy
    nx_new, ny_new = new_xy[0].shape
    xt = new_xy[0].ravel(order="F")
    yt = new_xy[1].ravel(order="F")

    print(f"  triangulating the old mesh ({(ox.shape[0] - 1) * (ox.shape[1] - 1) * 2} triangles) ...")
    trig = build_triangulation(ox, oy)

    outside = trig.get_trifinder()(xt, yt) < 0
    n_bad = int(outside.sum())
    fill = None
    if n_bad:
        idx, wts, dist = project_outside(trig, np.column_stack((xt[outside], yt[outside])))
        fill = (trig.triangles[idx], wts)
        print(
            f"  {n_bad} of {outside.size} new points ({100.0 * n_bad / outside.size:.3f}%) "
            f"lie outside the old mesh -> projected onto the nearest cell "
            f"(max distance {dist.max():.3e}, mean {dist.mean():.3e})"
        )

    out = {}
    for name in VARS:
        vals = old_fields[name].ravel(order="F")
        if method == "cubic":
            interp = CubicTriInterpolator(trig, vals, kind="geom")
        else:
            interp = LinearTriInterpolator(trig, vals)
        data = np.ma.filled(interp(xt, yt), np.nan)
        if fill is not None:
            data[outside] = np.einsum("ij,ij->i", vals[fill[0]], fill[1])
        if not np.isfinite(data).all():
            raise SystemExit(f"interpolation of '{name}' produced non-finite values")
        out[name] = data.reshape((nx_new, ny_new), order="F")

    return out, outside


# --------------------------------------------------------------------------
def report(tag, fields):
    print(f"  {tag}:")
    for name in VARS:
        a = fields[name]
        print(f"    {name}: min={a.min():14.6e}  max={a.max():14.6e}  mean={a.mean():14.6e}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Interpolate a OpenGLSA 2D flow field onto a new mesh.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--old-grid", default="grid.dat.old")
    ap.add_argument("--old-flow", default="flow.dat.old")
    ap.add_argument("--new-grid", default="grid.dat")
    ap.add_argument("--new-flow", default="flow.dat")
    ap.add_argument("--old-shape", type=parse_shape, help="NXxNY of the old mesh (guessed if omitted)")
    ap.add_argument("--new-shape", type=parse_shape, help="NXxNY of the new mesh (default: from control.nml)")
    ap.add_argument("--control", default="control.nml", help="namelist supplying nx_global/ny_global")
    ap.add_argument("--method", choices=("linear", "cubic"), default="linear")
    ap.add_argument("--istep", type=int, default=0, help="Istep written to the new flow file")
    ap.add_argument("--time", type=float, default=0.0, help="tt written to the new flow file")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output file")
    args = ap.parse_args(argv)

    for path in (args.old_grid, args.old_flow, args.new_grid):
        if not os.path.isfile(path):
            raise SystemExit(f"missing input file: {path}")
    if os.path.exists(args.new_flow) and not args.force:
        raise SystemExit(f"{args.new_flow} already exists (use --force to overwrite)")

    new_shape = args.new_shape
    if new_shape is None:
        new_shape = shape_from_control(args.control)
        if new_shape is not None:
            print(f"Reading new mesh shape from {args.control}: {new_shape[0]} x {new_shape[1]}")

    print("Reading meshes ...")
    old_xy = read_grid(args.old_grid, args.old_shape, label=args.old_grid)
    new_xy = read_grid(args.new_grid, new_shape, label=args.new_grid)
    old_shape = old_xy[0].shape
    new_shape = new_xy[0].shape
    print(f"  old mesh {old_shape[0]} x {old_shape[1]}   "
          f"x in [{old_xy[0].min():g}, {old_xy[0].max():g}]  "
          f"y in [{old_xy[1].min():g}, {old_xy[1].max():g}]")
    print(f"  new mesh {new_shape[0]} x {new_shape[1]}   "
          f"x in [{new_xy[0].min():g}, {new_xy[0].max():g}]  "
          f"y in [{new_xy[1].min():g}, {new_xy[1].max():g}]")

    print(f"Reading {args.old_flow} ...")
    istep, tt, old_fields = read_flow(args.old_flow, old_shape)
    print(f"  source field: Istep = {istep}, tt = {tt:g}")
    report("source range", old_fields)

    print(f"Interpolating ({args.method}) ...")
    new_fields, _ = interpolate(old_xy, old_fields, new_xy, method=args.method)
    report("result range", new_fields)

    for name in ("d", "T"):
        bad = int(np.count_nonzero(new_fields[name] <= 0.0))
        if bad:
            print(f"  warning: {bad} non-positive values in {name}")

    print(f"Writing {args.new_flow} (Istep = {args.istep}, tt = {args.time:g}) ...")
    write_flow(args.new_flow, new_fields, istep=args.istep, tt=args.time)
    print("Done!:>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
