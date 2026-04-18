import sys
from pathlib import Path
import numpy as np
import torch
import scipy.io as spio
import scipy.ndimage as ndimage
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
import os
import time
import multiprocessing as mp
from functools import partial
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
from datetime import datetime, timezone
import importlib
import skimage.io as skimio
import h5py
from scipy import stats, signal
from scipy.special import erf as _erf
import slap2_utils
import pandas as pd
import subprocess
from tqdm import tqdm
import re
from scipy.sparse.linalg import svds

sys.path.append('C:/Users/michael.xie/Documents/ophys-slap2-analysis/python')
# sys.path.append(str(Path(__file__).parent.parent))
import reconstruct

import cv2

from aind_data_schema.components.identifiers import Code
from aind_data_schema.core.processing import (
    DataProcess,
    Processing,
    ProcessName,
    ProcessStage,
    ResourceTimestamped,
    ResourceUsage,
)
from aind_data_schema_models.units import MemoryUnit
from aind_data_schema_models.system_architecture import OperatingSystem, CPUArchitecture

def to_serializable(val):
    if isinstance(val, (np.integer, np.int32, np.int64, np.uint8)):
        return int(val)
    elif isinstance(val, (np.floating, np.float32, np.float64)):
        return float(val)
    elif isinstance(val, np.ndarray):
        return val.tolist()
    return val

def nearest_interp(x, xp, yp):
    if len(xp) == 1:
        return yp
    
    x_bds = xp[:-1] / 2.0 + xp[1:] / 2.0
    idx = np.searchsorted(x_bds, x, side='left')
    idx = np.clip(idx, 0, len(xp) - 1)

    return yp[idx]

def fast_dilation(mask, kernel, iterations=1):
    if kernel is None:
        kernel = np.ones((7, 7), np.uint8)  # Structuring element for XY dilation
    # Fast path for the common 3x3 one-iteration case used in activity-map batching.
    # This avoids Python loops over (time, z) slices and is substantially faster.
    if (
        iterations == 1
        and kernel.shape == (3, 3)
        and np.all(kernel)
    ):
        m = mask.astype(bool, copy=False)
        p = np.pad(m, [(0, 0)] * (m.ndim - 2) + [(1, 1), (1, 1)], mode="constant", constant_values=False)
        out = np.zeros_like(m, dtype=bool)
        for dr in range(3):
            for dc in range(3):
                out |= p[..., dr:dr + m.shape[-2], dc:dc + m.shape[-1]]
        return out

    out = np.empty_like(mask, dtype=bool)

    # Generic path: iterate over all but last two dimensions.
    leading_shape = mask.shape[:-2]
    for idx in np.ndindex(leading_shape):
        out[idx] = cv2.dilate(
            mask[idx].astype(np.uint8, copy=False),
            kernel,
            iterations=iterations,
        ).astype(bool, copy=False)

    return out

def ref_pixs_to_drc(refPixs, dmdPixelsPerColumn, dmdPixelsPerRow):
    """Map flat reference pixel indices to DMD (depth D, column C, row R) integer indices."""
    refPixs = np.asarray(refPixs, dtype=np.int64)
    plane = int(dmdPixelsPerColumn) * int(dmdPixelsPerRow)
    npc = int(dmdPixelsPerColumn)
    refD = np.floor_divide(refPixs, plane).astype(np.int32)
    refC = np.floor_divide(refPixs - refD.astype(np.int64) * plane, npc).astype(np.int32)
    refR = np.mod(refPixs, npc).astype(np.int32)
    return refD, refC, refR

def _movmean_nan(x: np.ndarray, window: int) -> np.ndarray:
    """
    Centered moving mean that ignores NaNs, similar to MATLAB smoothdata(...,'movmean','omitnan').
    Uses partial windows at the edges (like MATLAB), returns NaN where the window has no valid data.
    """
    window = max(int(window), 1)
    if window == 1:
        return x.astype(float, copy=True)

    kernel = np.ones(window, dtype=float)
    valid = ~np.isnan(x)
    x_filled = np.where(valid, x, 0.0)

    num = np.convolve(x_filled, kernel, mode="same")
    den = np.convolve(valid.astype(float), kernel, mode="same")

    out = num / np.where(den == 0.0, np.nan, den)
    return out

def compute_f0(Fin, denoise_window: int, hull_window: int):
    """
    Python port of:

      function [F0] = computeF0(Fin, denoiseWindow, hullWindow)
      % 1) median filter to reduce noise
      % 2) rolling convex-hull-like min envelope on decimated grids
      % 3) discard doubtful samples near NaNs; smooth/fill; PCHIP back to full grid

    Args
    ----
    Fin : array_like, shape (T, ...) with time along axis 0 (NaNs allowed)
    denoise_window : int
        Median filter window (time samples).
    hull_window : int
        Window that controls the “convex hull”-like operation.

    Returns
    -------
    F0 : ndarray, same shape as Fin
    """
    F = np.asarray(Fin)
    orig_shape = F.shape
    if F.ndim == 1:
        F = F[:, None]
    else:
        F = F.reshape(F.shape[0], -1)
    
    T, C = F.shape
    if T < 4:
        return np.ones_like(Fin, dtype=float) * np.nanmean(Fin, axis=0, keepdims=True)

    hull_window = int(min(hull_window, T//4))
    delta_des = max(4.0, denoise_window / 6.0)
    
    sample_times = np.rint(np.linspace(0, T-1, num=int(np.ceil(T/delta_des)+1))).astype(int)
    n_samps_in_hull = int(np.ceil(hull_window / delta_des))

    F0 = pd.DataFrame(F).rolling(window=denoise_window, center=True, min_periods=1).median().to_numpy()
    
    origsz = F0.shape
    F0 = F0.reshape(origsz[0], -1)
    for cix in range(F0.shape[1]):
        if np.all(np.isnan(F0[:, cix])):
            continue
        
        F00 = np.full((sample_times.shape[0], n_samps_in_hull), np.nan)
        for dix in range(n_samps_in_hull, 0, -1):
            xi = sample_times[dix-1::n_samps_in_hull]
            F00[:,dix-1] = np.interp(sample_times,xi, F0[xi,cix], left=np.nan, right=np.nan)
        FF = np.nanmin(F00,axis=1)

        doubt = np.sum(~np.isnan(F00),axis=1) < int(np.ceil(n_samps_in_hull/2))
        if np.sum(~doubt)>2:
            FF[doubt] = np.nan

        win = 2*int(np.ceil(n_samps_in_hull/2.0)) + 1
        fill = _movmean_nan(FF, win)
        nan_mask = np.isnan(FF)
        FF[nan_mask] = fill[nan_mask]
        FF = _movmean_nan(FF, win)

        nan_mask = np.isnan(FF)
        if np.any(nan_mask):
            FF = np.interp(sample_times, sample_times[~nan_mask], FF[~nan_mask])

        pchip = PchipInterpolator(sample_times, FF, extrapolate=True)
        F0[:,cix] = pchip(np.arange(T))

    return F0.reshape(orig_shape)

def _gaussian_peaks_integrated(theta, yxdata):
    """Integrated isotropic 2D Gaussians over unit pixels on a regular grid.

    Port of gaussianPeaksIntegrated from getActImPeaks.m.
    theta:  (N, 4) array with columns [amp, mu_y, mu_x, sigma]
    yxdata: (M, 2) array with columns [y, x] — pixel centre coordinates
    Returns: (M,) array of predicted values.
    """
    x_int = yxdata[:, 1].astype(np.intp)
    y_int = yxdata[:, 0].astype(np.intp)
    x_min = int(x_int.min()); x_max = int(x_int.max())
    y_min = int(y_int.min()); y_max = int(y_int.max())
    x_idx = x_int - x_min
    y_idx = y_int - y_min

    A  = theta[:, 0]                                     # (N,)
    my = theta[:, 1]                                     # (N,)
    mx = theta[:, 2]                                     # (N,)
    s  = np.maximum(theta[:, 3], np.finfo(float).eps)    # (N,)

    c   = np.sqrt(np.pi / 2)
    rt2 = np.sqrt(2.0)

    xc = np.arange(x_min, x_max + 1, dtype=float)
    xL = (xc - 0.5)[:, np.newaxis]                      # (Wg, 1)
    xR = (xc + 0.5)[:, np.newaxis]                      # (Wg, 1)
    Ix = c * s[np.newaxis, :] * (
        _erf((xR - mx[np.newaxis, :]) / (rt2 * s[np.newaxis, :])) -
        _erf((xL - mx[np.newaxis, :]) / (rt2 * s[np.newaxis, :]))
    )  # (Wg, N)

    yc = np.arange(y_min, y_max + 1, dtype=float)
    yB = (yc - 0.5)[:, np.newaxis]                      # (Hg, 1)
    yT = (yc + 0.5)[:, np.newaxis]                      # (Hg, 1)
    Iy = c * s[np.newaxis, :] * (
        _erf((yT - my[np.newaxis, :]) / (rt2 * s[np.newaxis, :])) -
        _erf((yB - my[np.newaxis, :]) / (rt2 * s[np.newaxis, :]))
    )  # (Hg, N)

    img = (Iy * A[np.newaxis, :]) @ Ix.T   # (Hg, Wg)
    return img[y_idx, x_idx]


def _gaussian_peaks_integrated_val_jac(theta, yxdata):
    """Integrated Gaussians with analytical Jacobian (shared intermediates).

    theta:  (N, 4) [amp, mu_y, mu_x, sigma]
    yxdata: (M, 2) [y, x]
    Returns: (val, J) where val is (M,) and J is (M, 4*N).
    """
    x_int = yxdata[:, 1].astype(np.intp)
    y_int = yxdata[:, 0].astype(np.intp)
    M = len(x_int)
    N = theta.shape[0]

    x_min = int(x_int.min()); x_max = int(x_int.max())
    y_min = int(y_int.min()); y_max = int(y_int.max())
    x_idx = x_int - x_min
    y_idx = y_int - y_min

    A  = theta[:, 0]
    my = theta[:, 1]
    mx = theta[:, 2]
    s  = np.maximum(theta[:, 3], np.finfo(float).eps)

    c   = np.sqrt(np.pi / 2)
    rt2 = np.sqrt(2.0)
    inv_rt2s = 1.0 / (rt2 * s[np.newaxis, :])        # (1, N)

    xc = np.arange(x_min, x_max + 1, dtype=float)
    xL = (xc - 0.5)[:, np.newaxis]                   # (Wg, 1)
    xR = (xc + 0.5)[:, np.newaxis]                   # (Wg, 1)
    ux_L = (xL - mx[np.newaxis, :]) * inv_rt2s        # (Wg, N)
    ux_R = (xR - mx[np.newaxis, :]) * inv_rt2s        # (Wg, N)
    erf_ux_L = _erf(ux_L)
    erf_ux_R = _erf(ux_R)
    Ix = c * s[np.newaxis, :] * (erf_ux_R - erf_ux_L)  # (Wg, N)

    yc = np.arange(y_min, y_max + 1, dtype=float)
    yB = (yc - 0.5)[:, np.newaxis]                   # (Hg, 1)
    yT = (yc + 0.5)[:, np.newaxis]                   # (Hg, 1)
    uy_B = (yB - my[np.newaxis, :]) * inv_rt2s        # (Hg, N)
    uy_T = (yT - my[np.newaxis, :]) * inv_rt2s        # (Hg, N)
    erf_uy_B = _erf(uy_B)
    erf_uy_T = _erf(uy_T)
    Iy = c * s[np.newaxis, :] * (erf_uy_T - erf_uy_B)  # (Hg, N)

    # --- forward value ---
    img = (Iy * A[np.newaxis, :]) @ Ix.T              # (Hg, Wg)
    val = img[y_idx, x_idx]                            # (M,)

    # --- Jacobian (M, 4*N) ---
    Iy_m = Iy[y_idx, :]                               # (M, N)
    Ix_m = Ix[x_idx, :]                               # (M, N)

    # d/dA_n
    dval_dA = Iy_m * Ix_m                              # (M, N)

    # d(Iy)/d(my):  exp(-uy_B^2) - exp(-uy_T^2)   (simplification where c*2/(rt2*sqrt(pi)) = 1)
    exp_uyB2 = np.exp(-uy_B ** 2)                      # (Hg, N)
    exp_uyT2 = np.exp(-uy_T ** 2)                      # (Hg, N)
    dIy_dmy_m = (exp_uyB2 - exp_uyT2)[y_idx, :]       # (M, N)
    dval_dmy = A[np.newaxis, :] * dIy_dmy_m * Ix_m    # (M, N)

    # d(Ix)/d(mx):  exp(-ux_L^2) - exp(-ux_R^2)
    exp_uxL2 = np.exp(-ux_L ** 2)                      # (Wg, N)
    exp_uxR2 = np.exp(-ux_R ** 2)                      # (Wg, N)
    dIx_dmx_m = (exp_uxL2 - exp_uxR2)[x_idx, :]       # (M, N)
    dval_dmx = A[np.newaxis, :] * Iy_m * dIx_dmx_m    # (M, N)

    # d(Iy)/d(s) = Iy/s + sqrt(2)*(uy_B*exp(-uy_B^2) - uy_T*exp(-uy_T^2))
    # d(Ix)/d(s) = Ix/s + sqrt(2)*(ux_L*exp(-ux_L^2) - ux_R*exp(-ux_R^2))
    sqrt2 = rt2
    dIy_ds_m = (Iy / s[np.newaxis, :]
                + sqrt2 * (uy_B * exp_uyB2 - uy_T * exp_uyT2))[y_idx, :]
    dIx_ds_m = (Ix / s[np.newaxis, :]
                + sqrt2 * (ux_L * exp_uxL2 - ux_R * exp_uxR2))[x_idx, :]
    dval_ds = A[np.newaxis, :] * (dIy_ds_m * Ix_m + Iy_m * dIx_ds_m)

    J = np.empty((M, 4 * N), dtype=float)
    J[:, 0::4] = dval_dA
    J[:, 1::4] = dval_dmy
    J[:, 2::4] = dval_dmx
    J[:, 3::4] = dval_ds

    return val, J


def _lsq_curvefit(theta0, xdata, ydata, lb_flat, ub_flat, max_nfev=5000):
    """Bounded Levenberg-Marquardt solver for integrated-Gaussian curve fitting.

    Reimplements the subset of MATLAB lsqcurvefit (trust-region-reflective
    with default options) used by getActImPeaks.m, but with an analytical
    Jacobian so that each iteration needs only one forward+Jacobian evaluation
    instead of 4*N finite-difference evaluations.

    theta0:   (N, 4) initial parameters [amp, mu_y, mu_x, sigma]
    xdata:    (M, 2) pixel centre coordinates
    ydata:    (M,)   observed values (actIM(selPix) - mu_bg)
    lb_flat, ub_flat: (4*N,) bound vectors (row-major flattened)
    max_nfev: cap on function+Jacobian evaluations

    Returns: (N, 4) optimised parameters.
    """
    x = np.clip(theta0.ravel().copy(), lb_flat, ub_flat)

    val, J = _gaussian_peaks_integrated_val_jac(x.reshape(-1, 4), xdata)
    r = val - ydata
    cost = np.dot(r, r)
    nfev = 1

    lam = 1e-2          # initial damping  (MATLAB InitDamping default)

    for _ in range(400):                               # MATLAB MaxIter default
        if nfev >= max_nfev:
            break

        JtJ = J.T @ J
        Jtr = J.T @ r
        diag_JtJ = np.maximum(np.diag(JtJ), 1e-8)

        delta = np.linalg.solve(JtJ + lam * np.diag(diag_JtJ), -Jtr)
        x_new = np.clip(x + delta, lb_flat, ub_flat)

        val_new, J_new = _gaussian_peaks_integrated_val_jac(
            x_new.reshape(-1, 4), xdata)
        r_new = val_new - ydata
        cost_new = np.dot(r_new, r_new)
        nfev += 1

        if cost_new < cost:
            step_norm = np.max(np.abs(x_new - x))
            rel_cost_drop = (cost - cost_new) / max(cost, 1.0)

            x = x_new
            r = r_new
            J = J_new
            cost = cost_new
            lam = max(lam * 0.1, 1e-10)

            if step_norm < 1e-6 or rel_cost_drop < 1e-6:
                break
        else:
            lam = min(lam * 10.0, 1e10)

    return x.reshape(-1, 4)


def _detect_peaks_2d(act_im_2d, exclusion_mask, mu_bg, sigma_bg, peak_thresh, peak_th):
    """Per-plane peak detection with externally supplied background stats.

    Returns (N, 4) array [amplitude, mu_y, mu_x, sigma], or (0, 4) if none.
    """
    H, W = act_im_2d.shape
    empty = np.zeros((0, 4))
    AMP_SCALE = 1.0 / 0.75

    def _make_bounds(plocs):
        n = plocs.shape[0]
        lb = np.column_stack([np.zeros(n),
                              np.maximum(0, plocs[:, 0] - 1.5),
                              np.maximum(0, plocs[:, 1] - 1.5),
                              np.ones(n) * 0.35])
        ub = np.column_stack([np.full(n, np.inf),
                              np.minimum(H - 1, plocs[:, 0] + 1.5),
                              np.minimum(W - 1, plocs[:, 1] + 1.5),
                              np.full(n, 5.0)])
        return lb.ravel(), ub.ravel()

    def _peak_mask(tf):
        pIM = np.zeros((H, W), dtype=bool)
        if tf.shape[0] > 0:
            iy = np.clip(np.round(tf[:, 1]).astype(int), 0, H - 1)
            ix = np.clip(np.round(tf[:, 2]).astype(int), 0, W - 1)
            pIM[iy, ix] = True
        return pIM

    # --- Initial peak detection ---
    explored = act_im_2d.copy()
    explored[exclusion_mask | np.isnan(explored)] = -np.inf

    rank8 = ndimage.rank_filter(explored, rank=7, size=3)
    rank9 = ndimage.maximum_filter(explored, size=3)
    pTmp = (rank8 > peak_thresh) & (explored == rank9)

    if not np.any(pTmp):
        return empty

    pY, pX = np.where(pTmp)
    amp = act_im_2d[pY, pX] * AMP_SCALE
    n_peaks = len(pY)

    act_sel_pix = ndimage.binary_dilation(pTmp, structure=np.ones((9, 9)))
    act_sel_pix &= ~np.isnan(act_im_2d)

    thetaf = np.column_stack([amp, pY.astype(float), pX.astype(float),
                              0.5 * np.ones(n_peaks)])
    p_locs = np.column_stack([pY.astype(float), pX.astype(float)])

    # --- Initial least-squares fit ---
    sel_yx = np.column_stack(np.where(act_sel_pix)).astype(float)
    sel_vals = act_im_2d[act_sel_pix] - mu_bg

    lb_f, ub_f = _make_bounds(p_locs)
    thetaf = _lsq_curvefit(thetaf, sel_yx, sel_vals, lb_f, ub_f, max_nfev=5000)

    # --- Build fit image & residual ---
    buffer_mask = _peak_mask(thetaf)

    fit_im = np.zeros((H, W), dtype=float)
    fit_im[act_sel_pix] = _gaussian_peaks_integrated(thetaf, sel_yx)
    res_im = (act_im_2d - fit_im - mu_bg) / sigma_bg

    # --- Iterative residual peak finding (one new peak per CC per round) ---
    fit_support = fit_im > 1e-3
    reject_mask = np.zeros((H, W), dtype=bool)

    labeled_full, _ = ndimage.label(act_sel_pix)

    while True:
        # Build explored residual image
        e = res_im.copy()
        e[buffer_mask | exclusion_mask | reject_mask | ~fit_support] = -np.inf
        e[np.isnan(e)] = -np.inf

        # Find the best peak candidate in each CC simultaneously
        n_labels = int(labeled_full.max())
        if n_labels == 0:
            break

        label_ids = list(range(1, n_labels + 1))
        max_vals = ndimage.maximum(e, labeled_full, label_ids)
        max_pos = ndimage.maximum_position(e, labeled_full, label_ids)

        new_peaks = [(max_pos[i], label_ids[i])
                     for i in range(n_labels) if max_vals[i] > peak_th]
        if not new_peaks:
            break

        # Process each CC's new peak
        modified_ccs = []
        for (pY_new, pX_new), cc_label in new_peaks:
            print(f"peak of amp {act_im_2d[pY_new, pX_new]}")
            amp_new = act_im_2d[pY_new, pX_new] * AMP_SCALE

            n_before = thetaf.shape[0]
            thetaf = np.vstack([thetaf, [amp_new, float(pY_new), float(pX_new), 0.5]])
            p_locs = np.vstack([p_locs, [float(pY_new), float(pX_new)]])
            new_idx = n_before

            cc_mask = labeled_full == cc_label
            cc_yx = np.column_stack(np.where(cc_mask)).astype(float)
            cc_vals = act_im_2d[cc_mask] - mu_bg

            iy = np.clip(np.round(thetaf[:, 1]).astype(int), 0, H - 1)
            ix = np.clip(np.round(thetaf[:, 2]).astype(int), 0, W - 1)
            in_cc = cc_mask[iy, ix]

            lb_cc, ub_cc = _make_bounds(p_locs[in_cc])
            thetaf[in_cc] = _lsq_curvefit(thetaf[in_cc], cc_yx, cc_vals,
                                          lb_cc, ub_cc, max_nfev=5000)

            mu_y_new = thetaf[new_idx, 1]
            mu_x_new = thetaf[new_idx, 2]
            if abs(mu_y_new - round(mu_y_new)) < 1e-3 and abs(mu_x_new - round(mu_x_new)) < 1e-3:
                refit_mask = in_cc.copy()
                refit_mask[new_idx] = False
                if np.any(refit_mask):
                    lb_rf, ub_rf = _make_bounds(p_locs[refit_mask])
                    thetaf[refit_mask] = _lsq_curvefit(thetaf[refit_mask], cc_yx, cc_vals,
                                                       lb_rf, ub_rf, max_nfev=5000)
                thetaf = np.delete(thetaf, new_idx, axis=0)
                p_locs = np.delete(p_locs, new_idx, axis=0)
                reject_mask[pY_new, pX_new] = True

            modified_ccs.append((cc_mask, cc_yx))

        # Global update once per round (not once per peak)
        buffer_mask = _peak_mask(thetaf)
        for cc_mask, cc_yx in modified_ccs:
            fit_im[cc_mask] = _gaussian_peaks_integrated(thetaf, cc_yx)
        res_im = (act_im_2d - fit_im - mu_bg) / sigma_bg
        fit_support = fit_im > 1e-3

        # Rebuild support from current peak positions
        act_sel_pix = ndimage.binary_dilation(buffer_mask, structure=np.ones((9, 9)))
        act_sel_pix &= ~np.isnan(act_im_2d)
        labeled_full, _ = ndimage.label(act_sel_pix)

    # --- Remove small peaks (peakFuncOpt=2 threshold) ---
    if thetaf.shape[0] > 0:
        s = thetaf[:, 3]
        adj_thresh = peak_thresh / (np.pi / 2 * s**2 * _erf(1 / (np.sqrt(2) * s))**2)
        thetaf = thetaf[thetaf[:, 0] >= adj_thresh]

    return thetaf


def get_act_im_peaks(act_im, peak_th=3.0, exclusion_mask=None):
    """Find Gaussian peaks in a 3D (Z, H, W) activity image.

    Background statistics and the detection threshold are computed once
    across all planes so that sensitivity is uniform.  Each plane is then
    processed independently with the shared threshold.

    Parameters
    ----------
    act_im : (Z, H, W) array, may contain NaNs
    peak_th : float, threshold in MAD-normalised standard deviations
    exclusion_mask : None, or bool array broadcastable to (Z, H, W)

    Returns
    -------
    source_seeds : (N, 3) array with columns [z, mu_y, mu_x],
                   or (0, 3) if no peaks detected.
    """
    nZ, H, W = act_im.shape
    empty = np.zeros((0, 3))

    # Build per-plane exclusion masks
    if exclusion_mask is None:
        excl_planes = [np.zeros((H, W), dtype=bool)] * nZ
    elif exclusion_mask.ndim == 2:
        excl_planes = [exclusion_mask.astype(bool)] * nZ
    else:
        excl_planes = [exclusion_mask[z].astype(bool) for z in range(nZ)]

    # Global background statistics across all planes
    valid_vals = act_im[~np.isnan(act_im)]
    if valid_vals.size == 0:
        return empty

    mu_bg = float(np.nanmedian(act_im))
    sigma_bg = float(np.median(np.abs(valid_vals - np.median(valid_vals)))) / 0.6741891400433162
    if sigma_bg <= 0:
        return empty

    peak_thresh = mu_bg + peak_th * sigma_bg

    # Detect peaks plane-by-plane with the shared threshold
    source_seeds_list = []
    for z in range(nZ):
        thetaf_z = _detect_peaks_2d(act_im[z], excl_planes[z],
                                    mu_bg, sigma_bg, peak_thresh, peak_th)
        if thetaf_z.shape[0] > 0:
            z_col = np.full((thetaf_z.shape[0], 1), z, dtype=float)
            source_seeds_list.append(
                np.column_stack([z_col, thetaf_z[:, 1], thetaf_z[:, 2]]))

    return np.vstack(source_seeds_list) if source_seeds_list else empty


def get_trial_data(trial_info, DMDix, params, sampFreq, refStack, fastZ2RefZ, allSuperPixelIDs, dr, trialTable, all_channels=False):
    trialIx, keepTrial = trial_info
    
    if not keepTrial:
        return None

    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
    numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]

    nPixels = dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs

    source_fn = trialTable['filename'][DMDix,trialIx][0]
    firstLine = trialTable['firstLine'][DMDix,trialIx]
    lastLine = trialTable['lastLine'][DMDix,trialIx]

    importlib.reload(slap2_utils)
    if re.search(r'CYCLE\d+', source_fn):
        hDataFile = slap2_utils.MultiDataFiles(os.path.join(dr, source_fn))
    else:
        hDataFile = slap2_utils.DataFile(os.path.join(dr, source_fn))

    linesPerCycle = hDataFile.header['linesPerCycle']

    dt = 1/sampFreq/hDataFile.metaData.linePeriod_s

    DSframes = np.ceil(np.arange(firstLine, lastLine+1, dt))
    nDSframes= len(DSframes)

    # Pre-compute time windows for all frames
    dtRead = max(3 * dt, linesPerCycle)
    timeWindows = [np.arange(max(1,np.floor(DSframes[i]-dtRead)), 
                            min(np.ceil(DSframes[i]+dtRead),hDataFile.numCycles*linesPerCycle)+1) 
                    for i in range(nDSframes)]

    # Pre-compute line and cycle indices for all time windows
    lineIndices_all = [(tw - 1) % linesPerCycle + 1 for tw in timeWindows]
    cycleIndices_all = [np.floor((tw - 1) / linesPerCycle) + 1 for tw in timeWindows]

    # Collect all unique line-cycle combinations across all frames
    all_line_cycles = set()
    for DSframeIx in range(nDSframes):
        for li, ci in zip(lineIndices_all[DSframeIx], cycleIndices_all[DSframeIx]):
            all_line_cycles.add((int(li), int(ci)))

    # Get all line data at once
    all_lines = np.array([lc[0] for lc in all_line_cycles])
    all_cycles = np.array([lc[1] for lc in all_line_cycles])

    # Create a mapping from (line, cycle) to index in the cache
    line_cycle_to_idx = {(int(li), int(ci)): i for i, (li, ci) in enumerate(zip(all_lines, all_cycles))}

    start_line_data_time = time.time()
    print(f"Getting line data for {len(all_lines)} lines", end="")
    # Get all line data at once
    all_line_data = hDataFile.getLineData(all_lines, all_cycles, params['activityChannel'] if not all_channels else None)
    print(f" - completed in {time.time() - start_line_data_time:.3f} sec")

    data = np.zeros((numSuperPixels,nDSframes), dtype=np.float32)
    dataCt = np.zeros((numSuperPixels,nDSframes), dtype=np.float32)
    if all_channels:
        data2 = np.zeros((numSuperPixels,nDSframes), dtype=np.float32)
        dataCt2 = np.zeros((numSuperPixels,nDSframes), dtype=np.float32)
    # Initialize timing variables
    start_time = time.time()

    # Process each frame
    for DSframeIx in range(nDSframes):
        if DSframeIx % 100 == 0:
            avg_time = (time.time() - start_time) / DSframeIx if DSframeIx > 0 else 0
            print(f"{DSframeIx} of {nDSframes}, Average time per frame: {avg_time:.3f} sec")
        
        weights = np.exp(-np.abs(DSframes[DSframeIx] - timeWindows[DSframeIx]) / dt)

        # Get the line and cycle indices for this frame
        frame_line_indices = lineIndices_all[DSframeIx]
        frame_cycle_indices = cycleIndices_all[DSframeIx]
        
        # Process only valid lines
        valid_lines = [i for i, li in enumerate(frame_line_indices) 
                      if hDataFile.lineDataNumElements[int(li)-1] != 0]
        
        for i in valid_lines:
            line_idx = int(frame_line_indices[i])
            cycle_idx = int(frame_cycle_indices[i])
            
            # Get the cached line data
            cache_idx = line_cycle_to_idx[(line_idx, cycle_idx)]
            line_data = all_line_data[cache_idx]
            
            # Get positions and z-index
            positions = hDataFile.lineSuperPixelIDs[line_idx-1]
            zIdx = hDataFile.lineFastZIdxs[line_idx-1]
            
            # Compute lookup values and matches
            lookup_values = positions * 100 + zIdx
            matching_mask = np.isin(allSuperPixelIDs[f'DMD{DMDix+1}'], lookup_values)
            matching_indices = np.where(matching_mask)[0]
            
            if len(matching_indices) > 0:
                # Create a mapping from lookup values to their indices
                value_to_pos = dict(zip(lookup_values.astype(np.uint32), range(len(lookup_values))))
                # Get the positions in lineData for each matching index
                matched_positions = [value_to_pos[int(allSuperPixelIDs[f'DMD{DMDix+1}'][idx])] for idx in matching_indices]
                
                weight = weights[i]
                data[matching_indices, DSframeIx] += line_data[matched_positions, 0] * weight
                dataCt[matching_indices, DSframeIx] += weight
                if all_channels:
                    data2[matching_indices, DSframeIx] += line_data[matched_positions, 1] * weight
                    dataCt2[matching_indices, DSframeIx] += weight

    aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]
    # aData['DSframes'] = aData['DSframes'] * hDataFile.metaData.linePeriod_s
    
    if all_channels:
        return data/100, dataCt, aData, DSframes, data2/100, dataCt2
    else:
        return data/100, dataCt, aData, DSframes

def get_high_res_traces(trial_info, DMDix, params, sampFreq, refStack, subsampleMatrixInds, fastZ2RefZ, sparseHInds, sparseHVals, 
                allSuperPixelIDs, dr, trialTable, A_final, uniqueMotionDS, motIndsToKeepDS, median_z, psf, soma_sps):
    trialIx, keepTrial, backgroundDS = trial_info
    
    if not isinstance(A_final, torch.Tensor):
        A_final = torch.tensor(A_final, dtype=torch.float32)

    nSources = A_final.shape[1]

    if not keepTrial:
        return np.full((0,nSources),np.nan,dtype=np.float32), \
            np.full((0,nSources),np.nan,dtype=np.float32), \
                np.full((0,),np.nan), \
                    np.full((0,),np.nan), \
                        np.full((0,),np.nan), \
                            (np.full((0,),np.nan), np.full((0,),np.nan), np.full((0,),np.nan)), \
                                (np.full((0,),0,dtype=np.int16), np.full((0,),0,dtype=np.int16), np.full((0,),0,dtype=np.int16)), \
                                    np.full((0,len(soma_sps)),np.nan,dtype=np.float32)

    data_file = os.path.join(dr, f'trial_data_DMD{DMDix+1}_trial{trialIx}.npz')
    if os.path.exists(data_file):
        print(f'Loading existing trial data from {data_file}')
        data_arrays = np.load(data_file)
        dataNonNorm = data_arrays['dataNonNorm']
        dataCt = data_arrays['dataCt']
        frames = data_arrays['DSframes']
        data2NonNorm = data_arrays['data2NonNorm']
        dataCt2 = data_arrays['dataCt2']

        aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]
        aData['DSframes'] = data_arrays['aData_DSframes']
    else:
        dataNonNorm, dataCt, aData, frames, data2NonNorm, dataCt2 = get_trial_data(trial_info[:2], DMDix, params, sampFreq, refStack, fastZ2RefZ, allSuperPixelIDs, dr, trialTable, all_channels=True)
        # Save data arrays
        data_arrays = {
            'dataNonNorm': dataNonNorm,
            'dataCt': dataCt,
            'DSframes': frames,
            'aData_DSframes': aData['DSframes'],
            'data2NonNorm': data2NonNorm,
            'dataCt2': dataCt2
        }
        np.savez(data_file, **data_arrays)
        print(f'Saved trial data to {data_file}')
    
    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
    numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]

    data = dataNonNorm / dataCt
    data2 = data2NonNorm / dataCt2

    motionR = np.interp(frames, aData['DSframes'][0], aData['motionDSr'].T[0])
    motionC = np.interp(frames, aData['DSframes'][0], aData['motionDSc'].T[0])
    motionZ = np.interp(frames, aData['DSframes'][0], aData['motionDSz'].T[0])
    background = np.array([np.interp(frames, aData['DSframes'][0], backgroundDS[i]) for i in range(backgroundDS.shape[0])])

    onlineYshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineYshift'].T[0])
    onlineXshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineXshift'].T[0])
    onlineZshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineZshift'].T[0])

    uniqueMotion, motInds = np.unique(np.round(np.stack((motionR, motionC, motionZ), axis=1)), axis=0, return_inverse=True)
    
    framesToKeep = np.isin(motInds, np.flatnonzero(np.abs(uniqueMotion[:, 2] - median_z) <= 1.5))

    motInds = -1*np.ones((len(motInds),), dtype=np.int32)
    uniqueMotion, motInds[framesToKeep] = np.unique(np.round(np.stack((motionR, motionC), axis=1)[framesToKeep,:]),axis=0,return_inverse=True)

    motIndsToKeep = np.zeros_like(motIndsToKeepDS, dtype=np.int64)
    motIndsToKeep[:] = -1

    for i, motion_idx_DS in enumerate(motIndsToKeepDS):
        matches = np.flatnonzero(np.all(uniqueMotion[:,:2] == uniqueMotionDS[motion_idx_DS,:2],axis=1))
        if len(matches) > 0:
            motIndsToKeep[i] = matches[0]
    
    # background_spatial_components = background_spatial_components[:, motIndsToKeep != -1]
    motIndsToKeep = motIndsToKeep[motIndsToKeep != -1]

    framesToKeep = np.isin(motInds, motIndsToKeep)

    refD, refC, refR = ref_pixs_to_drc(subsampleMatrixInds[:, 0], dmdPixelsPerColumn, dmdPixelsPerRow)

    selPixMask = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), dtype=bool)
    for i in range(uniqueMotion.shape[0]):
        selPixMask[refD,refR + int(uniqueMotion[i,0]),refC + int(uniqueMotion[i,1])] = True
    selPixMask = ndimage.binary_dilation(selPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
    selPixIdxs = np.flatnonzero(selPixMask)

    if not isinstance(A_final, torch.Tensor):
        A_final = torch.tensor(A_final, dtype=torch.float32)

    nSources = A_final.shape[1]
    phi = torch.full((data.shape[1], nSources), np.nan, dtype=torch.float32)
    F0 = torch.full((data.shape[1], nSources), np.nan, dtype=torch.float32)

    residual = data - background
    
    for i, motion_idx in enumerate(motIndsToKeep):

        # Extract data for the most common motion mode
        motion_frames = (motInds == motion_idx).nonzero()[0]
        
        # data_tensor = torch.from_numpy(data[:, motion_frames].astype(np.float32))
        data_tensor = torch.from_numpy(residual[:, motion_frames].astype(np.float32))

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)
        sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

        # project image space (A) into superpixel space (X)
        X = torch.sparse.mm(H, A_final[selPixIdxs,:])

        # add background spatial component
        # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

        XtX = X.T @ X
        Xtd = X.T @ data_tensor  # This gives all time points at once
        
        # Add small regularization to ensure stability
        regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
        
        # Solve the system for all time points at once
        # We need to solve (X^T * X) * phi = X^T * data for each column of data
        phi[motion_frames,:] = torch.linalg.solve(
            regularized_XtX,
            Xtd
        ).T

        Xtbackground = X.T @ torch.from_numpy(background[:,motion_frames].astype(np.float32))
        F0[motion_frames,:] = torch.linalg.solve(
            XtX + 1e-10 * torch.eye(XtX.shape[0]),
            Xtbackground
        ).T

        # Xtbackground = X[:,:-1].T @ (X[:,-1].unsqueeze(-1) * phi[motion_frames,-1].unsqueeze(-1).T)
        # F0[motion_frames,:] = torch.linalg.solve(
        #     XtX[:-1,:-1] + 1e-10 * torch.eye(XtX.shape[0]-1),
        #     Xtbackground
        # ).T

    globalF = np.sum(data, axis=0)
    globalF[~framesToKeep] = np.nan

    F_soma = np.full((data.shape[1], len(soma_sps)), np.nan, dtype=np.float32)
    for i, roi_sps in enumerate(soma_sps):
        F_soma[:,i] = np.nansum(data2[roi_sps,:], axis=0)
    
    F_soma[~framesToKeep] = np.nan

    return phi.numpy(), F0.numpy(), frames, selPixIdxs, globalF, \
        (motionR, motionC, motionZ), (onlineYshifts, onlineXshifts, onlineZshifts), F_soma

def create_parameter_gui():
    root = tk.Tk()
    root.title("SLAP2 Analysis Parameters")
    
    # Create main frame
    main_frame = ttk.Frame(root, padding="10")
    main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    
    # Parameter variables
    param_vars = {}
    
    # Create parameter inputs
    row = 0
    
    # # Discard Initial (s)
    # ttk.Label(main_frame, text="Discard Initial (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    # param_vars['discardInitial_s'] = tk.DoubleVar(value=0.1)
    # ttk.Entry(main_frame, textvariable=param_vars['discardInitial_s'], width=15).grid(row=row, column=1, pady=2)
    # row += 1
    
    # Analyze Hz
    ttk.Label(main_frame, text="Analyze Hz:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['analyzeHz'] = tk.IntVar(value=100)
    ttk.Entry(main_frame, textvariable=param_vars['analyzeHz'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Glutamate Channel
    ttk.Label(main_frame, text="Glutamate Channel:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['activityChannel'] = tk.IntVar(value=1)
    ttk.Entry(main_frame, textvariable=param_vars['activityChannel'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Decay Tau (s)
    ttk.Label(main_frame, text="Decay Tau (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['decayTau_s'] = tk.DoubleVar(value=0.05)
    ttk.Entry(main_frame, textvariable=param_vars['decayTau_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Baseline Window (s)
    ttk.Label(main_frame, text="Baseline Window (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['baselineWindow_s'] = tk.DoubleVar(value=4)
    ttk.Entry(main_frame, textvariable=param_vars['baselineWindow_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Denoise Window (s)
    ttk.Label(main_frame, text="Denoise Window (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['denoiseWindow_s'] = tk.DoubleVar(value=1)
    ttk.Entry(main_frame, textvariable=param_vars['denoiseWindow_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # dXY
    ttk.Label(main_frame, text="dXY:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['dXY'] = tk.IntVar(value=5)
    ttk.Entry(main_frame, textvariable=param_vars['dXY'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Sparse Factor (log scale)
    ttk.Label(main_frame, text="Sparse Factor (log):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['sparse_fac_log'] = tk.DoubleVar(value=-3.0)
    ttk.Entry(main_frame, textvariable=param_vars['sparse_fac_log'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Operator
    ttk.Label(main_frame, text="Operator:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['operator'] = tk.StringVar(value='Maria Goeppert Mayer')
    ttk.Entry(main_frame, textvariable=param_vars['operator'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Max Workers
    ttk.Label(main_frame, text="Max Workers:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['max_workers'] = tk.IntVar(value=6)
    ttk.Entry(main_frame, textvariable=param_vars['max_workers'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Select Soma
    ttk.Label(main_frame, text="Select Soma?").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['select_soma'] = tk.IntVar(value=False)
    ttk.Entry(main_frame, textvariable=param_vars['select_soma'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Add some spacing
    row += 1
    
    # Result variable
    result = {'params': None, 'cancelled': True}
    
    def on_ok():
        try:
            result['params'] = {
                # 'discardInitial_s': param_vars['discardInitial_s'].get(),
                'analyzeHz': param_vars['analyzeHz'].get(),
                'activityChannel': param_vars['activityChannel'].get(),
                'decayTau_s': param_vars['decayTau_s'].get(),
                'baselineWindow_s': param_vars['baselineWindow_s'].get(),
                'dXY': param_vars['dXY'].get(),
                'sparse_fac': float(np.exp(param_vars['sparse_fac_log'].get())),
                'denoiseWindow_s': param_vars['denoiseWindow_s'].get(),
                'operator': param_vars['operator'].get(),
                'max_workers': param_vars['max_workers'].get(),
                'select_soma': param_vars['select_soma'].get()
            }
            result['cancelled'] = False
            root.destroy()
        except Exception as e:
            messagebox.showerror("Error", f"Invalid parameter values: {e}")
    
    def on_cancel():
        result['cancelled'] = True
        root.destroy()
    
    # Buttons
    button_frame = ttk.Frame(main_frame)
    button_frame.grid(row=row, column=0, columnspan=2, pady=10)
    
    ttk.Button(button_frame, text="OK", command=on_ok).pack(side=tk.LEFT, padx=5)
    ttk.Button(button_frame, text="Cancel", command=on_cancel).pack(side=tk.LEFT, padx=5)
    
    # Update window to calculate required size
    root.update_idletasks()
    
    # Get the required size
    width = main_frame.winfo_reqwidth()  # Add padding
    height = main_frame.winfo_reqheight()  # Add padding
    
    # Set window size to fit contents
    root.geometry(f"{width}x{height}")
    
    # Center the window
    root.transient()
    root.grab_set()
    root.mainloop()
    
    return result

def main():
    # load SLAP2 data folder
    dr = filedialog.askdirectory(initialdir = 'Z:\\scratch\\ophys\\Michael', \
                                        title = "Select data directory")
    # dr = 'Z:\\scratch\\ophys\\Michael\\slap2_integration+raster\\slap2_760268_2024-11-05_12-35-49\\fov1\\experiment1'
    print(dr)

    # Show GUI and get parameters
    gui_result = create_parameter_gui()
    
    if gui_result['cancelled']:
        print("Parameter selection cancelled by user")
        return

    start_time = datetime.now(timezone.utc).astimezone()

    # Extract metadata parameters
    params = gui_result['params']
    print("selected params:")
    print(params)

    # find trialTable.mat file in data directory
    trialTableFile = os.path.join(dr, 'trialTable.mat')
    if not os.path.exists(trialTableFile):
        raise FileNotFoundError(f"Trial table file not found at: {trialTableFile}")
    print(trialTableFile)

    ''' trialTable structure
    trialTable = 

    struct with fields:

                    refStack: {[1×1 struct]  [1×1 struct]}
                    filename: {2×16 cell}
                    firstLine: [2×16 double]
                    lastLine: [2×16 double]
            trialEndTimeFromPC: [7.3956e+05 7.3956e+05 … ] (1×16 double)
        trialStartTimeInferred: [7.3956e+05 7.3956e+05 … ] (1×16 double)
                trueTrialIx: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]
                        epoch: [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
                    lookupFile: 'Z:\scratch\ophys\Michael\slap2_integration+raster\slap2_760268_2024-11-05_12-35-49\fov1\experiment1\\integrationRegLookupTable.mat'
                    fnRegDSInt: {2×16 cell}
                    fnAdataInt: {2×16 cell}
                alignParamsInt: [1×1 struct]
                    fnRegDS: {2×16 cell}
                    fnAdata: {2×16 cell}
                alignParams: [1×1 struct]
    '''

    # Get the struct
    trialTable = spio.loadmat(trialTableFile)['trialTable'][0,0]

    align_params = {}
    tmp = trialTable['alignParamsInt'].reshape(-1)[0]
    for key in trialTable['alignParamsInt'].dtype.names:
        align_params[key] = to_serializable(np.squeeze(tmp[key]))

    nDMDs = trialTable['filename'].shape[0]
    nTrials = trialTable['filename'].shape[1]

    # Verify all required files exist
    keepTrials = np.ones((nDMDs, nTrials), dtype=bool)

    for trialIx in range(nTrials-1, -1, -1):
        for DMDix in range(nDMDs):                
            # Check alignment data files
            align_fn = os.path.splitext(os.path.basename(trialTable['fnAdataInt'][DMDix,trialIx][0]))[0]
            if not os.path.exists(os.path.join(dr, align_fn + '.mat')):
                print(f'Missing alignData file: {align_fn}')
                keepTrials[DMDix,trialIx] = False
                
            # Check source data files
            source_fn = trialTable['filename'][DMDix,trialIx][0]
            if not os.path.exists(os.path.join(dr, source_fn)):
                print(f'Missing source data file: {source_fn}')
                keepTrials[DMDix,trialIx] = False

    if not np.all(keepTrials):
        print(f'Files were missing for {np.sum(~keepTrials)} recordings; likely failed alignments. Proceeding without them.')

    if not np.any(keepTrials):
        raise RuntimeError('All trials were rejected due to missing alignment files!')

    firstValidTrial = np.where(np.all(keepTrials,axis=0))[0][0]

    # Create ExperimentSummary directory if it doesn't exist
    savedr = os.path.join(dr, 'ExperimentSummary')
    if not os.path.exists(savedr):
        os.makedirs(savedr)
    output_h5_filename = os.path.join(savedr, f'experiment_summary_{datetime.now().strftime("%Y%m%d-%H%M%S")}.h5')

    # Load aData file
    # aData = spio.loadmat(trialTable['fnAdataInt'][0,firstValidTrial][0])['aData'][0,0]

    # params['numChannels'] = aData['numChannels'][0,0]
    # params['alignHz'] = aData['alignHz'][0,0]

    # Get the lookup file path
    lookupFile = trialTable['lookupFile'][0]

    ''' Lookup table structure
    lookupTable = 

    struct with fields:

        likelihood_means: {[5-D single]  [5-D single]}
        allSuperPixelIDs: {[3786×1 uint32]  [2699×1 uint32]}
        sparseMaskInds: {[34102×2 uint32]  [24255×2 uint32]}
                    xPre: 25
                xPost: 25
                    yPre: 25
                yPost: 25
                    zPre: 10
                zPost: 10
    '''

    with h5py.File(lookupFile, 'r') as f:
        lt = f['lookupTable']
        refs = f['#refs#']
        
        # see if this can be made compatible for variable nDMDs
        # allSuperPixelIDs
        allSuperPixelIDs_refs = lt['allSuperPixelIDs'][:].flat
        allSuperPixelIDs = {'DMD1': refs[allSuperPixelIDs_refs[0]][:].T.astype(np.int32),
                            'DMD2': refs[allSuperPixelIDs_refs[1]][:].T.astype(np.int32)}
        
        # sparseMaskInds
        sparseMaskInds_refs = lt['sparseMaskInds'][:].flat
        sparseMaskInds = {'DMD1': refs[sparseMaskInds_refs[0]][:].T.astype(np.int32),
                        'DMD2': refs[sparseMaskInds_refs[1]][:].T.astype(np.int32)}
        
        fastZ2RefZ_refs = lt['fastZ2RefZ'][:].flat
        fastZ2RefZ = {'DMD1': refs[fastZ2RefZ_refs[0]][:].T.astype(np.int32),
                    'DMD2': refs[fastZ2RefZ_refs[1]][:].T.astype(np.int32)}

    del sparseMaskInds_refs, allSuperPixelIDs_refs, fastZ2RefZ_refs

    # Extract superpixel locations and print shapes to verify
    print("Shapes:")
    subsampleMatrixInds = {}
    for DMDix in range(nDMDs):
        numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
        subsampleMatrixInds[f'DMD{DMDix+1}'] = np.zeros((numSuperPixels,2), dtype=np.int32)
        for spIdx in range(numSuperPixels):
            currSpInds = np.where(sparseMaskInds[f'DMD{DMDix+1}'][:,1] == spIdx+1)[0]
            currSpOpenPixs = sparseMaskInds[f'DMD{DMDix+1}'][currSpInds,0]-1
            spRefPix = currSpOpenPixs[np.floor(len(currSpOpenPixs)/2).astype(int)]
            subsampleMatrixInds[f'DMD{DMDix+1}'][spIdx,0] = spRefPix
            subsampleMatrixInds[f'DMD{DMDix+1}'][spIdx,1] = spIdx+1

        print(f"allSuperPixelIDs DMD{DMDix+1}: {allSuperPixelIDs['DMD'+str(DMDix+1)].shape}")
        print(f"sparseMaskInds DMD{DMDix+1}: {sparseMaskInds['DMD'+str(DMDix+1)].shape}")

    # get integration reference stack
    refStack = {}
    for DMDix in range(nDMDs):
        pattern = f"**/*DMD{DMDix+1}_CONFIG2-REFERENCE*"
        matching_files = list(Path(dr).glob(pattern))
        if matching_files:
            first_file = str(matching_files[0])
            print(f"DMD{DMDix+1} ref stack file: {first_file}")
            numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
            refStackTmp = skimio.imread(first_file) / 100
            refStackTmp = refStackTmp.reshape(-1, numChannels, refStackTmp.shape[1], refStackTmp.shape[2]).transpose(1,0,2,3)
            refStack[f'DMD{DMDix+1}'] = refStackTmp
        else:
            pattern = f"**/*DMD{DMDix+1}-REFERENCE*"
            matching_files = list(Path(dr).glob(pattern))
            if matching_files:
                first_file = str(matching_files[0])
                print(f"DMD{DMDix+1} ref stack file: {first_file}")
                numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
                refStackTmp = skimio.imread(first_file) / 100
                refStackTmp = refStackTmp.reshape(-1, numChannels, refStackTmp.shape[1], refStackTmp.shape[2]).transpose(1,0,2,3)
                refStack[f'DMD{DMDix+1}'] = refStackTmp
            else:
                print(f"No matching files found for DMD{DMDix+1}")

    # Optional manual soma selection on reference images (scroll through fast-Z planes)
    soma_masks = {}
    soma_sps = {}
    params['alignHz'] = {}
    if params['select_soma']:
        print('Manually select soma mask on each DMD based on the reference image')
        print('Controls: E=edit/add ROIs on current plane, N/P=next/prev plane, ESC/Q=finish DMD')
        for DMDix in range(nDMDs):
            aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,firstValidTrial][0])['aData'][0,0]
            params['numChannels'] = aData['numChannels'][0,0]
            params['alignHz'][f'DMD{DMDix+1}'] = aData['alignHz'][0,0]

            avg_motionR = int(np.nanmedian(np.round(aData['motionDSr'].T[0])))
            avg_motionC = int(np.nanmedian(np.round(aData['motionDSc'].T[0])))
            avg_motionZ = int(np.nanmedian(np.round(aData['motionDSz'].T[0])))

            soma_sps[f'DMD{DMDix+1}'] = []
            ref = refStack[f'DMD{DMDix+1}']  # [channels, z, y, x]
            num_ref_z = ref.shape[1]
            yx_shape = (ref.shape[2], ref.shape[3])
            # choose channel with largest mean intensity
            ch_means = [np.nanmean(ref[c]) for c in range(ref.shape[0])]
            best_ch = int(np.argmax(ch_means)) if len(ch_means) > 0 else 0

            # Map fast-Z to ref-Z for viewing; adjust potential 1-based indexing
            z_map = np.array(fastZ2RefZ[f'DMD{DMDix+1}'] + avg_motionZ).reshape(-1) - 1
            num_fast_z = z_map.shape[0]
            
            sp_fastz = subsampleMatrixInds[f'DMD{DMDix+1}'][:,0] // (yx_shape[0]*yx_shape[1])
            sp_cols = avg_motionC + (subsampleMatrixInds[f'DMD{DMDix+1}'][:,0] - sp_fastz * (yx_shape[0]*yx_shape[1])) // yx_shape[0]
            sp_rows = avg_motionR + subsampleMatrixInds[f'DMD{DMDix+1}'][:,0] % yx_shape[0]
            sp_mask = np.zeros((num_fast_z, *yx_shape), dtype=bool)
            
            # clip to valid bounds
            valid = (sp_rows >= 0) & (sp_rows < yx_shape[0]) & (sp_cols >= 0) & (sp_cols < yx_shape[1]) & (sp_fastz >= 0) & (sp_fastz < num_fast_z)
            if np.any(valid):
                sp_mask[sp_fastz[valid], sp_rows[valid], sp_cols[valid]] = True

            # mask in fast-Z space (store labels per plane)
            mask_fastz = np.zeros((num_fast_z, *yx_shape), dtype=np.int8)

            window_name = f'Select soma ROI(s) DMD{DMDix+1}'
            cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
            cv2.resizeWindow(window_name, 800, 500)
            cv2.createTrackbar('z', window_name, 0, max(0, num_fast_z-1), lambda v: None)

            while True:
                curr_fz = int(np.clip(cv2.getTrackbarPos('z', window_name), 0, max(0, num_fast_z-1)))
                refz = int(np.clip(z_map[curr_fz], 0, max(0, num_ref_z-1)))

                plane = ref[best_ch, refz]
                im = np.nan_to_num(plane, nan=0.0)
                vmin = np.percentile(im, 1)
                vmax = np.percentile(im, 99.5)
                if not np.isfinite(vmin):
                    vmin = float(np.nanmin(im)) if np.any(np.isfinite(im)) else 0.0
                if not np.isfinite(vmax):
                    vmax = float(np.nanmax(im)) if np.any(np.isfinite(im)) else 1.0
                if vmax <= vmin:
                    vmax = vmin + 1.0
                im8 = np.clip((im - vmin) / (vmax - vmin), 0, 1)
                im8 = (im8 * 255).astype(np.uint8)

                # show overlay text
                disp = cv2.cvtColor(im8, cv2.COLOR_GRAY2BGR)
                text = f'FastZ {curr_fz+1}/{num_fast_z} (RefZ {refz+1}/{num_ref_z}) | E=edit, N/P=nav, ESC=done'
                cv2.putText(disp, text, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0,255,0), 2, cv2.LINE_AA)
                # overlay superpixel mask (green) on current plane
                if np.any(sp_mask[curr_fz]):
                    sp_color = np.zeros_like(disp)
                    sp_color[sp_mask[curr_fz]] = (0, 255, 0)
                    alpha_sp = 0.25
                    disp = cv2.addWeighted(disp, 1 - alpha_sp, sp_color, alpha_sp, 0)
                # overlay drawn soma mask (red)
                if (mask_fastz[curr_fz] > 0).any():
                    roi_color = np.zeros_like(disp)
                    roi_color[mask_fastz[curr_fz] > 0] = (0, 0, 255)
                    alpha_roi = 0.3
                    disp = cv2.addWeighted(disp, 1 - alpha_roi, roi_color, alpha_roi, 0)
                cv2.imshow(window_name, disp)

                key = cv2.waitKey(50) & 0xFF
                if key == ord('e'):
                    edit_name = f'Edit ROIs z={curr_fz+1}'
                    # prepare edit image with superpixel overlay for context
                    edit_disp = cv2.cvtColor(im8, cv2.COLOR_GRAY2BGR)
                    if np.any(sp_mask[curr_fz]):
                        sp_color = np.zeros_like(edit_disp)
                        sp_color[sp_mask[curr_fz]] = (0, 255, 0)
                        alpha_sp = 0.25
                        edit_disp = cv2.addWeighted(edit_disp, 1 - alpha_sp, sp_color, alpha_sp, 0)
                        rois = cv2.selectROIs(edit_name, edit_disp, showCrosshair=True, fromCenter=False)
                        cv2.resizeWindow(edit_name, 800, 500)
                        if rois is not None and len(rois) > 0:
                            next_label = int(np.max(mask_fastz)) + 1
                            for (x, y, w, h) in rois:
                                mask_fastz[curr_fz, y:y+h, x:x+w] = next_label
                                next_label += 1
                            print(f'Added {len(rois)} ROI(s) at fast-Z {curr_fz+1} for DMD{DMDix+1}')
                        cv2.destroyWindow(edit_name)
                elif key == ord('n'):
                    curr_fz = min(curr_fz + 1, num_fast_z - 1)
                    cv2.setTrackbarPos('z', window_name, curr_fz)
                elif key == ord('p'):
                    curr_fz = max(curr_fz - 1, 0)
                    cv2.setTrackbarPos('z', window_name, curr_fz)
                elif key in (27, ord('q')):
                    break

            cv2.destroyWindow(window_name)

            soma_masks[f'DMD{DMDix+1}'] = mask_fastz

            for roi in range(np.max(mask_fastz)):
                soma_sps[f'DMD{DMDix+1}'].append(np.flatnonzero(mask_fastz[sp_fastz, sp_rows, sp_cols] == roi + 1))

    # Get dilation 17 PSF or load in from file
    pattern = f"**/*DMD*-REFERENCE*"
    matching_files = list(Path(dr).glob(pattern))
    if matching_files:
        first_file = str(matching_files[0])
        first_file_name = os.path.splitext(os.path.basename(first_file))[0].split('_DMD')[0]
        psf_path = os.path.join(dr, first_file_name + '-PSF.tif')
    else:
        psf_path = os.path.join(dr, 'refStack-PSF.tif')
        print(f"No ref stack found")

    if not os.path.exists(psf_path):
        psf = {}
        for DMDix in range(nDMDs):
            print(f"Calculating PSF for DMD{DMDix+1}")
            psf[f'DMD{DMDix+1}'] = skimio.imread(os.path.join(str(Path(__file__).parent), 'psfs', 'dil-17.tif')) # todo: make dilation a parameter
        combined_dims = [0,0]
        for DMDix in range(nDMDs):
            combined_dims[0] = max(psf[f'DMD{DMDix+1}'].shape[0], combined_dims[0])
            combined_dims[1] = max(psf[f'DMD{DMDix+1}'].shape[1], combined_dims[1])

        psf_combined = np.zeros((nDMDs, combined_dims[0], combined_dims[1]), dtype=np.float32)
        for DMDix in range(nDMDs):
            psf_combined[DMDix] = np.pad(psf[f'DMD{DMDix+1}'], (((combined_dims[0] - psf[f'DMD{DMDix+1}'].shape[0])//2,(combined_dims[0] - psf[f'DMD{DMDix+1}'].shape[0])//2), ((combined_dims[1] - psf[f'DMD{DMDix+1}'].shape[1])//2,(combined_dims[1] - psf[f'DMD{DMDix+1}'].shape[1])//2)), constant_values = np.min(psf[f'DMD{DMDix+1}']))

        skimio.imsave(psf_path, psf_combined.astype(np.float32))
    else:
        psf_combined = skimio.imread(psf_path)
    
    psf_combined[0][psf_combined[0] < np.max(psf_combined[0])*np.exp(-3)] = 0
    psf_combined[1][psf_combined[1] < np.max(psf_combined[1])*np.exp(-3)] = 0

    psf = {
        'DMD1': psf_combined[0],
        'DMD2': psf_combined[1]
    }

    del psf_combined

    # Crop PSFs to remove boundary zeros
    for DMDix in range(nDMDs):
        # Find non-zero rows and columns
        non_zero_rows = np.any(psf[f'DMD{DMDix+1}'] != 0, axis=1)
        non_zero_cols = np.any(psf[f'DMD{DMDix+1}'] != 0, axis=0)
        
        # Get indices of first and last non-zero rows/cols
        row_start, row_end = np.where(non_zero_rows)[0][[0, -1]]
        col_start, col_end = np.where(non_zero_cols)[0][[0, -1]]
        
        # Crop the PSF
        psf[f'DMD{DMDix+1}'] = psf[f'DMD{DMDix+1}'][row_start:row_end+1, col_start:col_end+1]

    for DMDix in range(nDMDs): 
        # range(nDMDs-1, -1, -1):
        print(f'Processing DMD{DMDix+1}')

        dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
        dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
        numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
        numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]
        numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]

        nPixels = dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs

        refD, refC, refR = ref_pixs_to_drc(subsampleMatrixInds[f'DMD{DMDix+1}'][:, 0], dmdPixelsPerColumn, dmdPixelsPerRow)

        psf_cur = psf[f'DMD{DMDix+1}'].astype(np.float32, copy=False)
        psf_h, psf_w = psf_cur.shape
        filterSize = psf_h * psf_w
        plane_size = dmdPixelsPerColumn * dmdPixelsPerRow

        row_offsets = (np.arange(psf_h, dtype=np.int32) - (psf_h // 2)).reshape(-1, 1)
        col_offsets = (np.arange(psf_w, dtype=np.int32) - (psf_w // 2)).reshape(1, -1)
        row_offsets_flat = np.broadcast_to(row_offsets, (psf_h, psf_w)).ravel()
        col_offsets_flat = np.broadcast_to(col_offsets, (psf_h, psf_w)).ravel()
        psf_vals_flat = psf_cur.ravel()

        sparseHInds = np.zeros((2, numSuperPixels * filterSize), dtype=np.int32)
        sparseHVals = np.zeros((numSuperPixels * filterSize,), dtype=np.float32)
        sparseHInds[0] = np.repeat(subsampleMatrixInds[f'DMD{DMDix+1}'][:, 1] - 1, filterSize)
        for spIdx in range(subsampleMatrixInds[f'DMD{DMDix+1}'].shape[0]):
            start = spIdx * filterSize
            end = (spIdx + 1) * filterSize
            rows = refR[spIdx] + row_offsets_flat
            cols = refC[spIdx] + col_offsets_flat
            sparseHInds[1, start:end] = (
                refD[spIdx] * plane_size + rows * dmdPixelsPerRow + cols
            ).astype(np.int32, copy=False)
            sparseHVals[start:end] = psf_vals_flat
        non_zero_mask = sparseHVals != 0
        sparseHVals = sparseHVals[non_zero_mask]
        sparseHInds = sparseHInds[:, non_zero_mask]

        trial_info = [(i, keepTrials[DMDix,i]) for i in range(nTrials)]

        get_trial_data_partial = partial(
            get_trial_data,
            DMDix=DMDix,
            params=params,
            sampFreq=params['alignHz'][f'DMD{DMDix+1}'],
            refStack=refStack,
            fastZ2RefZ=fastZ2RefZ,
            allSuperPixelIDs=allSuperPixelIDs,
            dr=dr,
            trialTable=trialTable,
            all_channels=True
        )

        data_file = os.path.join(dr, f'lowres_data_DMD{DMDix+1}.npz')
        
        if os.path.exists(data_file):
            print(f'Loading existing low resolution data from {data_file}')
            with np.load(data_file) as data_arrays:
                lowResData = data_arrays['lowResData']
                lowResDataCt = data_arrays['lowResDataCt'] 
                lowResMotionR = data_arrays['lowResMotionR']
                lowResMotionC = data_arrays['lowResMotionC']
                lowResMotionZ = data_arrays['lowResMotionZ']
                lowResTrialID = data_arrays['lowResTrialID']
                lowResData2 = data_arrays['lowResData2']
                lowResDataCt2 = data_arrays['lowResDataCt2']
        else:
            with mp.Pool(processes=min(params['max_workers'],min(mp.cpu_count(),len(trial_info)))) as pool:
                results = list(pool.imap(get_trial_data_partial, trial_info))

            lowResData = np.concatenate([r[0] for r in results if r is not None], axis=1)
            lowResDataCt = np.concatenate([r[1] for r in results if r is not None], axis=1)
            lowResMotionR = np.concatenate([r[2]['motionDSr'] for r in results if r is not None], axis=0)
            lowResMotionC = np.concatenate([r[2]['motionDSc'] for r in results if r is not None], axis=0)
            lowResMotionZ = np.concatenate([r[2]['motionDSz'] for r in results if r is not None], axis=0)
            lowResTrialID = np.concatenate([np.ones_like(r[3])*i for i,r in enumerate(results) if r is not None], axis=0)
            lowResData2 = np.concatenate([r[4] for r in results if r is not None], axis=1)
            lowResDataCt2 = np.concatenate([r[5] for r in results if r is not None], axis=1)

            # Save data arrays
            data_arrays = {
                'lowResData': lowResData,
                'lowResDataCt': lowResDataCt,
                'lowResMotionR': lowResMotionR,
                'lowResMotionC': lowResMotionC,
                'lowResMotionZ': lowResMotionZ,
                'lowResTrialID': lowResTrialID,
                'lowResData2': lowResData2,
                'lowResDataCt2': lowResDataCt2
            }
            np.savez(data_file, **data_arrays)
            print(f'Saved low resolution data to {data_file}')

        lowResDataNorm = lowResData / lowResDataCt
        lowResData2Norm = lowResData2 / lowResDataCt2
        vIM = 1 / lowResDataCt
        del lowResData, lowResData2, lowResDataCt, lowResDataCt2
        
        uniqueMotion, motInds = np.unique(np.round(np.concatenate((lowResMotionR,lowResMotionC,lowResMotionZ),axis=1)),axis=0,return_inverse=True)
        # Calculate the median from the full lowResMotionZ array, then filter uniqueMotion directly
        median_z = np.median(lowResMotionZ)
        # Boolean mask: True if within 1.5 um from median Z and has >100 frames
        bin_counts = np.bincount(motInds, minlength=uniqueMotion.shape[0])
        keep_mask = (np.abs(uniqueMotion[:, 2] - median_z) <= 1.5) & (bin_counts > 100)
        motIndsToKeep = np.nonzero(keep_mask)[0]

        framesToKeep = np.isin(motInds, motIndsToKeep)

        # after filtering out frames with large Z motion, we can consider there to be no Z motion in the remaining frames

        mean_im = np.full((2,numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), np.nan, dtype=np.float32)
        most_common_mot = np.argmax(np.bincount(motInds))
        mean_im[0,refD,refR + int(uniqueMotion[most_common_mot,0]),refC + int(uniqueMotion[most_common_mot,1])] = np.nanmean(lowResDataNorm[:,framesToKeep], axis=1)
        mean_im[1,refD,refR + int(uniqueMotion[most_common_mot,0]),refC + int(uniqueMotion[most_common_mot,1])] = np.nanmean(lowResData2Norm[:,framesToKeep], axis=1)
        del lowResData2Norm

        motIndsYX = -1*np.ones((len(motInds),), dtype=np.int32)
        uniqueMotionToKeepYX, motIndsYX[framesToKeep] = np.unique(np.round(np.concatenate((lowResMotionR,lowResMotionC),axis=1)[framesToKeep,:]),axis=0,return_inverse=True)

        selPixMask = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), dtype=bool)
        for i in range(len(uniqueMotionToKeepYX)):
            selPixMask[refD,refR + int(uniqueMotionToKeepYX[i,0]),refC + int(uniqueMotionToKeepYX[i,1])] = True
        selPixMask = ndimage.binary_dilation(selPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
        selPixIdxs = np.flatnonzero(selPixMask)

        # Get pixel coordinates for selected pixels
        # [z, y, x] coordinates, TODO: should update the name to selPixCoords
        pixel_coords = np.empty((len(selPixIdxs), 3), dtype=np.int32)
        pixel_coords[:,0] = selPixIdxs // (dmdPixelsPerColumn * dmdPixelsPerRow)
        remainder = selPixIdxs % (dmdPixelsPerColumn * dmdPixelsPerRow)
        pixel_coords[:,1] = remainder // dmdPixelsPerRow
        pixel_coords[:,2] = remainder % dmdPixelsPerRow
        
        # Convert to torch tensor
        pixel_coords_tensor = torch.tensor(pixel_coords, dtype=torch.float32)
        
        # Convert PSF to torch tensor and extract values where valid
        psf_tensor = torch.from_numpy(psf[f'DMD{DMDix+1}']).float()
        psf_tensor = psf_tensor / torch.sum(psf_tensor)

        # ----- Create expanded PSF with 5x resolution, maintaining center alignment -----
        ex_fac = 2
        psf_shape = psf[f'DMD{DMDix+1}'].shape
        psf_center = (psf_shape[0] // 2, psf_shape[1] // 2)
        psf_tensor_expanded = torch.zeros((psf_shape[0]*ex_fac, psf_shape[1]*ex_fac), dtype=torch.float32) # TODO: consider whether to keep expanded PSF
        # Calculate center points
        expanded_center_y = (psf_shape[0]*ex_fac - 1) / 2
        expanded_center_x = (psf_shape[1]*ex_fac - 1) / 2
        # Create coordinate grids centered around zero
        orig_y = torch.linspace(-expanded_center_y, expanded_center_y, psf_shape[0])
        orig_x = torch.linspace(-expanded_center_x, expanded_center_x, psf_shape[1])
        expanded_y = torch.linspace(-expanded_center_y, expanded_center_y, psf_shape[0]*ex_fac)
        expanded_x = torch.linspace(-expanded_center_x, expanded_center_x, psf_shape[1]*ex_fac)
        # Interpolate PSF values while maintaining center alignment
        interp_spline = RectBivariateSpline(orig_y.numpy(), orig_x.numpy(), psf_tensor.numpy())
        psf_tensor_expanded = torch.from_numpy(interp_spline(expanded_y.numpy(), expanded_x.numpy())).float()
        # Normalize expanded PSF
        psf_tensor_expanded = psf_tensor_expanded / torch.sum(psf_tensor_expanded)
        psf_shape_expanded = psf_tensor_expanded.shape
        psf_center_expanded = (psf_shape_expanded[0] // 2, psf_shape_expanded[1] // 2)

        # build D matrices for each plane (convolution matrix for PSF)
        D = [None] * numFastZs
        D_expanded = [None] * numFastZs
        interp_data = np.full((selPixIdxs.shape[0],lowResDataNorm.shape[1]), np.nan, dtype=np.float32)
        plane_size = dmdPixelsPerColumn * dmdPixelsPerRow
        sel_pix_z = selPixIdxs // plane_size
        sel_pix_remainder = selPixIdxs % plane_size
        z_sel_indices = [np.flatnonzero(sel_pix_z == z) for z in range(numFastZs)]
        refD_z_indices = [np.flatnonzero(refD == z) for z in range(numFastZs)]

        # interp data is for estimating background fluorescence
        def build_interp_data(data, refR, refC, sel_pixels_2d,
                        H=800, W=1280, dtype=np.float32):
            # ---- normalize dtypes once
            refR = np.asarray(refR, dtype=np.int32)
            refC = np.asarray(refC, dtype=np.int32)
            data = np.asarray(data, dtype=dtype)

            r0 = int(max(0, refR.min()) + uniqueMotionToKeepYX[:, 0].min())
            r1 = int(min(H, refR.max()) + uniqueMotionToKeepYX[:, 0].max())
            c0 = int(max(0, refC.min()) + uniqueMotionToKeepYX[:, 1].min())
            c1 = int(min(W, refC.max()) + uniqueMotionToKeepYX[:, 1].max())

            out = np.full((sel_pixels_2d.shape[0], data.shape[1]), np.nan, dtype=dtype)

            # Cache selected pixels by column once.
            sel_cols = sel_pixels_2d[:, 1]
            sel_rows = sel_pixels_2d[:, 0]
            unique_sel_cols = np.unique(sel_cols)
            selPixIdxs_by_col_all = {int(c): np.flatnonzero(sel_cols == c) for c in unique_sel_cols}

            # ---- frames per motion bucket (avoid (motInds==idx) every time)
            frames_by_motion = [np.flatnonzero(motIndsYX == idx) for idx in range(len(uniqueMotionToKeepYX))]

            for m_idx, frames in enumerate(frames_by_motion):
                if frames.size == 0:
                    continue

                # shift once per motion group
                sR = refR + int(uniqueMotionToKeepYX[m_idx, 0])
                sC = refC + int(uniqueMotionToKeepYX[m_idx, 1])

                # columns that actually land in our window and have any samples
                in_win = (sR >= r0) & (sR <= r1) & (sC >= c0) & (sC <= c1)
                if not np.any(in_win):
                    continue
                cols = np.unique(sC[in_win])

                rows_by_col = {}
                idxs_by_col = {}
                single_points = []

                for c in cols:
                    c_int = int(c)
                    if c_int not in selPixIdxs_by_col_all:
                        continue
                    mask = (sC == c)
                    data_ix = np.flatnonzero(mask)
                    if data_ix.size == 1:
                        rr = int(sR[data_ix[0]])
                        col_sel = selPixIdxs_by_col_all[c_int]
                        rr_match = col_sel[sel_rows[col_sel] == rr]
                        if rr_match.size > 0 and (r0 <= rr <= r1):
                            single_points.append((rr_match, int(data_ix[0])))
                        continue

                    rows = sR[data_ix]
                    order = rows.argsort(kind="mergesort")
                    rows_by_col[c_int] = rows[order]
                    idxs_by_col[c_int] = data_ix[order]

                # Interpolate per column, but batch all frames in this motion bucket.
                for c_int, rows in rows_by_col.items():
                    colSelPixIdxs = selPixIdxs_by_col_all[c_int]
                    target_rows = sel_rows[colSelPixIdxs]
                    src_idx = idxs_by_col[c_int]

                    # Vectorized linear interpolation with NaN outside range.
                    pos = np.searchsorted(rows, target_rows, side='left')
                    in_bounds = pos < rows.size
                    exact = np.zeros_like(pos, dtype=bool)
                    exact[in_bounds] = rows[pos[in_bounds]] == target_rows[in_bounds]

                    if np.any(exact):
                        exact_rows = np.flatnonzero(exact)
                        exact_src = src_idx[pos[exact]]
                        out[np.ix_(colSelPixIdxs[exact_rows], frames)] = data[exact_src][:, frames]

                    lo = pos - 1
                    hi = pos
                    interp_mask = (~exact) & (lo >= 0) & (hi < rows.size)
                    if not np.any(interp_mask):
                        continue

                    interp_rows = np.flatnonzero(interp_mask)
                    lo_v = lo[interp_mask]
                    hi_v = hi[interp_mask]
                    row_lo = rows[lo_v].astype(dtype, copy=False)
                    row_hi = rows[hi_v].astype(dtype, copy=False)
                    denom = row_hi - row_lo
                    nonzero = denom != 0
                    if not np.any(nonzero):
                        continue

                    interp_rows = interp_rows[nonzero]
                    lo_v = lo_v[nonzero]
                    hi_v = hi_v[nonzero]
                    row_lo = row_lo[nonzero]
                    row_hi = row_hi[nonzero]
                    alpha = ((target_rows[interp_rows].astype(dtype) - row_lo) / (row_hi - row_lo))[:, None]

                    vals_lo = data[src_idx[lo_v]][:, frames]
                    vals_hi = data[src_idx[hi_v]][:, frames]
                    out[np.ix_(colSelPixIdxs[interp_rows], frames)] = (1.0 - alpha) * vals_lo + alpha * vals_hi

                # Keep single points (no interpolation), like original behavior.
                for rr_match, src_pos in single_points:
                    out[np.ix_(rr_match, frames)] = data[src_pos, frames]

            return out, (r0, r1), (c0, c1)

        for z in tqdm(range(numFastZs),desc='Computing D matrices for each plane'):
            z_idxs = z_sel_indices[z]
            if z_idxs.size == 0:
                D[z] = torch.zeros((0, 0), dtype=torch.float32)
                D_expanded[z] = torch.zeros((0, 0), dtype=torch.float32)
                continue

            # Pre-allocate result matrix
            D[z] = torch.zeros((z_idxs.size, z_idxs.size), dtype=torch.float32)

            # Convert selected pixel indices to 2D coordinates (vectorized)
            sel_pixels_2d = np.column_stack([
                    sel_pix_remainder[z_idxs] // dmdPixelsPerRow,  # rows
                    sel_pix_remainder[z_idxs] % dmdPixelsPerRow    # cols
                ])

            # Vectorized computation using broadcasting
            src_rows = sel_pixels_2d[:, 0][np.newaxis, :]  # Shape: (1, n_sources)
            src_cols = sel_pixels_2d[:, 1][np.newaxis, :]  # Shape: (1, n_sources)
            tgt_rows = sel_pixels_2d[:, 0][:, np.newaxis]  # Shape: (n_targets, 1)
            tgt_cols = sel_pixels_2d[:, 1][:, np.newaxis]  # Shape: (n_targets, 1)

            # Calculate relative positions
            rel_rows = tgt_rows - src_rows + psf_center[0]  # Shape: (n_targets, n_sources)
            rel_cols = tgt_cols - src_cols + psf_center[1]  # Shape: (n_targets, n_sources)

            # Create mask for valid PSF indices
            valid_mask = ((rel_rows >= 0) & (rel_rows < psf_shape[0]) & 
                        (rel_cols >= 0) & (rel_cols < psf_shape[1]))

            D[z][valid_mask] = psf_tensor[rel_rows[valid_mask], rel_cols[valid_mask]]

            D_expanded[z] = torch.zeros_like(D[z], dtype=torch.float32)

            rel_rows_expanded = tgt_rows - src_rows + psf_center_expanded[0]
            rel_cols_expanded = tgt_cols - src_cols + psf_center_expanded[1]

            valid_mask_expanded = ((rel_rows_expanded >= 0) & (rel_rows_expanded < psf_shape_expanded[0]) & 
                                (rel_cols_expanded >= 0) & (rel_cols_expanded < psf_shape_expanded[1]))

            D_expanded[z][valid_mask_expanded] = psf_tensor_expanded[rel_rows_expanded[valid_mask_expanded], rel_cols_expanded[valid_mask_expanded]]

            ref_z_idxs = refD_z_indices[z]
            interp_data[z_idxs], _, _ = build_interp_data(lowResDataNorm[ref_z_idxs], refR[ref_z_idxs], refC[ref_z_idxs], sel_pixels_2d)

        baseline_window = int(params['alignHz'][f'DMD{DMDix+1}'] * params['baselineWindow_s'])
        valid = ~np.isnan(interp_data)
        np.nan_to_num(interp_data, copy=False, nan=0.0)

        # Use running means and scale to window sums for a faster NaN-aware baseline.
        sum_vals = ndimage.uniform_filter1d(interp_data, size=baseline_window, axis=1, mode='nearest') * baseline_window
        valid_f = valid.astype(np.float32, copy=False)
        count_vals = ndimage.uniform_filter1d(valid_f, size=baseline_window, axis=1, mode='nearest') * baseline_window

        interp_data_background = np.empty_like(sum_vals, dtype=interp_data.dtype)
        interp_data_background.fill(np.nan)
        _ = np.divide(sum_vals, count_vals, out=interp_data_background, where=count_vals > 0)

        del sum_vals, count_vals, valid_f

        # interp_data[nan_mask] = np.nanmedian(interp_data[~nan_mask])
        # interp_data_background = ndimage.uniform_filter1d(interp_data, size=baseline_window, axis=1, mode='nearest')
        # interp_data_background[nan_mask] = np.nan

        del interp_data
            
        background = np.full_like(lowResDataNorm, np.nan, dtype=np.float32)
        
        for motion_idx in tqdm(range(len(uniqueMotionToKeepYX)),desc='Computing background for all motion indices'):
            motion_frames = np.flatnonzero(motIndsYX == motion_idx)
            dR, dC = int(uniqueMotionToKeepYX[motion_idx, 0]), int(uniqueMotionToKeepYX[motion_idx, 1])
            sD = refD
            sR = refR + dR
            sC = refC + dC

            shifted_indices = sD * (dmdPixelsPerColumn * dmdPixelsPerRow) + sR * dmdPixelsPerRow + sC

            # Build shifted->selected pixel mapping once per motion index.
            bg_idxs = np.searchsorted(selPixIdxs, shifted_indices)
            idxs_mask = bg_idxs < selPixIdxs.size
            idxs_mask[idxs_mask] &= (selPixIdxs[bg_idxs[idxs_mask]] == shifted_indices[idxs_mask])

            if np.any(idxs_mask):
                # Fill all frames for this motion index in one shot.
                valid_rows = np.flatnonzero(idxs_mask)
                background[np.ix_(valid_rows, motion_frames)] = interp_data_background[np.ix_(bg_idxs[valid_rows], motion_frames)]

        # nan_mask = np.isnan(background)
        # # Impute NaN values using low-rank matrix completion (soft-impute)
        # # Initialize NaN values with global mean
        # global_mean = np.nanmean(background)
        # background_imputed = background.copy()
        # background_imputed[nan_mask] = global_mean
        
        # # Iterative soft-impute algorithm for matrix completion
        # max_iter = 20
        # rank = min(50, min(background.shape) - 1)
        # tol = 1e-4
        
        # for iter_idx in tqdm(range(max_iter), desc='Imputing background'):
        #     background_old = background_imputed.copy()
            
        #     # Truncated SVD on current estimate (much faster for large matrices)
        #     U, s, Vt = svds(background_imputed, k=rank)
            
        #     # Sort by descending singular values (svds returns ascending order)
        #     idx = np.argsort(s)[::-1]
        #     U = U[:, idx]
        #     s = s[idx]
        #     Vt = Vt[idx, :]
            
        #     # Reconstruct with truncated rank
        #     background_approx = U @ np.diag(s) @ Vt
            
        #     # Update only the missing values (keep observed values fixed)
        #     background_imputed[nan_mask] = background_approx[nan_mask]
            
        #     # Check convergence
        #     change = np.linalg.norm(background_imputed[nan_mask] - background_old[nan_mask])
        #     if change < tol * np.linalg.norm(background_imputed[nan_mask]):
        #         print(f"Converged at iteration {iter_idx+1} with change: {change:.8f}")
        #         break
        # background = background_imputed

        # del background_imputed, background_old, U, s, Vt, background_approx

        residualRaw = lowResDataNorm - background
        residual = residualRaw / ((background+0.5) * vIM)
        del interp_data_background

        # precompute H matrices
        H_mots = [None] * len(uniqueMotionToKeepYX)
        base_sparse_rows = sparseHInds[0]
        base_sparse_cols = sparseHInds[1]
        motion_shifts = (
            uniqueMotionToKeepYX[:, 0].astype(np.int64, copy=False) * dmdPixelsPerRow
            + uniqueMotionToKeepYX[:, 1].astype(np.int64, copy=False)
        )
        sparseHIndsShiftedSelPix = np.empty((2, base_sparse_cols.shape[0]), dtype=np.int64)
        sparseHIndsShiftedSelPix[0] = base_sparse_rows
        for i, pix_shift in enumerate(motion_shifts):
            shifted_cols = base_sparse_cols + pix_shift
            sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs, shifted_cols)
            H_mots[i] = torch.sparse_coo_tensor(
                sparseHIndsShiftedSelPix,
                sparseHVals,
                (numSuperPixels, selPixIdxs.shape[0]),
                dtype=torch.float32,
            )
        
        # Pre-compute z-specific column masks/remaps once (independent of motion index).
        sel_pix_z = selPixIdxs // plane_size
        z_masks_torch = [torch.from_numpy(sel_pix_z == z) for z in range(numFastZs)]
        z_col_idxs_torch = [torch.nonzero(z_mask, as_tuple=False).squeeze(1) for z_mask in z_masks_torch]
        z_remaps_torch = []
        for z in range(numFastZs):
            remap = torch.full((selPixIdxs.shape[0],), -1, dtype=torch.long)
            if z_col_idxs_torch[z].numel() > 0:
                remap[z_col_idxs_torch[z]] = torch.arange(z_col_idxs_torch[z].numel(), dtype=torch.long)
            z_remaps_torch.append(remap)

        # Compute rho on-the-fly per motion to avoid storing all HD_mots in memory.
        rho = np.full((len(selPixIdxs),lowResDataNorm.shape[1]), np.nan, dtype=np.float32)
        z_col_idxs_np = [z_idx.numpy() for z_idx in z_col_idxs_torch]
        for i in tqdm(range(len(uniqueMotionToKeepYX)), desc='Computing rho for all motion indices'):
            motion_frames = np.flatnonzero(motIndsYX == i)
            if motion_frames.size == 0:
                continue

            H_mot = H_mots[i].coalesce()
            H_idxs = H_mot.indices()
            H_vals = H_mot.values()
            nrows, _ = H_mot.size()
            residual_motion = residual[:, motion_frames].T

            validPixMask = np.zeros((numFastZs, dmdPixelsPerColumn, dmdPixelsPerRow), dtype=bool)
            validPixMask[refD, refR + int(uniqueMotionToKeepYX[i, 0]), refC + int(uniqueMotionToKeepYX[i, 1])] = True
            validPixMask = ndimage.binary_dilation(validPixMask, structure=np.ones((1, psf[f'DMD{DMDix+1}'].shape[0], psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
            validPixMask = ndimage.binary_erosion(validPixMask, structure=np.ones((1, psf[f'DMD{DMDix+1}'].shape[0], psf[f'DMD{DMDix+1}'].shape[1] * 2 - 1), dtype=bool))
            validPixIdxs = np.flatnonzero(validPixMask)
            valid_lookup = np.zeros(numFastZs * plane_size, dtype=bool)
            valid_lookup[validPixIdxs] = True
            valid_sel_cols = valid_lookup[selPixIdxs]

            for z in range(numFastZs):
                new_ncols = int(z_col_idxs_torch[z].numel())
                if new_ncols == 0:
                    continue

                zMask = z_masks_torch[z]
                keep_mask = zMask[H_idxs[1]]
                if keep_mask.sum().item() == 0:
                    continue

                remap = z_remaps_torch[z]
                new_rows = H_idxs[0, keep_mask]
                new_cols = remap[H_idxs[1, keep_mask]]
                new_idxs = torch.stack([new_rows, new_cols], dim=0)
                new_vals = H_vals[keep_mask]
                H_sub = torch.sparse_coo_tensor(new_idxs, new_vals, (nrows, new_ncols), dtype=torch.float32).coalesce()

                HD = torch.sparse.mm(H_sub, D[z])
                HD = HD / torch.sum(HD, dim=0, keepdim=True)
                HD_expanded = torch.sparse.mm(H_sub, D_expanded[z])
                HD_expanded = HD_expanded / torch.sum(HD_expanded, dim=0, keepdim=True)
                HD_diff = HD - HD_expanded

                z_cols = z_col_idxs_np[z]
                z_valid_mask = valid_sel_cols[z_cols]
                if not np.any(z_valid_mask):
                    continue
                z_valid_cols = z_cols[z_valid_mask]
                local_valid_cols = remap[torch.from_numpy(z_valid_cols)].long()

                rho[np.ix_(z_valid_cols,motion_frames)] = (residual_motion @ HD_diff[:, local_valid_cols].numpy()).T

        
        nanCt = np.mean(np.isnan(rho),axis=1)
        rho[nanCt > 0.5] = np.nan

        decayTau_frames = params['decayTau_s'] * params['alignHz'][f'DMD{DMDix+1}']
        k1d = np.exp(
            np.linspace(-np.ceil(decayTau_frames * 3), 0, int(np.ceil(decayTau_frames * 3) + 1))
            / decayTau_frames
        )
        k1d = k1d / np.sum(k1d)
        k2d = np.expand_dims(k1d, 0)
        # NaN-aware smoothing along time only (1xL kernel). Process row chunks so
        # np.nan_to_num / isfinite never build full-array boolean temps (OOM risk).
        n_time = rho.shape[1]
        bytes_per_row = max(1, n_time) * np.dtype(np.float32).itemsize * 3
        row_chunk = max(64, min(4096, int(200_000_000 // bytes_per_row)))
        # --- Temporarily disabled: robust rolling-MAD normalization (divide by MAD, scale to Gaussian std) ---
        # mad_eps = np.float32(1e-12)
        # mad_gauss_scale = np.float32(0.6741891400433162)
        # mad_win = max(len(k1d), 3)
        # pad_mad = mad_win // 2
        # mad_bytes_per_row = (n_time + 2 * pad_mad) * np.dtype(np.float32).itemsize
        # Sub-batch rolling MAD so padded (b,T) stays in cache / avoids huge temporaries.
        # mad_subbatch = max(8, min(row_chunk, int(100_000_000 // max(mad_bytes_per_row * 3, 1))))
        for r0 in tqdm(range(0, rho.shape[0], row_chunk), desc='Smoothing rho'):
            r1 = min(r0 + row_chunk, rho.shape[0])
            rc = rho[r0:r1].copy()
            row_has_data = np.any(np.isfinite(rc), axis=1)
            if not np.any(row_has_data):
                rho[r0:r1] = rc
                continue

            rc_valid = rc[row_has_data]
            rho_num = np.nan_to_num(rc_valid, nan=0.0)
            rho_den = signal.convolve(
                np.isfinite(rc_valid).astype(np.float32), k2d, mode='same'
            )
            rho_num = signal.convolve(rho_num, k2d, mode='same')
            valid_den = rho_den > 0.75
            np.divide(rho_num, rho_den, out=rho_num, where=valid_den)
            rho_num[~valid_den] = np.nan
            rc[row_has_data] = rho_num
            rho[r0:r1] = rc

        # Process frames more efficiently by avoiding unnecessary allocations
        n_frames = rho.shape[1]
        batch_size = min(1000, n_frames)  # Smaller batch size to reduce memory
        num_batches = int(np.ceil(n_frames / batch_size))
        
        # Pre-compute 2D coordinates once
        spatial_coords = np.unravel_index(selPixIdxs, (numFastZs,dmdPixelsPerColumn, dmdPixelsPerRow))
        depth_indices = spatial_coords[0][None, :].astype(np.intp, copy=False)
        row_indices = spatial_coords[1][None, :].astype(np.intp, copy=False)  # Shape: (1, len(selPixIdxs))
        col_indices = spatial_coords[2][None, :].astype(np.intp, copy=False)  # Shape: (1, len(selPixIdxs))
        
        # Pre-compute dilation structure
        dilation_struct = np.ones((3,3), dtype=np.uint8)
        
        # Initialize output array
        act_im = np.zeros((numFastZs,dmdPixelsPerColumn, dmdPixelsPerRow), dtype=np.float32)
        temporal_pad = 1

        # Crop computation to a padded ROI that contains all selected pixels.
        # This skips large all-NaN spatial regions that can never be populated.
        row_min = int(np.min(row_indices))
        row_max = int(np.max(row_indices))
        col_min = int(np.min(col_indices))
        col_max = int(np.max(col_indices))
        row_start = max(0, row_min - 1)
        row_stop = min(dmdPixelsPerColumn, row_max + 2)  # exclusive
        col_start = max(0, col_min - 1)
        col_stop = min(dmdPixelsPerRow, col_max + 2)  # exclusive
        roi_h = row_stop - row_start
        roi_w = col_stop - col_start

        row_indices_roi = (row_indices - row_start).astype(np.intp, copy=False)
        col_indices_roi = (col_indices - col_start).astype(np.intp, copy=False)

        # Preallocate reusable buffers for the largest needed batch (including temporal padding)
        prealloc_size = int(min(batch_size + 2*temporal_pad, n_frames))
        batch_rho = np.empty((prealloc_size, numFastZs, roi_h, roi_w), dtype=np.float32)
        nan_mask = np.empty_like(batch_rho, dtype=bool)
        interior_shape = (
            max(prealloc_size - 2, 0),
            numFastZs,
            max(roi_h - 2, 0),
            max(roi_w - 2, 0),
        )
        local_maxima_core = np.empty(interior_shape, dtype=bool)
        compare_tmp_core = np.empty(interior_shape, dtype=bool)
        batch_rho_pow2_core = np.empty(interior_shape, dtype=np.float32)
        time_indices_pre = np.arange(prealloc_size)[:, None]
        profile_activity_map = bool(params.get("profile_activity_map", False))
        if profile_activity_map:
            t_fill_scatter = 0.0
            t_dilate = 0.0
            t_localmax = 0.0
            t_accumulate = 0.0

        for batch_idx in tqdm(range(num_batches), desc="Creating activity map"):
            # Calculate batch bounds
            batch_start = batch_idx * batch_size 
            batch_end = min(batch_start + batch_size, n_frames)
            padded_start = max(0, batch_start - temporal_pad)
            padded_end = min(n_frames, batch_end + temporal_pad)
            curr_size = padded_end - padded_start
            
            # Views into preallocated buffers for the current batch
            br = batch_rho[:curr_size]
            nm = nan_mask[:curr_size]

            # Fill batch data and record which voxels were written.
            # rho is (pixels, time), so batch over axis 1 and transpose to (time, pixels).
            # Initializing nm=True marks unwritten voxels as invalid without scanning full br.
            if profile_activity_map:
                t0 = time.perf_counter()
            br.fill(0)
            nm.fill(True)
            time_indices = time_indices_pre[:curr_size]
            batch_vals = rho[:, padded_start:padded_end].T
            br[time_indices, depth_indices, row_indices_roi, col_indices_roi] = batch_vals
            nm[time_indices, depth_indices, row_indices_roi, col_indices_roi] = np.isnan(batch_vals)
            if profile_activity_map:
                t_fill_scatter += time.perf_counter() - t0

            # Handle NaN values
            if profile_activity_map:
                t0 = time.perf_counter()
            br[nm] = 0
            core_t = curr_size - 2
            if core_t <= 0:
                continue

            dilated_nan_mask = fast_dilation(nm, dilation_struct)
            if profile_activity_map:
                t_dilate += time.perf_counter() - t0

            if profile_activity_map:
                t0 = time.perf_counter()
            center = br[1:-1, :, 1:-1, 1:-1]
            lmc = local_maxima_core[:core_t]
            ctmp = compare_tmp_core[:core_t]
            b3c = batch_rho_pow2_core[:core_t]

            np.greater(center, br[:-2, :, 1:-1, 1:-1], out=lmc)
            np.greater(center, br[2:, :, 1:-1, 1:-1], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            np.greater(center, br[1:-1, :, :-2, 1:-1], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            np.greater(center, br[1:-1, :, 2:, 1:-1], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            np.greater(center, br[1:-1, :, 1:-1, :-2], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            np.greater(center, br[1:-1, :, 1:-1, 2:], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            np.logical_not(dilated_nan_mask[1:-1, :, 1:-1, 1:-1], out=ctmp)
            np.logical_and(lmc, ctmp, out=lmc)
            if profile_activity_map:
                t_localmax += time.perf_counter() - t0

            if profile_activity_map:
                t0 = time.perf_counter()
            np.multiply(center, center, out=b3c)
            np.multiply(b3c, lmc, out=b3c, casting='unsafe')
            act_im[:, row_start + 1:row_stop - 1, col_start + 1:col_stop - 1] += np.sum(b3c, axis=0, dtype=np.float32)
            if profile_activity_map:
                t_accumulate += time.perf_counter() - t0

        # del rho, br, batch_rho
        if profile_activity_map:
            total_profiled = t_fill_scatter + t_dilate + t_localmax + t_accumulate
            if total_profiled > 0:
                print(
                    "[profile_activity_map] "
                    f"fill_scatter={t_fill_scatter:.3f}s ({100*t_fill_scatter/total_profiled:.1f}%), "
                    f"dilate={t_dilate:.3f}s ({100*t_dilate/total_profiled:.1f}%), "
                    f"localmax={t_localmax:.3f}s ({100*t_localmax/total_profiled:.1f}%), "
                    f"accumulate={t_accumulate:.3f}s ({100*t_accumulate/total_profiled:.1f}%)"
                )

        ##
        nan_mask = np.full_like(act_im, True, dtype=bool)
        validSelPix = np.flatnonzero(nanCt <= 0.5)
        nan_mask[np.unravel_index(selPixIdxs[validSelPix], nan_mask.shape)] = False
 
        # nan_mask = (act_im == 0) | (np.isnan(act_im))
        act_im[nan_mask] = np.nan

        med_act_im = ndimage.generic_filter(act_im, np.nanmedian, size=(1,11,11))
        act_im = act_im - med_act_im
        act_im[nan_mask] = np.nan
        
        # act_im[nan_mask] = np.nanmedian(act_im.flatten())
        # act_im_filt = ndimage.gaussian_filter(act_im, sigma=[0, 0.5, 0.5])
        # act_im_filt[nan_mask] = np.nan
        # act_im[nan_mask] = np.nan

        # act_im = act_im_filt

        # act_im_hp = act_im - ndimage.gaussian_filter(act_im, sigma=[10*i for i in psf_shape])

        excl_mask = None
        if params['select_soma'] and (f'DMD{DMDix+1}' in soma_masks):
            excl_mask = soma_masks[f'DMD{DMDix+1}'] > 0

        source_seeds = get_act_im_peaks(act_im, peak_th=params.get('peakth', 3.0),
                                        exclusion_mask=excl_mask)

        nSources = source_seeds.shape[0]

        def sel_pix_gaussian_profile(gaussian_params):
            z_planes = gaussian_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]
            y_means = gaussian_params[:, 1].unsqueeze(0)  # Shape: [1, nSources]
            x_means = gaussian_params[:, 2].unsqueeze(0)  # Shape: [1, nSources]
            y_sigmas = gaussian_params[:, 3].unsqueeze(0)  # Shape: [1, nSources]
            x_sigmas = gaussian_params[:, 4].unsqueeze(0)  # Shape: [1, nSources]
            corr_coef = torch.tanh(gaussian_params[:, 5].unsqueeze(0))  # Shape: [1, nSources]

            # Center the coordinates
            y_centered = (pixel_coords_tensor[:, 1].unsqueeze(1) - y_means)  # Shape: [len(selPixIdxs), nSources]
            x_centered = (pixel_coords_tensor[:, 2].unsqueeze(1) - x_means)  # Shape: [len(selPixIdxs), nSources]

            # Compute terms for bivariate Gaussian with correlation
            z_score_y = y_centered / y_sigmas
            z_score_x = x_centered / x_sigmas
            
            # Full bivariate Gaussian formula with correlation
            exponent = (-1 / (2 * (1 - corr_coef**2))) * (
                z_score_y**2 - 
                2 * corr_coef * z_score_x * z_score_y + 
                z_score_x**2
            )

            profile = torch.exp(exponent)
            profile = profile * ((torch.sqrt(z_score_y**2 + z_score_x**2) <= 3) & (pixel_coords_tensor[:, 0].unsqueeze(1) == z_planes)).float()

            return profile / torch.sum(profile,dim=0,keepdim=True)

        def sel_pix_patch_profile(patch_params):
            z_planes = patch_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]
            y_means = patch_params[:, 1].unsqueeze(0)  # Shape: [1, nSources]
            x_means = patch_params[:, 2].unsqueeze(0)  # Shape: [1, nSources]
            y_radii = patch_params[:, 3].unsqueeze(0)  # Shape: [1, nSources]
            x_radii = patch_params[:, 4].unsqueeze(0)  # Shape: [1, nSources]

            # Center the coordinates
            y_centered = (pixel_coords_tensor[:, 1].unsqueeze(1) - y_means)  # Shape: [len(selPixIdxs), nSources]
            x_centered = (pixel_coords_tensor[:, 2].unsqueeze(1) - x_means)  # Shape: [len(selPixIdxs), nSources]

            profile = (torch.abs(y_centered) < y_radii) & (torch.abs(x_centered) < x_radii) & (pixel_coords_tensor[:, 0].unsqueeze(1) == z_planes)

            return profile
        
        source_params = torch.cat([
            torch.tensor(source_seeds, dtype=torch.float32), # z, y, x means
            torch.ones(nSources, 2, dtype=torch.float32), # y, x sigmas
            torch.zeros(nSources, 1, dtype=torch.float32) # invtanh of correlation / tilt
        ], dim=1)

        A_patches = torch.zeros((nPixels, nSources), dtype=torch.bool)
        A_patches[selPixIdxs,:] = sel_pix_patch_profile(torch.cat([torch.tensor(source_seeds, dtype=torch.float32),
                                                            params['dXY'] * torch.ones(nSources, 2, dtype=torch.float32)],
                                                            dim=1))

        X_support_mots = [None] * uniqueMotionToKeepYX.shape[0]
        for i in range(uniqueMotionToKeepYX.shape[0]):
            X_support_mots[i] = torch.sparse.mm(H_mots[i], A_patches[selPixIdxs,:].float()) > 0

        # initialize sources as gaussian blobs
        A = torch.zeros((nPixels, nSources), dtype=torch.float32)
        A[selPixIdxs,:] = sel_pix_gaussian_profile(source_params * torch.tensor([1, 1, 1, 1, 1, 1]))
        A[~A_patches] = 0
        sources_total_mass = torch.sum(A, dim=0, keepdim=True)
        sources_total_mass[sources_total_mass <= 0] = 1
        A /= sources_total_mass

        # NMF parameters
        outer_loop_iters = 10
        
        mult_nmf_max_iters = 10
        nmf_tol = 1e-6

        # Gaussian fitting optimization parameters
        learning_rate = 0.01
        num_epochs = uniqueMotionToKeepYX.shape[0] * 5
        gd_tol = 1e-4
        
        # phi_lowRes = torch.zeros(lowResData.shape[1], nSources+1, dtype=torch.float32)
        # phi_lowRes[:] = np.nan

        phi_lowRes = torch.full((lowResDataNorm.shape[1], nSources), np.nan, dtype=torch.float32)

        X_mots = [None] * uniqueMotionToKeepYX.shape[0]
        overall_losses = [0] * (outer_loop_iters+1)

        data_for_nmf = torch.from_numpy(residual.astype(np.float32, copy=False))
        # data_for_nmf = torch.from_numpy(signal.convolve2d(lowResDataNorm,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same').astype(np.float32))
        # background_spatial_components = torch.from_numpy(background_spatial_components.astype(np.float32)) / torch.norm(background_spatial_components)
        
        for outer_loop_iter in range(outer_loop_iters): # outer loop

            # get phi_lowRes for current spatial profiles
            for i in range(uniqueMotionToKeepYX.shape[0]):
                motion_idx = i
                motion_frames = (motIndsYX == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # print warning for blank spatial components
                norms = torch.norm(X,dim=0,keepdim=True)
                if torch.any(norms == 0):
                    print(f"Warning: {np.nonzero(norms.numpy().squeeze() == 0)[0]} norms are zero for motion {i}")

                # least squares to fit phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T

                loss_contribution = torch.sum((data_for_nmf[:, motion_frames] - X @ phi_lowRes[motion_frames,:].T) ** 2).item()
                # print(f"Loss contribution for motion {motion_idx}: {loss_contribution}")

                overall_losses[outer_loop_iter] += loss_contribution
            print(f"Overall loss for outer loop iteration {outer_loop_iter}: {overall_losses[outer_loop_iter]}")

            shuffled_indices = torch.randperm(uniqueMotionToKeepYX.shape[0])
            for idx in tqdm(shuffled_indices,desc=f'Multiplicative NMF per motion displacement'): # loop over all motion displacements in random order
                i = idx.item()
                motion_idx = i
                motion_frames = (motIndsYX == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                X[~X_support_mots[i]] = 0

                # add background spatial component
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # normalize spatial components
                X = X / torch.norm(X,dim=0,keepdim=True)

                # projected least squares to initialize phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]  # This gives all time points at once
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
                phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:], min=0)
            
                # Initialize variables for NMF
                error_values = []
                prev_reconstruction_error = float('inf')

                for iter_idx in range(mult_nmf_max_iters): # iterative fit in sp space with NMF
                    # Multiplicative update for NMF
                    # Update phi (temporal components) using multiplicative update rule
                    # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
                    numerator = X.T @ data_for_nmf[:, motion_frames]
                    denominator = (X.T @ X) @ phi_lowRes[motion_frames,:].T + 1e-10
                    phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:] * (numerator / denominator).T, min=0)
                    
                    # Update X (spatial components) using multiplicative update rule
                    # X = X * (data * phi^T) / (X * phi * phi^T + epsilon)
                    numerator = data_for_nmf[:, motion_frames] @ phi_lowRes[motion_frames,:]
                    denominator = X @ (phi_lowRes[motion_frames,:].T @ phi_lowRes[motion_frames,:]) + 1e-10
                    X = torch.clamp(X * (numerator / denominator), min=0)

                    # Apply spatial support constraint
                    X = X * X_support_mots[i].float()
                    # X[:,:-1] = X[:,:-1] * X_support_mots[i].float()
                    # X[:,-1] = background_spatial_components[:,i]
                    
                    # Normalize X and phi to avoid scaling ambiguity
                    norms = torch.norm(X, dim=0, keepdim=True)
                    X = torch.where(norms > 0, X / norms, X)
                    phi_lowRes[motion_frames,:] = torch.where(norms > 0, phi_lowRes[motion_frames,:] * norms, phi_lowRes[motion_frames,:])
                    
                    # Calculate reconstruction error
                    # reconstruction = X @ phi_lowRes[motion_frames,:].T
                    # current_error = torch.mean((data_for_nmf[:, motion_frames] - reconstruction)**2).item()
                    # error_values.append(current_error)
                    
                    # # Check for convergence
                    # if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    #     break
                        
                    # prev_reconstruction_error = current_error

                    if iter_idx % 3 == 0:   # sparsify step
                        X_max = torch.max(X, dim=0, keepdim=True)[0]
                        X = torch.where(X_max > 0, X / X_max, X)
                        X = torch.clamp(X, min=params['sparse_fac']) - params['sparse_fac']

                        # re-normalize X after sparsification
                        norms = torch.norm(X, dim=0, keepdim=True)
                        X = torch.where(norms > 0, X / norms, X)
                
                X_mots[i] = X
            
            # Initialize Adam optimizer
            optim_loc_params = source_params[:,1:3].clone().requires_grad_(True)
            optim_scale_params = source_params[:,3:5].clone().requires_grad_(True)
            optim_tilt_params = source_params[:,5].unsqueeze(1).clone().requires_grad_(True)
            optimizer = torch.optim.Adam([{'params': optim_loc_params, 'lr': 10 * learning_rate},
                                            {'params': optim_scale_params, 'lr': 0.1 * learning_rate},
                                            {'params': optim_tilt_params, 'lr': 0.1 * learning_rate}])
            
            # Gradient descent optimization for Gaussian parameters
            losses = []
            for epoch in tqdm(range(num_epochs),desc='Fitting Gaussian spatial profiles'): # gradient descent to fit pixel space profile to sp profile
                # Zero gradients
                optimizer.zero_grad()
                
                # calculate spatial profile
                A_step_sel_pix = sel_pix_gaussian_profile(torch.cat([source_params[:,0].unsqueeze(1), optim_loc_params, optim_scale_params, optim_tilt_params], dim=1))
                A_step_sel_pix = A_step_sel_pix * A_patches[selPixIdxs,:].float()
                sources_total_mass = torch.sum(A_step_sel_pix, dim=0, keepdim=True)
                sources_total_mass[sources_total_mass <= 0] = 1
                A_step_sel_pix /= sources_total_mass
                # A_step_sel_pix = torch.where(sources_total_mass > 0, A_step_sel_pix / sources_total_mass, A_step_sel_pix)
                
                if epoch % uniqueMotionToKeepYX.shape[0] == 0:
                    shuffled_indices = torch.randperm(uniqueMotionToKeepYX.shape[0])

                i = shuffled_indices[epoch % uniqueMotionToKeepYX.shape[0]].item()
                motion_idx = i

                # Compute superpixel spatial profile
                X_step = torch.sparse.mm(H_mots[i], A_step_sel_pix)
                norms = torch.norm(X_step, dim=0, keepdim=True)
                X_step = torch.where(norms > 0, X_step / norms, X_step)
            
                # Compute loss across all sources
                # loss = torch.sum((X_step - X_mots[i][:,:-1]) ** 2)
                loss = torch.sum((X_step - X_mots[i]) ** 2)

                loss.backward()
                optimizer.step()
                
                # Apply constraints after optimizer step
                with torch.no_grad():
                    # losses.append(loss.item())
                    optim_loc_params.clamp_(min=torch.from_numpy(source_seeds[:,1:3] - params['dXY']).float(), max=torch.from_numpy(source_seeds[:,1:3] + params['dXY']).float())
                    optim_scale_params[:,0].clamp_(min=0.3, max=5)
                    optim_scale_params[:,1].clamp_(min=0.3, max=5)

                if epoch % uniqueMotionToKeepYX.shape[0] == uniqueMotionToKeepYX.shape[0]-1:
                    total_loss = 0
                    for i in range(uniqueMotionToKeepYX.shape[0]):
                        X_step = torch.sparse.mm(H_mots[i], A_step_sel_pix)
                        norms = torch.norm(X_step, dim=0, keepdim=True)
                        X_step = torch.where(norms > 0, X_step / norms, X_step)

                        # total_loss += torch.sum((X_step - X_mots[i][:,:-1]) ** 2)
                        total_loss += torch.sum((X_step - X_mots[i]) ** 2)

                    # Check for convergence
                    losses.append(total_loss.item())
                    if epoch // uniqueMotionToKeepYX.shape[0] > 0 and abs(losses[-1] - losses[-2]) < gd_tol:
                        print(f"Converged at epoch {epoch+1} with loss: {loss.item():.8f}")
                        break

            # Update parameters
            source_params = torch.cat([source_params[:,0].unsqueeze(1), optim_loc_params.detach(), optim_scale_params.detach(), optim_tilt_params.detach()], dim=1)

            A = torch.zeros((nPixels, nSources), dtype=torch.float32)            
            A[selPixIdxs,:] = sel_pix_gaussian_profile(source_params * torch.tensor([1, 1, 1, 1, 1, 1]))
            A[~A_patches] = 0
            sources_total_mass = torch.sum(A, dim=0, keepdim=True)
            sources_total_mass[sources_total_mass <= 0] = 1
            A /= sources_total_mass
            # A = torch.where(sources_total_mass > 0, A / sources_total_mass, A)

            # get phi_lowRes for current spatial profiles
            for i in range(uniqueMotionToKeepYX.shape[0]):
                motion_idx = i
                motion_frames = (motIndsYX == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # print warning for blank spatial components
                norms = torch.norm(X,dim=0,keepdim=True)
                if torch.any(norms == 0):
                    print(f"Warning: {np.nonzero(norms.numpy().squeeze() == 0)[0]} norms are zero for motion {i}")

                # least squares to fit phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
            
            # sort sources by variance
            sortorder = np.argsort(-np.nansum((phi_lowRes[:,:nSources].numpy()-np.nanmean(phi_lowRes[:,:nSources].numpy(),axis=0))**2,axis=0))
            source_params = source_params[sortorder,:]
            source_seeds = source_seeds[sortorder,:]
            A = A[:,sortorder]
            A_patches = A_patches[:,sortorder]
            X_support_mots = [support[:,sortorder] for support in X_support_mots]
            # sortorder = np.append(sortorder,max(sortorder)+1)
            phi_lowRes = phi_lowRes[:,sortorder]

            if (outer_loop_iter+1) % 4 == 3: # prune sources
                residual = torch.zeros_like(data_for_nmf, dtype=torch.float32)
                residual[:] = np.nan
                for i in range(uniqueMotionToKeepYX.shape[0]):
                    motion_idx = i
                    motion_frames = (motIndsYX == motion_idx).nonzero()[0]
                    X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                    # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)
                    residual[:,motion_frames] = data_for_nmf[:,motion_frames] - X @ phi_lowRes[motion_frames,:].T

                varExpSource = torch.zeros(nSources, dtype=torch.float32)
                varResidual = torch.zeros(nSources, dtype=torch.float32)
                for j in range(nSources-1,-1,-1):
                    # find pixels that correspond to half mass of spatial profile
                    # sorted_vals = torch.sort(A[selPixIdxs,j], descending=True)[0]
                    # cumsum_vals = torch.cumsum(sorted_vals, dim=0)
                    # total_mass = cumsum_vals[-1]
                    # half_mass_idx = torch.searchsorted(cumsum_vals, total_mass * 0.5)
                    # thresh = sorted_vals[half_mass_idx]
                    # validPixs = (A[selPixIdxs,j] >= thresh).nonzero()[:,0]
                    # tmp_A = A[selPixIdxs,j].unsqueeze(1).clone()
                    # tmp_A[~validPixs] = 0

                    for i in range(uniqueMotionToKeepYX.shape[0]):
                        motion_idx = i
                        motion_frames = (motIndsYX == motion_idx).nonzero()[0]

                        # X = torch.sparse.mm(H_mots[i], tmp_A)
                        X = torch.sparse.mm(H_mots[i], A[selPixIdxs,j].unsqueeze(1))

                        contributingPixs = (X > 0).nonzero()[:,0]
                        varExpSource[j] += torch.sum(torch.sum(X[contributingPixs] * phi_lowRes[motion_frames,j].unsqueeze(-1).T,dim=0)**2)
                        varResidual[j] += torch.sum(torch.sum(residual[contributingPixs][:,motion_frames],dim=0)**2)
                
                SNR = varExpSource / varResidual
                
                SNR_cut = 1/3
                keepSources = (SNR > SNR_cut).nonzero()[:,0]
                print(f"Keeping {keepSources.shape[0]} of {nSources} sources")
                nSources = keepSources.shape[0]
                source_params = source_params[keepSources,:]
                source_seeds = source_seeds[keepSources,:]
                A = A[:,keepSources]
                A_patches = A_patches[:,keepSources]
                X_support_mots = [support[:,keepSources] for support in X_support_mots]
                # keepSources = np.append(keepSources,phi_lowRes.shape[1]-1)
                phi_lowRes = phi_lowRes[:,keepSources]
        
        # get phi_lowRes for current spatial profiles
        for i in range(uniqueMotionToKeepYX.shape[0]):
            motion_idx = i
            motion_frames = (motIndsYX == motion_idx).nonzero()[0]

            # project image space (A) into superpixel space (X)
            X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
            # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

            overall_losses[outer_loop_iters] += torch.sum((data_for_nmf[:, motion_frames] - X @ phi_lowRes[motion_frames,:].T) ** 2).item()
        print(f"Final overall loss: {overall_losses[outer_loop_iters]}")

        trial_info = [(i, keepTrials[DMDix,i], background[:,lowResTrialID == i]) for i in range(nTrials)]

        get_high_res_traces_partial = partial(get_high_res_traces,
                                              DMDix=DMDix,
                                              params=params,
                                              sampFreq=params['analyzeHz'],
                                              refStack=refStack,
                                              subsampleMatrixInds=subsampleMatrixInds[f'DMD{DMDix+1}'],
                                              fastZ2RefZ=fastZ2RefZ,
                                              sparseHInds=sparseHInds,
                                              sparseHVals=sparseHVals,
                                              allSuperPixelIDs=allSuperPixelIDs,
                                              dr=dr, 
                                              trialTable=trialTable, 
                                              A_final=A,
                                              uniqueMotionDS=uniqueMotionToKeepYX,
                                              motIndsToKeepDS=np.arange(uniqueMotionToKeepYX.shape[0]),
                                              median_z=median_z,
                                              psf=psf,
                                              soma_sps=soma_sps[f'DMD{DMDix+1}'])

        with mp.Pool(processes=min(params['max_workers'],min(mp.cpu_count(),len(trial_info)))) as pool:
            results = list(pool.imap(get_high_res_traces_partial, trial_info))

        # Save data to HDF5 file
        with h5py.File(output_h5_filename, 'a') as f:
            # Delete group if it exists
            group_name = f'DMD{DMDix+1}'
            if group_name in f:
                del f[group_name]
            
            # Create group and add datasets
            dmd_group = f.create_group(group_name)
            
            # Create subgroups for trial data
            source_group = dmd_group.create_group('sources')

            spatial_group = source_group.create_group('spatial')

            # spatial_group.create_dataset('source_params', data=source_params.numpy())
            spatial_group.create_dataset('profiles', data=A.numpy().reshape(dmdPixelsPerRow,dmdPixelsPerColumn,-1).transpose(2,0,1))
            spatial_group.create_dataset('coords', data=source_params[:,:3].numpy())

            temporal_group = source_group.create_group('temporal')

            dF = np.concatenate([r[0] for r in results], axis=0)
            F0 = np.concatenate([r[1] for r in results], axis=0)
            F = dF + F0

            F0 = compute_f0(F, int(np.ceil(params['denoiseWindow_s']*params['analyzeHz'])), int(np.ceil(params['baselineWindow_s']*params['analyzeHz'])))
            dF = F - F0
            dFF = dF / np.clip(F0, 1e-4, None)
            temporal_group.create_dataset('dF', data=dF)
            temporal_group.create_dataset('dFF', data=dFF)
            temporal_group.create_dataset('F0', data=F0)
            # temporal_group.create_dataset('F', data=F)

            # trial_start_idxs = np.concatenate([[0], np.cumsum([len(r[2]) for r in results])[:-1]])
            trial_num_frames = np.concatenate([[len(r[2])] for r in results])

            frame_group = dmd_group.create_group('frame_info')
            # frame_group.create_dataset('trial_start_idxs', data=trial_start_idxs)
            frame_group.create_dataset('trial_num_frames', data=trial_num_frames)
            frame_group.create_dataset('discard_frames', data=np.any(np.isnan(F), axis=1))
            frame_group.create_dataset('frame_line_idxs', data=np.concatenate([r[2] for r in results], axis=0))

            frame_group.create_dataset('offlineXshifts', data=np.concatenate([r[5][1] for r in results], axis=0))
            frame_group.create_dataset('offlineYshifts', data=np.concatenate([r[5][0] for r in results], axis=0))
            frame_group.create_dataset('offlineZshifts', data=np.concatenate([r[5][2] for r in results], axis=0))

            frame_group.create_dataset('onlineXshifts', data=np.concatenate([r[6][1] for r in results], axis=0))
            frame_group.create_dataset('onlineYshifts', data=np.concatenate([r[6][0] for r in results], axis=0))
            frame_group.create_dataset('onlineZshifts', data=np.concatenate([r[6][2] for r in results], axis=0))

            # frame_group.create_dataset('selPixIdxs', data=[r[2] for r in results])

            globalF = np.concatenate([r[4] for r in results], axis=0)
            global_group = dmd_group.create_group('global')
            global_group.create_dataset('F', data=globalF)

            visualizations = dmd_group.create_group('visualizations')
            visualizations.create_dataset('act_im', data=act_im)
            visualizations.create_dataset('mean_im', data=np.expand_dims(mean_im, axis=0))
            
            if params['select_soma'] and (f'DMD{DMDix+1}' in soma_masks):
                user_roi_group = dmd_group.create_group('user_rois')
                masks_by_roi = np.zeros((numFastZs, *yx_shape, np.max(soma_masks[f'DMD{DMDix+1}'])), dtype=bool)
                for roi in range(masks_by_roi.shape[3]):
                    masks_by_roi[:,:,:,roi] = soma_masks[f'DMD{DMDix+1}'] == roi + 1
                user_roi_group.create_dataset('mask', data=masks_by_roi)
                user_roi_group.create_dataset('F', data=np.concatenate([r[7] for r in results], axis=0))

            ref_stack_channels = trialTable['refStack'][0,DMDix]['channels'][0,0].squeeze().reshape((-1,))
            num_ref_stack_channels = len(ref_stack_channels)
            
            ref_stack = trialTable['refStack'][0,DMDix]['IM'][0,0].T
            ref_stack_reshaped = np.zeros((num_ref_stack_channels,ref_stack.shape[0] // num_ref_stack_channels,ref_stack.shape[1], ref_stack.shape[2]), dtype=np.float32)
            for i in range(ref_stack.shape[0]):
                c = i % num_ref_stack_channels
                ref_stack_reshaped[c,i//num_ref_stack_channels,:,:] = ref_stack[i,:,:]
            del ref_stack
            
            ref_stack_dataset = visualizations.create_dataset('ref_stack', data=ref_stack_reshaped)
            ref_stack_dataset.attrs['channels'] = ref_stack_channels

        print(f"Added DMD{DMDix+1} data to {output_h5_filename}")

    end_time = datetime.now(timezone.utc).astimezone()

    try:
        version = subprocess.check_output(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"], text=True
        ).strip()
    except:
        version = "Unknown"

    p = Processing.create_with_sequential_process_graph(
        pipelines=[
            Code(
                name="SLAP2 band scanning processing pipeline (Python)",
                url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/",
                version=version,
            ),
        ],
        data_processes=[
            DataProcess(
                process_type=ProcessName.VIDEO_MOTION_CORRECTION,
                pipeline_name="SLAP2 band scanning processing pipeline (Python)",
                experimenters=[align_params.get('operator', 'Unknown')],
                stage=ProcessStage.PROCESSING,
                start_date_time=datetime.fromisoformat(align_params.get('startTime', datetime.now(timezone.utc).astimezone().isoformat())),
                end_date_time=datetime.fromisoformat(align_params.get('endTime', datetime.now(timezone.utc).astimezone().isoformat())),
                output_path="..",
                code=Code(
                    url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/blob/main/matlab/preprocessing/integrationRegistration.m",
                    version=version,
                    parameters=align_params
                ),
            ),
            DataProcess(
                process_type=ProcessName.VIDEO_ROI_TIMESERIES_EXTRACTION,
                experimenters=[params.get('operator', 'Unknown')],
                stage=ProcessStage.PROCESSING,
                start_date_time=start_time,
                end_date_time=end_time,
                output_path=".",
                pipeline_name="SLAP2 band scanning processing pipeline (Python)",
                code=Code(
                    url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/blob/main/python/slap2analysis/extractSLAP2IntegrationSources.py",
                    version=version,
                    parameters={key: to_serializable(val) for key, val in params.items()}
                ),
            ),
        ],
    )

    serialized = p.model_dump_json()
    deserialized = Processing.model_validate_json(serialized)
    p.write_standard_file(Path(savedr),suffix='_integration.json')

if __name__ == '__main__':
    main()