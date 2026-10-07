"""
_shared.style — shared style helpers for Claude Python renderers.

Three things are standardised across every figure:

  1.  Font family   : Arial / Helvetica / DejaVu Sans (fallbacks)
  2.  Font sizes    : 5 cm-readable
        title 14 | row label 13 | col label 12 | cell value 9 | star 12
        cbar label 12 | cbar tick 10
  3.  Diverging palette (NES heatmaps) : deeper blue-white-red
        #053061 → #FFFFFF → #67001F, clipped to ±3.

The module is a drop-in for any new script:
    from _shared.style import apply_style, NES_CMAP, NES_NORM, FONTS, stars
"""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

CM = 1 / 2.54
DPI = 600

FONTS = {"title": 14, "value": 20, "star": 18, "row": 18, "col": 14, "caption": 9, "cbar_lbl": 13, "cbar_tick": 11, "group": 14, "mod": 14, "leg": 11}


def apply_style(base: float = 10.0):
    """Configure matplotlib defaults (Arial, linewidths, savefig)."""
    mpl.rcParams.update({
        "font.family":     "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size":       base,
        "axes.linewidth":  0.6,
        "savefig.dpi":     DPI,
        "savefig.bbox":    "tight",
        "pdf.fonttype":    42,
        "ps.fonttype":     42,
    })


# Diverging blue-white-red palette — fixed to user's reference image.
# Do NOT change without explicit user instruction.
NES_CMAP = LinearSegmentedColormap.from_list(
    "div_bwr", ["#2166AC", "#FFFFFF", "#B2182B"], N=256
)
NES_NORM = TwoSlopeNorm(vcenter=0.0, vmin=-2, vmax=2)


def stars(p) -> str:
    if p is None or (isinstance(p, float) and p != p):
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""
