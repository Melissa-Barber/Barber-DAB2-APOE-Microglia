"""NES-space PCA v2 (collapsed MG panel) — SAME AESTHETICS as the original
four_separate biplots, applied to the corrected feature set:
   • Haney2024_LDAM_markers DROPPED as a feature (circularity: it is also a
     contrast source) -> replaced by independent Marschallinger2020_LDAM_up
   • LXR_cholesterol_efflux_MG ADDED
   • RELN axis = single RELN_DAB1 / RELN_DAB2 (the two used throughout MS)
   • same 28 contrasts (APOE2 excluded); scaled (correlation) PCA

Outputs (NEW filenames; nothing old overwritten):
   Fig_NES_PCA_biplot_v2collapsed.png            (PC1 x PC2)
   Fig_NES_PCA_biplot_PC1_PC3_v2collapsed.png
   Fig_NES_PCA_biplot_PC2_PC3_v2collapsed.png
   Fig_NES_PCA_PC1_loadings_v2collapsed.png  (+PC2,PC3)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import FancyArrowPatch, Ellipse
from adjustText import adjust_text

FGSEA   = Path(sys.argv[1])
OUT_BASE = Path(sys.argv[2])
(OUT_BASE / "figures").mkdir(parents=True, exist_ok=True)
(OUT_BASE / "data").mkdir(parents=True, exist_ok=True)

CM, DPI = 1/2.54, 600
F_TITLE, F_SUB, F_PANEL, F_LABEL, F_TICK, F_LEG, F_LBL_PT, F_ARROW, F_BAR = \
    16, 11, 14, 13, 11, 11, 9, 11, 12

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial","Helvetica","Liberation Sans","DejaVu Sans"],
    "font.size": 11, "axes.linewidth": 0.7,
    "savefig.dpi": DPI, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# ── feature set (corrected) ──────────────────────────────────────────────
MG_SETS = ["MG_HM","MG_DAM","MG_HLA","MG_CRM","MG_IRM",
           "Marschallinger2020_LDAM_up","Victor2022_APOE4_lipid_iMG",
           "LXR_cholesterol_efflux_MG"]
RELN_SETS = ["RELN_DAB1","RELN_DAB2"]
ALL_SETS = RELN_SETS + MG_SETS

GROUP_MAP = {
    "BIONMG_DAB1_KD_vs_Ctrl_E3":"BIONMG KD","BIONMG_DAB1_KD_vs_Ctrl_E4":"BIONMG KD",
    "BIONMG_DAB1_KD_vs_Ctrl_pooled":"BIONMG KD","BIONMG_DAB2_KD_vs_Ctrl_E3":"BIONMG KD",
    "BIONMG_DAB2_KD_vs_Ctrl_E4":"BIONMG KD","BIONMG_DAB2_KD_vs_Ctrl_pooled":"BIONMG KD",
    "BIONMG_E4_vs_E3_Ctrl":"BIONMG genotype","BIONMG_LRT_Genotype_x_siRNA":"BIONMG genotype",
    "Coray_E44_vs_E33_LDneg":"iPSC-MG APOE","Coray_E44_vs_E33_LDpos":"iPSC-MG APOE",
    "Coray_LDpos_vs_LDneg_E33":"iPSC-MG APOE","Coray_LDpos_vs_LDneg_E44":"iPSC-MG APOE",
    "MIC_E4_vs_E3":"iPSC-MG APOE","BitBio_E4_vs_E3":"iPSC-MG APOE",
    "Coray_snMG_AD33_vs_Ctrl":"Coray AD snRNA","Coray_snMG_AD44_vs_AD33":"Coray AD snRNA",
    "Coray_snMG_AD44_vs_Ctrl":"Coray AD snRNA","MG_ADN_vs_Ctrl":"AD Progression Atlas",
    "MG_EC_PathAdv_vs_Early":"AD Progression Atlas","MG_ITG_PathAdv_vs_Early":"AD Progression Atlas",
    "MG_PFC_PathAdv_vs_Early":"AD Progression Atlas","MG_V1_PathAdv_vs_Early":"AD Progression Atlas",
    "MG_V2_PathAdv_vs_Early":"AD Progression Atlas","Galatro_aged_vs_young":"Developmental MG",
    "Kracht_GW_early_vs_late":"Developmental MG","Kracht_GW_late_vs_early":"Developmental MG",
    "Ostrem_P14_vs_E13":"Developmental MG","Ostrem_P14_vs_E15":"Developmental MG",
}
LABEL_MAP = {
    "BIONMG_DAB1_KD_vs_Ctrl_E3":"DAB1 KD E3","BIONMG_DAB1_KD_vs_Ctrl_E4":"DAB1 KD E4",
    "BIONMG_DAB1_KD_vs_Ctrl_pooled":"DAB1 KD pool","BIONMG_DAB2_KD_vs_Ctrl_E3":"DAB2 KD E3",
    "BIONMG_DAB2_KD_vs_Ctrl_E4":"DAB2 KD E4","BIONMG_DAB2_KD_vs_Ctrl_pooled":"DAB2 KD pool",
    "BIONMG_E4_vs_E3_Ctrl":"E4 vs E3","BIONMG_LRT_Genotype_x_siRNA":"LRT geno×siRNA",
    "Coray_E44_vs_E33_LDneg":"E44 vs E33 LD−","Coray_E44_vs_E33_LDpos":"E44 vs E33 LD+",
    "Coray_LDpos_vs_LDneg_E33":"LD+ vs LD− E33","Coray_LDpos_vs_LDneg_E44":"LD+ vs LD− E44",
    "MIC_E4_vs_E3":"MIC E4 vs E3","BitBio_E4_vs_E3":"BitBio E4 vs E3",
    "Coray_snMG_AD33_vs_Ctrl":"snMG AD33","Coray_snMG_AD44_vs_AD33":"snMG AD44 vs 33",
    "Coray_snMG_AD44_vs_Ctrl":"snMG AD44","MG_ADN_vs_Ctrl":"MG ADN",
    "MG_EC_PathAdv_vs_Early":"EC PathAdv","MG_ITG_PathAdv_vs_Early":"ITG PathAdv",
    "MG_PFC_PathAdv_vs_Early":"PFC PathAdv","MG_V1_PathAdv_vs_Early":"V1 PathAdv",
    "MG_V2_PathAdv_vs_Early":"V2 PathAdv","Galatro_aged_vs_young":"Aged vs young",
    "Kracht_GW_early_vs_late":"GW early vs late","Kracht_GW_late_vs_early":"GW late vs early",
    "Ostrem_P14_vs_E13":"P14 vs E13","Ostrem_P14_vs_E15":"P14 vs E15",
}
GROUP_COLORS = {"AD Progression Atlas":"#D6604D","Coray AD snRNA":"#B2182B",
                "iPSC-MG APOE":"#F4A582","BIONMG genotype":"#74ADD1",
                "BIONMG KD":"#2166AC","Developmental MG":"#4DAC26"}
GROUP_MARKERS = {"AD Progression Atlas":"o","Coray AD snRNA":"^","iPSC-MG APOE":"s",
                 "BIONMG genotype":"D","BIONMG KD":"o","Developmental MG":"v"}

# ── PCA (scaled / correlation; NA->0 to match v2 outputs) ─────────────────
fg = pd.read_csv(FGSEA, low_memory=False)
sub = fg[fg["contrast_id"].isin(GROUP_MAP) & fg["pathway"].isin(ALL_SETS)][
    ["contrast_id","pathway","NES"]]
wide = sub.groupby(["contrast_id","pathway"])["NES"].mean().unstack("pathway")
wide["RELN→DAB1"] = wide["RELN_DAB1"]
wide["RELN→DAB2"] = wide["RELN_DAB2"]

PCA_FEATURES = ["RELN→DAB1","RELN→DAB2"] + MG_SETS
X = wide[PCA_FEATURES].fillna(0.0)
Xs = StandardScaler().fit_transform(X.values)
pca = PCA(n_components=3, random_state=0)
scores = pca.fit_transform(Xs)
pct = 100 * pca.explained_variance_ratio_

scores_df = pd.DataFrame(scores, columns=["PC1","PC2","PC3"], index=X.index)
scores_df["Group"] = scores_df.index.map(GROUP_MAP)
scores_df["label"] = [LABEL_MAP.get(i, i) for i in scores_df.index]
loadings = pd.DataFrame(pca.components_.T, columns=["PC1","PC2","PC3"], index=PCA_FEATURES)
print(f"PC1={pct[0]:.1f}%  PC2={pct[1]:.1f}%  PC3={pct[2]:.1f}%")

FEAT_LABEL = {"RELN→DAB1":"RELN→DAB1","RELN→DAB2":"RELN→DAB2",
              "MG_HM":"HM","MG_DAM":"DAM","MG_HLA":"HLA","MG_CRM":"CRM","MG_IRM":"IRM",
              "Marschallinger2020_LDAM_up":"LDAM","Victor2022_APOE4_lipid_iMG":"APOE4-lipid",
              "LXR_cholesterol_efflux_MG":"LXR-chol"}
ARROW_COLOR = {"RELN":"#2A6F4D","MG":"#8B1A1A"}
def cat_for(f): return "RELN" if f in ("RELN→DAB1","RELN→DAB2") else "MG"
DISPLAY = {"RELN→DAB1":"RELN∩DAB1","RELN→DAB2":"RELN∩DAB2",
           "MG_HM":"MG_HM","MG_DAM":"MG_DAM","MG_HLA":"MG_HLA","MG_CRM":"MG_CRM","MG_IRM":"MG_IRM",
           "Marschallinger2020_LDAM_up":"Marschallinger LDAM",
           "Victor2022_APOE4_lipid_iMG":"Victor2022 APOE4-lipid",
           "LXR_cholesterol_efflux_MG":"LXR cholesterol-efflux"}

def render_biplot(pc_x, pc_y, out_name, title_suffix="", label_overrides=None,
                  dab1_label_offset=(0.45,-1.30), dab2_label_offset=(0.45,0.75),
                  ad_progression_offset=(2.6,1.2), ad_progression_scales=(0.85,0.55),
                  apoe_theta_deg=-25.0, apoe_length_scale=1.95,
                  show_poles=True, show_apoe=True, show_hypothesis=True):
    fig, ax = plt.subplots(figsize=(22.0*CM, 22.0*CM))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.91, bottom=0.13)
    ax.axhline(0, color="#bbb", lw=0.6, linestyle="--", zorder=0)
    ax.axvline(0, color="#bbb", lw=0.6, linestyle="--", zorder=0)

    sx, sy = scores_df[pc_x].values, scores_df[pc_y].values
    score_range = max(abs(sx).max(), abs(sy).max()) * 0.55
    load_norm = np.sqrt(loadings[pc_x]**2 + loadings[pc_y]**2).max()
    scale = score_range / load_norm
    xmin, xmax = sx.min()-2.4, sx.max()+3.4
    ymin, ymax = sy.min()-2.0, sy.max()+2.6

    ax.add_patch(plt.Rectangle((xmin,0), -xmin, ymax, color="#E3EBF5", alpha=0.45, zorder=-2, ec="none"))
    ax.add_patch(plt.Rectangle((0,0), xmax, ymax, color="#FCE6E2", alpha=0.45, zorder=-2, ec="none"))
    ax.add_patch(plt.Rectangle((xmin,ymin), -xmin, -ymin, color="#F2F2F2", alpha=0.50, zorder=-2, ec="none"))
    ax.add_patch(plt.Rectangle((0,ymin), xmax, -ymin, color="#FFF5E0", alpha=0.45, zorder=-2, ec="none"))

    def add_group_ellipse(groups, color, edge, label, label_offset=(0,0),
                          width_scale=1.0, height_scale=1.0):
        sel = scores_df[scores_df["Group"].isin(groups)]
        if len(sel) < 2: return
        xy = sel[[pc_x, pc_y]].values; mu = xy.mean(axis=0)
        if len(sel) >= 3:
            cov = np.cov(xy.T); ev, evec = np.linalg.eigh(cov)
            order = ev.argsort()[::-1]; ev, evec = ev[order], evec[:,order]
            angle = np.degrees(np.arctan2(*evec[:,0][::-1]))
            width, height = 2*2.0*np.sqrt(np.clip(ev, 0.05, None))
        else:
            d = np.linalg.norm(xy[1]-xy[0]); width, height = max(1.5, d*1.3), 1.2
            angle = np.degrees(np.arctan2(*(xy[1]-xy[0])[::-1]))
        width *= width_scale; height *= height_scale
        ax.add_patch(Ellipse(xy=mu, width=width, height=height, angle=angle,
                             facecolor=color, alpha=0.16, edgecolor=edge,
                             linewidth=1.6, linestyle="--", zorder=-1))
        ax.text(mu[0]+label_offset[0], mu[1]+height/2*1.05+label_offset[1], label,
                ha="center", va="bottom", color=edge, style="italic",
                fontsize=F_BAR+1, weight="bold", zorder=10,
                bbox=dict(facecolor="white", edgecolor=edge, boxstyle="round,pad=0.30",
                          linewidth=0.9, alpha=0.95))

    add_group_ellipse(["AD Progression Atlas","Coray AD snRNA"], "#D7263D", "#B2182B",
                      "AD progression", label_offset=ad_progression_offset,
                      width_scale=ad_progression_scales[0], height_scale=ad_progression_scales[1])

    def add_contrast_ellipse(substr, color, edge, label, label_dx=0.0, label_dy=-0.7):
        mask = scores_df.index.to_series().str.contains(substr)
        sel = scores_df[mask]
        if len(sel) < 2: return
        xy = sel[[pc_x, pc_y]].values; mu = xy.mean(axis=0)
        if len(sel) >= 3:
            cov = np.cov(xy.T); ev, evec = np.linalg.eigh(cov)
            order = ev.argsort()[::-1]; ev = np.clip(ev[order],0.05,None); evec = evec[:,order]
            angle = np.degrees(np.arctan2(*evec[:,0][::-1])); width, height = 2*2.0*np.sqrt(ev)
        else:
            d = np.linalg.norm(xy[1]-xy[0]); width, height = max(1.5, d*1.3), 1.2
            angle = np.degrees(np.arctan2(*(xy[1]-xy[0])[::-1]))
        ax.add_patch(Ellipse(xy=mu, width=width, height=height, angle=angle,
                             facecolor=color, alpha=0.18, edgecolor=edge,
                             linewidth=1.6, linestyle="--", zorder=-1))
        ax.text(mu[0]+label_dx, mu[1]+label_dy, label, ha="center", va="center",
                color=edge, style="italic", fontsize=F_BAR+1, weight="bold", zorder=10,
                bbox=dict(facecolor="white", edgecolor=edge, boxstyle="round,pad=0.30",
                          linewidth=0.9, alpha=0.95))

    add_contrast_ellipse("BIONMG_DAB2_KD_vs_Ctrl", "#5BA0CF", "#1F4E79", "DAB2-KD",
                         label_dx=dab2_label_offset[0], label_dy=dab2_label_offset[1])
    add_contrast_ellipse("BIONMG_DAB1_KD_vs_Ctrl", "#9DBED9", "#3A6FA0", "DAB1-KD",
                         label_dx=dab1_label_offset[0], label_dy=dab1_label_offset[1])

    inner = 0.04
    if show_poles:
        ax.text(xmin+(xmax-xmin)*inner, ymax-(ymax-ymin)*inner, "DAB-KD pole",
                ha="left", va="top", fontsize=F_BAR, color="#1F4E79", style="italic", weight="bold")
        ax.text(xmax-(xmax-xmin)*inner, ymax-(ymax-ymin)*inner, "AD / RELN-active pole",
                ha="right", va="top", fontsize=F_BAR, color="#B2182B", style="italic", weight="bold")
        ax.text(xmin+(xmax-xmin)*inner, ymin+(ymax-ymin)*inner, "Developmental / aging pole",
                ha="left", va="bottom", fontsize=F_BAR, color="#555", style="italic", weight="bold")
        ax.text(xmax-(xmax-xmin)*inner, ymin+(ymax-ymin)*inner, "Homeostatic / DAM-family pole",
                ha="right", va="bottom", fontsize=F_BAR, color="#8C7B00", style="italic", weight="bold")

    DEFAULT_LABEL_OVERRIDES = {
        "Victor2022_APOE4_lipid_iMG": ("above", 0.0, 0.55),
        "RELN→DAB2": ("above", -0.4, 1.10),
        "RELN→DAB1": ("above", -0.4, 0.55),
    }
    # Arrows drawn first; feature labels collected for automatic de-collision
    # (adjustText) with leader lines tethering each label to its arrow tip.
    feat_texts, tip_x, tip_y = [], [], []
    # RELN labels are pinned OFF the arrow shaft (perpendicular, up-left) with a
    # leader line, and excluded from the auto-repel so they never sit on the
    # green arrows. Different magnitudes stagger DAB1 vs DAB2.
    RELN_OFF = {"RELN→DAB2": 1.55, "RELN→DAB1": 0.85}
    for feat, row in loadings.iterrows():
        x1, y1 = row[pc_x]*scale, row[pc_y]*scale
        col = ARROW_COLOR[cat_for(feat)]
        lw = 2.6 if cat_for(feat) == "RELN" else 1.6
        ax.add_patch(FancyArrowPatch((0,0), (x1,y1), arrowstyle="-|>", mutation_scale=16,
                                     lw=lw, color=col, zorder=4))
        if feat in RELN_OFF:
            nrm = np.hypot(x1, y1) or 1.0
            perpx, perpy = -y1/nrm, x1/nrm
            if perpx > 0: perpx, perpy = -perpx, -perpy   # bias to upper-left
            off = RELN_OFF[feat]
            ax.annotate(FEAT_LABEL.get(feat, feat), xy=(x1, y1),
                        xytext=(x1 + perpx*off, y1 + perpy*off),
                        ha="center", va="center", color=col, fontsize=F_ARROW,
                        fontweight="bold", zorder=7,
                        bbox=dict(facecolor="white", edgecolor=col, boxstyle="round,pad=0.20",
                                  linewidth=0.7, alpha=0.96),
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.6, alpha=0.85))
        else:
            feat_texts.append(ax.text(x1*1.08, y1*1.08, FEAT_LABEL.get(feat, feat),
                    ha="center", va="center", color=col, fontsize=F_ARROW, fontweight="bold",
                    zorder=6, bbox=dict(facecolor="white", edgecolor=col, boxstyle="round,pad=0.20",
                                        linewidth=0.7, alpha=0.96)))
            tip_x.append(x1); tip_y.append(y1)

    apoe_row = loadings.loc["Victor2022_APOE4_lipid_iMG"]
    base_x, base_y = apoe_row[pc_x]*scale, apoe_row[pc_y]*scale
    base_len = np.sqrt(base_x**2 + base_y**2)
    target_len = base_len * apoe_length_scale
    theta = np.radians(apoe_theta_deg)
    dirx = (base_x*np.cos(theta) - base_y*np.sin(theta)) / base_len
    diry = (base_x*np.sin(theta) + base_y*np.cos(theta)) / base_len
    ax_apoe_x, ax_apoe_y = dirx*target_len, diry*target_len
    if show_apoe:
        ax.add_patch(FancyArrowPatch((0,0), (ax_apoe_x, ax_apoe_y), arrowstyle="-|>",
                                     mutation_scale=24, lw=3.0, color="#E08531", alpha=0.92, zorder=3.5))
        ax.text(ax_apoe_x*1.22, ax_apoe_y*1.22, "APOE4", ha="center", va="center", color="#A05810",
                fontsize=F_ARROW+2, fontweight="bold", zorder=6,
                bbox=dict(facecolor="white", edgecolor="#E08531", boxstyle="round,pad=0.30",
                          linewidth=1.0, alpha=0.95))

    kd_mask = scores_df.index.to_series().str.contains("BIONMG_DAB[12]_KD")
    kd_sel = scores_df[kd_mask]
    if show_hypothesis and len(kd_sel) >= 2:
        kd_x, kd_y = kd_sel[[pc_x, pc_y]].mean().values
        mid_x, mid_y = (kd_x+ax_apoe_x)/2, (kd_y+ax_apoe_y)/2
        pull = 0.30
        start_x, start_y = kd_x+(mid_x-kd_x)*pull, kd_y+(mid_y-kd_y)*pull
        end_x, end_y = ax_apoe_x+(mid_x-ax_apoe_x)*pull, ax_apoe_y+(mid_y-ax_apoe_y)*pull
        ax.add_patch(FancyArrowPatch((start_x,start_y), (end_x,end_y),
                                     connectionstyle="arc3,rad=0.30", arrowstyle="<->",
                                     mutation_scale=18, linestyle=(0,(4,3)), color="#444",
                                     lw=1.8, alpha=0.85, shrinkA=2, shrinkB=2, zorder=4))
        mx, my = (start_x+end_x)/2, (start_y+end_y)/2
        dx, dy = end_x-start_x, end_y-start_y
        chord = np.sqrt(dx*dx+dy*dy); nx, ny = dy/chord, -dx/chord
        ax.text(mx+nx*0.18*chord, my+ny*0.18*chord, "hypothesis", ha="center", va="center",
                color="#444", fontsize=F_ARROW, style="italic", zorder=6,
                bbox=dict(facecolor="white", edgecolor="#888", boxstyle="round,pad=0.25",
                          linewidth=0.7, alpha=0.95))

    for grp, color in GROUP_COLORS.items():
        sel = scores_df[scores_df["Group"] == grp]
        if sel.empty: continue
        ax.scatter(sel[pc_x], sel[pc_y], s=160, marker=GROUP_MARKERS[grp], facecolor=color,
                   edgecolor="white", linewidth=1.5, alpha=0.95, zorder=3, label=grp)

    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)

    # Auto de-collide the feature labels: repel from each other and from all
    # data points, keep a thin leader line back to each arrow tip.
    px = list(scores_df[pc_x].values); py = list(scores_df[pc_y].values)
    adjust_text(feat_texts, x=px, y=py, target_x=tip_x, target_y=tip_y, ax=ax,
                expand=(1.35, 1.6), force_text=(0.5, 0.9), force_static=(0.25, 0.4),
                force_pull=(0.02, 0.02), ensure_inside_axes=True, min_arrow_len=6,
                arrowprops=dict(arrowstyle="-", color="0.45", lw=0.6, alpha=0.85))

    pc_idx = {"PC1":0, "PC2":1, "PC3":2}
    pcx_var, pcy_var = pct[pc_idx[pc_x]], pct[pc_idx[pc_y]]
    def disp(feats): return " · ".join(FEAT_LABEL.get(f,f) for f in feats)
    px_pos = loadings.sort_values(pc_x, ascending=False).head(2)
    px_neg = loadings.sort_values(pc_x, ascending=True).head(2)
    py_pos = loadings.sort_values(pc_y, ascending=False).head(2)
    py_neg = loadings.sort_values(pc_y, ascending=True).head(2)
    ax.set_xlabel(f"← {disp(px_neg.index)}                 {pc_x} ({pcx_var:.1f}% variance)"
                  f"                 {disp(px_pos.index)} →", fontsize=F_LABEL, fontweight="bold")
    ax.set_ylabel(f"← {disp(py_neg.index)}                 {pc_y} ({pcy_var:.1f}% variance)"
                  f"                 {disp(py_pos.index)} →", fontsize=F_LABEL, fontweight="bold")
    ax.tick_params(axis="both", labelsize=F_TICK)
    fig.suptitle(f"Microglial state landscape:  {pc_x} × {pc_y}{title_suffix}",
                 fontsize=F_TITLE, fontweight="bold", y=0.965)
    fig.text(0.5, 0.925, f"n = {len(scores_df)} contrasts (APOE2 excl.; Haney dropped as feature)  |  "
             f"{pc_x} = {pcx_var:.1f}%   {pc_y} = {pcy_var:.1f}%", ha="center",
             fontsize=F_SUB, style="italic", color="#444")
    handles = [plt.Line2D([],[], marker=GROUP_MARKERS[g], color="w", markerfacecolor=c,
                          markeredgecolor="white", markeredgewidth=1.2, markersize=12,
                          linestyle="", label=g) for g,c in GROUP_COLORS.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5,-0.08), ncol=3,
              frameon=False, fontsize=F_LEG, handletextpad=0.4, columnspacing=1.4)
    out = OUT_BASE/"figures"/out_name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print(f"saved {out}")

def render_loading(pc, var_pct, out_name, title):
    df = loadings.reset_index().rename(columns={"index":"feature"})
    df["display"] = df["feature"].map(DISPLAY)
    df = df.sort_values(pc, ascending=True).reset_index(drop=True)
    POS, NEG = "#B2182B", "#2166AC"
    colors = [POS if v >= 0 else NEG for v in df[pc]]
    fig, ax = plt.subplots(figsize=(15.0*CM, 14.5*CM))
    fig.subplots_adjust(left=0.36, right=0.94, top=0.85, bottom=0.22)
    ax.barh(np.arange(len(df)), df[pc], color=colors, edgecolor="white", linewidth=0.8, height=0.78)
    ax.axvline(0, color="#222", lw=0.7)
    ax.set_yticks(np.arange(len(df))); ax.set_yticklabels(df["display"], fontsize=F_BAR+1)
    ax.tick_params(axis="x", labelsize=F_TICK)
    ax.set_xlabel(f"{pc} loading", fontsize=F_LABEL+1, fontweight="bold")
    ax.set_title(title, fontsize=F_PANEL, fontweight="bold", color="#111", pad=10, loc="left")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#999"); ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_color("#444"); ax.spines["bottom"].set_linewidth(0.7)
    leg_h = [plt.Rectangle((0,0),1,1,color=POS), plt.Rectangle((0,0),1,1,color=NEG)]
    ax.legend(leg_h, ["Up-regulated (+ loading)","Down-regulated (− loading)"],
              loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              frameon=False, fontsize=F_LEG, handletextpad=0.5, columnspacing=2.0)
    fig.suptitle(f"{pc} ({var_pct:.1f}% variance) — gene-set loadings",
                 fontsize=F_TITLE, fontweight="bold", y=0.97)
    out = OUT_BASE/"figures"/out_name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print(f"saved {out}")

render_biplot("PC1","PC2","Fig_NES_PCA_biplot_v2collapsed.png")
render_biplot("PC1","PC3","Fig_NES_PCA_biplot_PC1_PC3_v2collapsed.png","  (supplementary view)",
    label_overrides={
        "Victor2022_APOE4_lipid_iMG":("above",0.6,0.35),
        "RELN→DAB2":("above",-0.55,0.85), "RELN→DAB1":("above",-0.55,0.35),
        "Marschallinger2020_LDAM_up":("offset",0.55,-0.30),
        "LXR_cholesterol_efflux_MG":("offset",0.45,0.30),
        "MG_HM":("above",0.0,0.45), "MG_CRM":("above",-0.45,0.30),
        "MG_DAM":("offset",0.45,-0.20), "MG_HLA":("offset",0.55,-0.10), "MG_IRM":("offset",0.30,-0.25),
    },
    dab2_label_offset=(0.45,1.10), dab1_label_offset=(0.45,-1.10),
    ad_progression_offset=(2.8,0.9), apoe_theta_deg=-45.0, apoe_length_scale=2.10)
render_biplot("PC2","PC3","Fig_NES_PCA_biplot_PC2_PC3_v2collapsed.png","  (supplementary view)",
    label_overrides={
        "Victor2022_APOE4_lipid_iMG":("above",0.55,0.35),
        "RELN→DAB2":("above",0.0,0.55), "RELN→DAB1":("above",0.0,-0.55),
        "Marschallinger2020_LDAM_up":("offset",0.55,-0.10),
        "LXR_cholesterol_efflux_MG":("offset",0.45,0.30),
        "MG_HM":("above",0.0,0.45), "MG_CRM":("above",-0.55,0.35),
        "MG_DAM":("offset",0.45,-0.20), "MG_HLA":("offset",0.55,-0.10), "MG_IRM":("offset",-0.20,-0.35),
    },
    dab2_label_offset=(-1.5,0.7), dab1_label_offset=(-1.5,-0.8),
    ad_progression_offset=(1.5,0.5), ad_progression_scales=(0.95,0.55),
    apoe_theta_deg=-35.0, apoe_length_scale=1.55)
# Cleaned PC2×PC3: pole labels + APOE4 guidance arrow + hypothesis arc removed,
# leaving a de-cluttered correlation biplot (loading arrows + labels only).
render_biplot("PC2","PC3","Fig_NES_PCA_biplot_PC2_PC3_v2collapsed_CLEAN.png",
    "  (correlation view)",
    ad_progression_offset=(1.6,0.6), ad_progression_scales=(0.95,0.55),
    dab2_label_offset=(-1.6,0.8), dab1_label_offset=(-1.6,-0.9),
    show_poles=False, show_apoe=False, show_hypothesis=False)
render_loading("PC1", pct[0], "Fig_NES_PCA_PC1_loadings_v2collapsed.png", "PC1 loadings")
render_loading("PC2", pct[1], "Fig_NES_PCA_PC2_loadings_v2collapsed.png", "PC2 loadings")
render_loading("PC3", pct[2], "Fig_NES_PCA_PC3_loadings_v2collapsed.png", "PC3 loadings")
print("DONE")
