#!/usr/bin/env python3
"""RELN signalling schematic v4 — integrin arm: ITGB1 + ITGB8 with DAB2 adaptor (DAB1 removed)."""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
CM, DPI = 1/2.54, 600
FONTS = {"title": 15, "row": 13, "value": 11}
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","Liberation Sans","DejaVu Sans"],"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
COL_BLUE_DK="#1F4E79"; COL_BLUE_LT="#D6E4F0"; COL_GREEN_DK="#1F7A4C"; COL_GREEN_LT="#D6EBDE"
COL_ORANGE_DK="#C46A1F"; COL_ORANGE_LT="#FBE4D0"; COL_APOE="#B2182B"; COL_RELN="#1F4E79"
COL_RELN_GR="#2A6F4E"; COL_MEMBR="#F2D17C"; COL_PURPLE_DK="#6A3D9A"; COL_PURPLE_LT="#E6DCF2"
fig_w_cm, fig_h_cm = 23.0, 14.0
fig = plt.figure(figsize=(fig_w_cm*CM, fig_h_cm*CM))
ax = fig.add_axes([0,0,1,1]); ax.set_xlim(0,fig_w_cm); ax.set_ylim(0,fig_h_cm); ax.axis("off")
MEM_Y, MEM_H = 8.2, 0.9
ax.add_patch(Rectangle((0.4,MEM_Y), fig_w_cm-0.8, MEM_H, facecolor=COL_MEMBR, edgecolor="#A0834A", linewidth=0.6, zorder=2))
REC_H, REC_W = 1.3, 1.8; REC_Y = MEM_Y - 0.4
def receptor(x,label,sublabel=None,fill=COL_BLUE_DK,edge="#222",linestyle="-",w=REC_W,h=REC_H,label_color="white",zorder=4):
    ax.add_patch(FancyBboxPatch((x,REC_Y),w,h,boxstyle="round,pad=0.02,rounding_size=0.12",facecolor=fill,edgecolor=edge,linewidth=1.0,linestyle=linestyle,zorder=zorder))
    y_lab=REC_Y+h/2+(0.18 if sublabel else 0.0)
    ax.text(x+w/2,y_lab,label,ha="center",va="center",fontsize=FONTS["row"]-1,weight="bold",color=label_color,zorder=10)
    if sublabel: ax.text(x+w/2,y_lab-0.50,sublabel,ha="center",va="center",fontsize=FONTS["value"]-0.5,color=label_color,style="italic",zorder=10)
def dab_circle(cx,cy,label,fill,r=0.55,fs=None,zorder=5):
    ax.add_patch(Circle((cx,cy),r,facecolor=fill,edgecolor="#222",linewidth=0.7,zorder=zorder))
    ax.text(cx,cy+0.02,label,ha="center",va="center",fontsize=fs if fs else FONTS["row"]-1,weight="bold",color="white",zorder=zorder+5)
def block_header(x,w,text,color): ax.text(x+w/2,MEM_Y+MEM_H+0.5,text,ha="center",va="bottom",fontsize=FONTS["row"],weight="bold",color=color,zorder=10)
ADAPT_Y = MEM_Y - 1.55
# Canonical (blue)
CAN_X, CAN_W = 1.0, 4.4
receptor(CAN_X+0.2,"VLDLR",fill=COL_BLUE_DK); receptor(CAN_X+2.4,"LRP8",sublabel="ApoER2",fill=COL_BLUE_DK)
block_header(CAN_X,CAN_W,"Canonical",COL_BLUE_DK)
DAB1_X = CAN_X+CAN_W/2+0.35
dab_circle(DAB1_X,ADAPT_Y,"DAB1",COL_BLUE_DK)
ax.text(DAB1_X+0.72,ADAPT_Y-0.02,"pY",ha="left",va="center",fontsize=FONTS["value"]-0.5,style="italic",color="#222",zorder=12)
D2C_X = DAB1_X+1.70
dab_circle(D2C_X,ADAPT_Y,"DAB2",COL_ORANGE_DK,r=0.38,fs=FONTS["value"]-3.0)
ax.add_patch(FancyArrowPatch((CAN_X+2.4+REC_W/2,REC_Y),(D2C_X-0.05,ADAPT_Y+0.30),arrowstyle="-|>",mutation_scale=11,color=COL_ORANGE_DK,linewidth=1.2,linestyle=(0,(3,2.5)),connectionstyle="arc3,rad=-0.30",zorder=3))
# ===== Integrin / adhesion (green) — ITGB1 + ITGB8, DAB2 adaptor (DAB1 removed) =====
INT_X, INT_W = 7.4, 4.4; INT_CX = INT_X+INT_W/2
ITGB1_X = INT_CX-REC_W/2
receptor(ITGB1_X,"ITGB1",sublabel="α3",fill=COL_GREEN_DK)
block_header(INT_X,INT_W,"Integrin / adhesion",COL_GREEN_DK)
DAB2I_X = INT_CX
dab_circle(DAB2I_X,ADAPT_Y,"DAB2",COL_ORANGE_DK)
ITGB1_CX = ITGB1_X+REC_W/2
# Lipoprotein / APOE (orange)
LIP_X, LIP_W = 13.8, 7.2
GHOST_X=14.1; LRP1_X,LRP2_X = LIP_X+2.4, LIP_X+4.8
block_header(LIP_X,LIP_W,"Lipoprotein / APOE",COL_ORANGE_DK)
receptor(GHOST_X,"LRP8",fill=COL_PURPLE_LT,edge=COL_PURPLE_DK,linestyle="--",label_color=COL_PURPLE_DK)
receptor(LRP1_X,"LRP1",fill=COL_ORANGE_DK); receptor(LRP2_X,"LRP2",fill=COL_ORANGE_DK)
DAB2_X = LIP_X+LIP_W/2+0.3
dab_circle(DAB2_X,ADAPT_Y,"DAB2",COL_ORANGE_DK)
def tether(x0,y0,x1,y1,col,lw=1.5,ls="-"): ax.plot([x0,x1],[y0,y1],color=col,linewidth=lw,linestyle=ls,zorder=3,solid_capstyle="round")
VLDLR_CX=CAN_X+0.2+REC_W/2; LRP8c_CX=CAN_X+2.4+REC_W/2; LRP1_CX=LRP1_X+REC_W/2; LRP2_CX=LRP2_X+REC_W/2
tether(DAB1_X,ADAPT_Y,VLDLR_CX,REC_Y,COL_BLUE_DK)
tether(DAB1_X,ADAPT_Y,LRP8c_CX,REC_Y,COL_BLUE_DK)
tether(DAB2I_X,ADAPT_Y,ITGB1_CX,REC_Y,COL_ORANGE_DK)   # DAB2 -> ITGB1 (beta1-integrin endocytosis; Teckchandani 2012)
tether(DAB2_X,ADAPT_Y,LRP1_CX,REC_Y,COL_ORANGE_DK); tether(DAB2_X,ADAPT_Y,LRP2_CX,REC_Y,COL_ORANGE_DK)
# Ligand labels
LIG_Y=12.3; APOE_X=CAN_X+CAN_W/2
ax.text(APOE_X,LIG_Y,"APOE4",ha="center",va="center",color=COL_APOE,fontsize=FONTS["title"],weight="bold",zorder=10)
APOE_TIP_Y=MEM_Y+MEM_H+1.70
ax.plot([APOE_X,APOE_X],[LIG_Y-0.45,APOE_TIP_Y+0.10],color=COL_APOE,linewidth=1.8,zorder=3,solid_capstyle="round")
BAR_HW=0.55
ax.plot([APOE_X-BAR_HW,APOE_X+BAR_HW],[APOE_TIP_Y+0.10,APOE_TIP_Y+0.10],color=COL_APOE,linewidth=2.6,zorder=3,solid_capstyle="round")
ax.text(APOE_X-BAR_HW-0.10,(LIG_Y+APOE_TIP_Y)/2+0.10,"competes with\nRELN",ha="right",va="center",fontsize=FONTS["value"]-1.5,color=COL_APOE,style="italic",zorder=10)
RELN_X=(CAN_X+CAN_W+INT_X)/2
ax.text(RELN_X,LIG_Y,"RELN",ha="center",va="center",color=COL_RELN,fontsize=FONTS["title"],weight="bold",zorder=10)
ax.add_patch(FancyArrowPatch((RELN_X-0.5,LIG_Y-0.45),(CAN_X+CAN_W-0.6,MEM_Y+MEM_H+0.85),arrowstyle="-|>",mutation_scale=12,color=COL_RELN,linewidth=1.5,zorder=3))
ax.add_patch(FancyArrowPatch((RELN_X+0.5,LIG_Y-0.45),(INT_CX,MEM_Y+MEM_H+1.30),arrowstyle="-|>",mutation_scale=12,color=COL_RELN_GR,linewidth=1.5,zorder=3))
APOE2_X=LIP_X+LIP_W/2+1.0
ax.text(APOE2_X,LIG_Y,"APOE",ha="center",va="center",color=COL_APOE,fontsize=FONTS["title"],weight="bold",zorder=10)
ax.add_patch(FancyArrowPatch((APOE2_X,LIG_Y-0.45),((LRP1_X+LRP2_X+REC_W)/2,MEM_Y+MEM_H+0.85),arrowstyle="-|>",mutation_scale=14,color=COL_APOE,linewidth=1.6,zorder=3))
ax.add_patch(FancyArrowPatch((RELN_X+0.6,LIG_Y-0.45),(GHOST_X+REC_W/2,MEM_Y+MEM_H+0.85),arrowstyle="-|>",mutation_scale=12,color=COL_RELN,linewidth=1.3,linestyle=(0,(3,2.5)),connectionstyle="arc3,rad=-0.35",zorder=3))
ax.add_patch(FancyArrowPatch((GHOST_X+REC_W/2,REC_Y),(DAB2_X-0.45,ADAPT_Y+0.5),arrowstyle="-|>",mutation_scale=11,color=COL_PURPLE_DK,linewidth=1.2,linestyle=(0,(3,2.5)),connectionstyle="arc3,rad=0.15",zorder=3))
D2E_X=GHOST_X+REC_W/2-1.15
dab_circle(D2E_X,ADAPT_Y,"DAB2",COL_ORANGE_DK,r=0.38,fs=FONTS["value"]-3.0)
ax.add_patch(FancyArrowPatch((GHOST_X+REC_W/2,REC_Y),(D2E_X+0.05,ADAPT_Y+0.30),arrowstyle="-|>",mutation_scale=11,color=COL_ORANGE_DK,linewidth=1.2,linestyle=(0,(3,2.5)),connectionstyle="arc3,rad=0.30",zorder=3))
# Downstream
DS_Y_1, DS_H = 5.3, 0.9; DS_Y_2, DS_H2 = 3.55, 2.55
def downstream_box(x,w,y,h,fill_lt,fill_dk): ax.add_patch(FancyBboxPatch((x,y-h),w,h,boxstyle="round,pad=0.02,rounding_size=0.08",facecolor=fill_lt,edgecolor=fill_dk,linewidth=0.9,zorder=4))
def box_title(cx,y,h,txt,color): ax.text(cx,y-h/2,txt,ha="center",va="center",fontsize=FONTS["row"]-1,weight="bold",color=color,zorder=10)
def stacked(cx,y_top,lines,fs=None):
    fs=fs if fs else FONTS["value"]-0.5; dy=0.52
    for i,(txt,col) in enumerate(lines): ax.text(cx,y_top-0.35-i*dy,txt,ha="center",va="center",fontsize=fs,weight="bold",color=col,zorder=10)
def connector(x,y0,y1,col,ls="-"): ax.add_patch(FancyArrowPatch((x,y0),(x,y1),arrowstyle="-|>",mutation_scale=10,color=col,linewidth=1.2,linestyle=ls,zorder=3))
CBX=CAN_X-0.4; CBW=CAN_W+0.8; CBcx=CBX+CBW/2
SF_cx,SF_w=1.85,2.5; CK_cx,CK_w=4.45,2.3
ax.add_patch(FancyArrowPatch((DAB1_X,ADAPT_Y-0.55),(SF_cx,DS_Y_1+0.05),arrowstyle="-|>",mutation_scale=10,color=COL_BLUE_DK,linewidth=1.2,zorder=3))
ax.add_patch(FancyArrowPatch((DAB1_X,ADAPT_Y-0.55),(CK_cx,DS_Y_1+0.05),arrowstyle="-|>",mutation_scale=10,color=COL_BLUE_DK,linewidth=1.2,zorder=3))
downstream_box(SF_cx-SF_w/2,SF_w,DS_Y_1,DS_H,COL_BLUE_LT,COL_BLUE_DK); box_title(SF_cx,DS_Y_1,DS_H,"SRC/FYN →\nPI3K-AKT",COL_BLUE_DK)
downstream_box(CK_cx-CK_w/2,CK_w,DS_Y_1,DS_H,COL_BLUE_LT,COL_BLUE_DK); box_title(CK_cx,DS_Y_1,DS_H,"CRKL / SFK",COL_BLUE_DK)
downstream_box(CBX,CBW,DS_Y_2,DS_H2,"#FFFFFF",COL_BLUE_DK)
stacked(CBcx,DS_Y_2,[("Cytoskeleton","#222"),("Motility","#222"),("p-Tau phos.","#222"),("ABCA1 / lipids","#222"),("↓ LRP8",COL_ORANGE_DK)])
ax.add_patch(FancyArrowPatch((SF_cx,DS_Y_1-DS_H),(SF_cx,DS_Y_2+0.04),arrowstyle="-|>",mutation_scale=9,color=COL_BLUE_DK,linewidth=1.0,zorder=3))
ax.add_patch(FancyArrowPatch((CK_cx,DS_Y_1-DS_H),(CK_cx,DS_Y_2+0.04),arrowstyle="-|>",mutation_scale=9,color=COL_BLUE_DK,linewidth=1.0,zorder=3))
# Integrin downstream (now from DAB2)
connector(DAB2I_X,ADAPT_Y-0.55,DS_Y_1+0.05,COL_GREEN_DK)
IBX=INT_X-0.4; IBW=INT_W+0.8; IBcx=IBX+IBW/2
downstream_box(IBX,IBW,DS_Y_1,DS_H,COL_GREEN_LT,COL_GREEN_DK); box_title(IBcx,DS_Y_1,DS_H,"Focal-adhesion\ndisassembly",COL_GREEN_DK)
downstream_box(IBX,IBW,DS_Y_2,DS_H2,"#FFFFFF",COL_GREEN_DK)
stacked(IBcx,DS_Y_2-0.55,[("integrin internalisation",COL_ORANGE_DK),("↓ surface β1 integrin",COL_ORANGE_DK),("↓ FAK / adhesion",COL_GREEN_DK),("altered migration",COL_GREEN_DK)])
connector(IBcx,DS_Y_1-DS_H,DS_Y_2+0.04,COL_GREEN_DK)
# Lipoprotein downstream
CLA_cx,CLA_w=19.3,4.2; PI3_cx,PI3_w=15.3,3.2
ax.add_patch(FancyArrowPatch((DAB2_X,ADAPT_Y-0.55),(CLA_cx,DS_Y_1+0.05),arrowstyle="-|>",mutation_scale=10,color=COL_ORANGE_DK,linewidth=1.2,zorder=3))
ax.add_patch(FancyArrowPatch((DAB2_X,ADAPT_Y-0.55),(PI3_cx,DS_Y_1+0.05),arrowstyle="-|>",mutation_scale=10,color=COL_PURPLE_DK,linewidth=1.4,linestyle=(0,(4,3)),zorder=3))
downstream_box(CLA_cx-CLA_w/2,CLA_w,DS_Y_1,DS_H,COL_ORANGE_LT,COL_ORANGE_DK); box_title(CLA_cx,DS_Y_1,DS_H,"Clathrin\nendocytosis",COL_ORANGE_DK)
downstream_box(CLA_cx-CLA_w/2,CLA_w,DS_Y_2,DS_H2,"#FFFFFF",COL_ORANGE_DK)
stacked(CLA_cx,DS_Y_2,[("APOE / Aβ uptake","#222"),("Cholesterol","#222"),("LDLR / NPC1L1 sorting",COL_ORANGE_DK),("lipid handling","#222"),("↓ LRP1 / LRP2",COL_ORANGE_DK)])
ax.add_patch(FancyArrowPatch((CLA_cx,DS_Y_1-DS_H),(CLA_cx,DS_Y_2+0.04),arrowstyle="-|>",mutation_scale=9,color=COL_ORANGE_DK,linewidth=1.0,zorder=3))
downstream_box(PI3_cx-PI3_w/2,PI3_w,DS_Y_1,DS_H,COL_PURPLE_LT,COL_PURPLE_DK); box_title(PI3_cx,DS_Y_1,DS_H,"PI3K-AKT",COL_PURPLE_DK)
downstream_box(PI3_cx-PI3_w/2,PI3_w,DS_Y_2,DS_H2,"#FFFFFF",COL_PURPLE_DK)
stacked(PI3_cx,DS_Y_2,[("NF-κB",COL_PURPLE_DK),("↑ VCAM1 / ICAM1",COL_PURPLE_DK),("pro-inflammatory",COL_PURPLE_DK),("↓ LRP8",COL_ORANGE_DK)])
ax.add_patch(FancyArrowPatch((PI3_cx,DS_Y_1-DS_H),(PI3_cx,DS_Y_2+0.04),arrowstyle="-|>",mutation_scale=9,color=COL_PURPLE_DK,linewidth=1.0,zorder=3))
fig.savefig("RELN_signalling_schematic_v4fix.png",dpi=DPI,bbox_inches="tight",facecolor="white")
print("saved")
