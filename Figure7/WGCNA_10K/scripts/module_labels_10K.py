#!/usr/bin/env python3
"""Single source of truth for module names and colours across every figure.

WGCNA colour names are assigned by module size and carry no meaning; the same
word means different gene sets in the two networks (pooled brown rises with
ageing, DAB1 brown falls). Network names are provenance, not biology: the
"DAB1 network" is the one built from the DAB1-arm samples, and several of its
modules are not DAB1-responsive.

So figures use CONTENT labels, with provenance kept for the legend only:
    N39            pooled network, all 39 samples, batch-corrected
    N25            DAB1-arm samples (n=25), batch + APOE genotype regressed out

Colours are assigned by biological family, so related modules across the two
networks share a hue and can be tracked by eye:
    vascular / ECM / development   teal-green
    lysosomal / DAM / TREM2        purple-magenta
    inflammatory NF-kB / IFN       red-orange
    translation / metabolic        blue
    chromatin / cell cycle         brown-grey

LABEL[(network, module)] -> (content label, family, colour, provenance tag)
"""

FAMILY_COLOUR = {
    "vascular":     "#159AA8",
    "ecm":          "#2E9B3C",
    "lysosomal":    "#7E3FAE",
    "dam":          "#B5479B",
    "inflammatory": "#D6453B",
    "interferon":   "#E8823C",
    "translation":  "#2C6FB3",
    "metabolic":    "#6C8FB8",
    "chromatin":    "#8A6A4F",
    "other":        "#9A9A9A",
}

_L = [
    # network, module, content label, family, provenance
    ("pooled", "turquoise",  "Vascular / tube morphogenesis",     "vascular",     "N39"),
    ("pooled", "purple",     "Lysosomal–phagocytic · DAM/TREM2",  "lysosomal",    "N39"),
    ("pooled", "salmon",     "Innate immune–lysosomal · HM",      "lysosomal",    "N39"),
    ("pooled", "greenyellow","Lysosomal–lipid · DAM",             "dam",          "N39"),
    ("pooled", "red",        "Cation transport; TREM2-LAM · DAM", "dam",          "N39"),
    ("pooled", "brown",      "TNFα–NFκB / metabolic · CRM",       "inflammatory", "N39"),
    ("pooled", "green",      "TNFα–NFκB / circadian",             "inflammatory", "N39"),
    ("pooled", "pink",       "Inflammatory / NOD-like · CRM",     "inflammatory", "N39"),
    ("pooled", "magenta",    "Complement / innate immune · HM",   "inflammatory", "N39"),
    ("pooled", "black",      "Antiviral–IFN · IRM",               "interferon",   "N39"),
    ("pooled", "tan",        "Antigen presentation · HLA",        "interferon",   "N39"),
    ("pooled", "blue",       "mTORC1 / cholesterol · RM/HLA",     "metabolic",    "N39"),
    ("pooled", "yellow",     "Nucleosome / chromatin",            "chromatin",    "N39"),
    ("pooled", "cyan",       "Unannotated",                       "other",        "N39"),
    ("DAB1",   "blue",       "Lysosomal–innate immune · DAM/HM",  "lysosomal",    "N25"),
    ("DAB1",   "turquoise",  "IFN / TNFα–NFκB · IRM/CRM",         "interferon",   "N25"),
    ("DAB1",   "pink",       "Innate immune / neutrophil",        "inflammatory", "N25"),
    ("DAB1",   "red",        "Translation / ribosome · RM",       "translation",  "N25"),
    ("DAB1",   "green",      "ECM / EMT",                         "ecm",          "N25"),
    ("DAB1",   "yellow",     "ECM / EMT",                         "ecm",          "N25"),
    ("DAB1",   "brown",      "L1CAM / development",               "vascular",     "N25"),
    ("DAB1",   "black",      "Cell cycle",                        "chromatin",    "N25"),
    ("DAB1",   "magenta",    "Glycosylation / collagen",          "ecm",          "N25"),
]

# True WGCNA colour, so swatches match the panel c/d bar charts exactly.
WGCNA_COLOUR = {"turquoise":"#1FB6B6","blue":"#2C6FB3","brown":"#7F4F2D","yellow":"#C9A800",
    "green":"#2E9B3C","red":"#D62728","black":"#333333","pink":"#E377A8","magenta":"#C026B3",
    "purple":"#7E3FAE","greenyellow":"#6F8F00","tan":"#B08A4F","salmon":"#E8735E","cyan":"#159AA8"}

LABEL = {(n, m): (lab, fam, FAMILY_COLOUR[fam], prov) for n, m, lab, fam, prov in _L}

def label(nw, m):   return LABEL.get((nw, m), (m, "other", FAMILY_COLOUR["other"], "?"))[0]
def colour(nw, m):  return WGCNA_COLOUR.get(m, "#999999")   # true WGCNA colour
def famcolour(nw, m): return LABEL.get((nw, m), (m, "other", FAMILY_COLOUR["other"], "?"))[2]
def family(nw, m):  return LABEL.get((nw, m), (m, "other", FAMILY_COLOUR["other"], "?"))[1]
def prov(nw, m):    return LABEL.get((nw, m), (m, "other", FAMILY_COLOUR["other"], "?"))[3]
def tag(nw, m):
    """Short identifier for a figure row: provenance + original colour."""
    return f"{prov(nw, m)}-{m}"

if __name__ == "__main__":
    for (n, m), (lab, fam, col, pv) in LABEL.items():
        print(f"{pv}-{m:12s} {fam:13s} {col}  {lab}")
