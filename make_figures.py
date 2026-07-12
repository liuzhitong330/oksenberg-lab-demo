"""
make_figures.py — Oksenberg Lab Demo
Figure 1: Dual bar chart sorted by Plp1/Bdnf ratio — myelin load vs neuroprotective
          capacity across MS-relevant brain regions. Highlights brainstem and cerebellum
          as having the most imbalanced myelin-to-BDNF ratio.
Figure 2: Schematic connecting the DYSF-ZNF638 severity locus (Oksenberg/IMSGC 2023 Nature)
          to the regional CNS resilience framework.
"""

import csv, os

OUT = os.path.dirname(os.path.abspath(__file__))

rows = []
with open(os.path.join(OUT, "ms_genes.tsv")) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

# ── Figure 1: Dual bar chart sorted by Plp1/Bdnf ratio ───────────────────────
sorted_rows = sorted(rows,
    key=lambda r: float(r["Plp1"]) / max(float(r["Bdnf"]), 1e-6),
    reverse=True)

FW, FH = 700, 420
PAD_L = 72; PAD_R = 155; PAD_T = 85; PAD_B = 55
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

N = len(sorted_rows)
row_h = AH / N
bar_gap = row_h * 0.18
bar_h = (row_h - bar_gap * 3) / 2

all_vals = [float(r["Plp1"]) for r in rows] + [float(r["Bdnf"]) for r in rows]
xmax = max(all_vals) * 1.1

def bx(v): return PAD_L + v / xmax * AW
def by_p(i): return PAD_T + i * row_h + bar_gap
def by_b(i): return PAD_T + i * row_h + bar_gap * 2 + bar_h

HIGH_IMPACT = {"BS", "Cereb", "Thal"}   # MS high-disability regions
COLORS = {
    "BS":    "#c0392b",
    "Cereb": "#c0392b",
    "Thal":  "#e67e22",
    "Ctx":   "#7f8c8d",
    "Hipp":  "#7f8c8d",
    "CP":    "#1a5c8a",
    "CC":    "#1a5c8a",
}

bars = ""
# Gridlines
for v in [0.04, 0.08, 0.12]:
    gx = bx(v)
    bars += (f'<line x1="{gx:.1f}" y1="{PAD_T}" x2="{gx:.1f}" y2="{PAD_T+AH}" '
             f'stroke="#eee" stroke-width="1"/>'
             f'<text x="{gx:.1f}" y="{PAD_T+AH+18}" text-anchor="middle" '
             f'font-size="9" fill="#bbb">{v:.2f}</text>')

for i, row in enumerate(sorted_rows):
    reg = row["region"]
    pv = float(row["Plp1"])
    bv = float(row["Bdnf"])
    ratio = pv / max(bv, 1e-6)
    col = COLORS.get(reg, "#95a5a6")
    hi = reg in HIGH_IMPACT

    if hi:
        bars += (f'<rect x="{PAD_L}" y="{PAD_T + i*row_h:.1f}" '
                 f'width="{AW}" height="{row_h:.1f}" fill="#fff5f5" rx="0"/>')

    bars += (f'<text x="{PAD_L-8}" y="{PAD_T + i*row_h + row_h/2 + 4:.1f}" '
             f'text-anchor="end" font-size="10.5" fill="{col}" '
             f'font-weight="{"700" if hi else "400"}">{reg}</text>')

    # Plp1 bar
    bars += (f'<rect x="{PAD_L}" y="{by_p(i):.1f}" width="{bx(pv)-PAD_L:.1f}" '
             f'height="{bar_h:.1f}" fill="{col}" opacity="0.85" rx="2"/>'
             f'<text x="{bx(pv)+4:.1f}" y="{by_p(i)+bar_h/2+3.5:.1f}" '
             f'font-size="8" fill="{col}">{pv:.3f}</text>')

    # Bdnf bar
    bars += (f'<rect x="{PAD_L}" y="{by_b(i):.1f}" width="{bx(bv)-PAD_L:.1f}" '
             f'height="{bar_h:.1f}" fill="{col}" opacity="0.35" rx="2"/>'
             f'<text x="{bx(bv)+4:.1f}" y="{by_b(i)+bar_h/2+3.5:.1f}" '
             f'font-size="8" fill="{col}">{bv:.4f}</text>')

    # Ratio on right
    ratio_x = PAD_L + AW + 14
    ratio_y = PAD_T + i * row_h + row_h / 2 + 4
    bars += (f'<text x="{ratio_x}" y="{ratio_y:.1f}" font-size="10" '
             f'fill="{col}" font-weight="{"700" if hi else "400"}">'
             f'{ratio:.0f}×</text>')

axes = (f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" '
        f'stroke="#ccc" stroke-width="1.2"/>'
        f'<line x1="{PAD_L}" y1="{PAD_T+AH}" x2="{PAD_L+AW}" y2="{PAD_T+AH}" '
        f'stroke="#ccc" stroke-width="1.2"/>'
        f'<text x="{PAD_L+AW/2:.0f}" y="{PAD_T+AH+36}" text-anchor="middle" '
        f'font-size="10" fill="#555">Expression density (Allen Mouse Brain Atlas ISH)</text>'
        f'<text x="{PAD_L+AW+14}" y="{PAD_T-10}" font-size="9" fill="#555" font-weight="700">'
        f'Plp1 / Bdnf</text>'
        f'<text x="{PAD_L+AW+14}" y="{PAD_T+0}" font-size="9" fill="#555">ratio →</text>')

leg_y = PAD_T + AH + 46
legend = (
    f'<rect x="{PAD_L}" y="{leg_y}" width="14" height="10" fill="#555" opacity="0.85" rx="1"/>'
    f'<text x="{PAD_L+18}" y="{leg_y+9}" font-size="9" fill="#555">Plp1 (myelin gene)</text>'
    f'<rect x="{PAD_L+135}" y="{leg_y}" width="14" height="10" fill="#555" opacity="0.35" rx="1"/>'
    f'<text x="{PAD_L+153}" y="{leg_y+9}" font-size="9" fill="#555">Bdnf (neuroprotection)</text>'
    f'<rect x="{PAD_L+330}" y="{leg_y-1}" width="12" height="12" fill="#fff5f5" '
    f'stroke="#c0392b" stroke-width="1" rx="1"/>'
    f'<text x="{PAD_L+346}" y="{leg_y+9}" font-size="9" fill="#c0392b">High-disability MS region</text>'
)

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="20" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    Brainstem and Cerebellum: Most Myelin, Least Neuroprotective Reserve
  </text>
  <text x="{FW//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas ISH · Plp1 (myelin) vs Bdnf (neuroprotection) · sorted by Plp1/Bdnf ratio
  </text>
  <text x="{FW//2}" y="55" text-anchor="middle" font-size="10" fill="#444">
    Cerebellum (41×) and brainstem (29×) have the highest myelin-to-BDNF imbalance — least capacity for neuroprotective repair
  </text>
  {axes}{bars}{legend}
</svg>"""

with open(os.path.join(OUT, "ms_resilience.svg"), "w") as f:
    f.write(svg1)
print("Wrote ms_resilience.svg")


# ── Figure 2: DYSF/CNS resilience severity model ─────────────────────────────
FW2, FH2 = 700, 380

svg2_body = ""

# Two pathways: Immune (susceptibility) vs CNS (severity)
# Layout: top row = susceptibility factors, bottom = severity/resilience

IMM_COL = "#2980b9"
CNS_COL = "#c0392b"
BOX_W = 130; BOX_H = 48

def box(cx, cy, label, col, sublabel=""):
    bx0 = cx - BOX_W//2; by0 = cy - BOX_H//2
    out = (f'<rect x="{bx0}" y="{by0}" width="{BOX_W}" height="{BOX_H}" '
           f'rx="6" fill="{col}" opacity="0.1"/>'
           f'<rect x="{bx0}" y="{by0}" width="{BOX_W}" height="{BOX_H}" '
           f'rx="6" fill="none" stroke="{col}" stroke-width="1.5"/>')
    lines = label.split("\n")
    for j, line in enumerate(lines):
        ly = cy - (len(lines)-1)*7 + j*14
        out += (f'<text x="{cx}" y="{ly}" text-anchor="middle" '
                f'font-size="9.5" font-weight="700" fill="{col}">{line}</text>')
    if sublabel:
        out += (f'<text x="{cx}" y="{cy+BOX_H//2+14}" text-anchor="middle" '
                f'font-size="8" fill="#888">{sublabel}</text>')
    return out

def arrow(x1, y1, x2, y2, col, dash=False):
    d = f'stroke-dasharray="5,3"' if dash else ""
    dx = x2-x1; dy = y2-y1
    L = (dx**2+dy**2)**0.5
    ux = dx/L; uy = dy/L
    ax = x2-ux*8; ay = y2-uy*8
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{ax:.1f}" y2="{ay:.1f}" '
            f'stroke="{col}" stroke-width="1.8" {d}/>'
            f'<polygon points="{x2:.1f},{y2:.1f} {ax-uy*4:.1f},{ay+ux*4:.1f} '
            f'{ax+uy*4:.1f},{ay-ux*4:.1f}" fill="{col}"/>')

ROW1_Y = 90   # susceptibility row
ROW2_Y = 220  # severity row
ROW3_Y = 320  # outcome

# Susceptibility track (top)
svg2_body += box(110, ROW1_Y, "HLA-DRB1*15:01\n+ 200 SNPs", IMM_COL, "immune susceptibility")
svg2_body += arrow(175, ROW1_Y, 265, ROW1_Y, IMM_COL)
svg2_body += box(340, ROW1_Y, "CNS immune\nattack", IMM_COL, "demyelination / lesion")
svg2_body += arrow(405, ROW1_Y, 495, ROW1_Y, IMM_COL)
svg2_body += box(575, ROW1_Y, "Lesion\nlocation", IMM_COL, "determines symptoms")

# Severity/resilience track (bottom)
svg2_body += box(110, ROW2_Y, "DYSF-ZNF638\nlocus", CNS_COL, "rs10191329 (Nature 2023)")
svg2_body += arrow(175, ROW2_Y, 265, ROW2_Y, CNS_COL)
svg2_body += box(340, ROW2_Y, "CNS membrane\nresilience", CNS_COL, "repair capacity after damage")
svg2_body += arrow(405, ROW2_Y, 495, ROW2_Y, CNS_COL)
svg2_body += box(575, ROW2_Y, "Disability\naccumulation", CNS_COL, "−3.7 yr to walking aid\nin risk allele homozygotes")

# Vertical connection between tracks
svg2_body += arrow(340, ROW1_Y+24, 340, ROW2_Y-24, "#7f8c8d", dash=True)
svg2_body += (f'<text x="340" y="{(ROW1_Y+ROW2_Y)//2+4}" text-anchor="middle" '
              f'font-size="8.5" fill="#7f8c8d">regional</text>'
              f'<text x="340" y="{(ROW1_Y+ROW2_Y)//2+16}" text-anchor="middle" '
              f'font-size="8.5" fill="#7f8c8d">expression</text>')

# Regional resilience box
REG_BOX_Y = ROW2_Y + 40
svg2_body += (
    f'<rect x="170" y="{REG_BOX_Y}" width="340" height="40" rx="4" '
    f'fill="#fff8f0" stroke="#e67e22" stroke-width="1.2"/>'
    f'<text x="340" y="{REG_BOX_Y+15}" text-anchor="middle" font-size="9.5" '
    f'fill="#e67e22" font-weight="700">Regional myelin-to-BDNF balance</text>'
    f'<text x="340" y="{REG_BOX_Y+30}" text-anchor="middle" font-size="8.5" fill="#555">'
    f'Cereb (41×) · BS (29×) · Thal (16×) — high myelin load, low neuroprotective reserve</text>'
)

# Key finding box at bottom
KF_Y = FH2 - 55
svg2_body += (
    f'<rect x="30" y="{KF_Y}" width="{FW2-60}" height="42" rx="4" '
    f'fill="#fff5f5" stroke="#c0392b" stroke-width="1.2"/>'
    f'<text x="{FW2//2}" y="{KF_Y+14}" text-anchor="middle" font-size="9.5" '
    f'fill="#c0392b" font-weight="700">IMSGC / Oksenberg lab key finding (Nature 2023):</text>'
    f'<text x="{FW2//2}" y="{KF_Y+28}" text-anchor="middle" font-size="8.5" fill="#444">'
    f'DYSF-ZNF638 is the first MS severity locus — heritability enriched in CNS tissue, not immune tissue.</text>'
    f'<text x="{FW2//2}" y="{KF_Y+40}" text-anchor="middle" font-size="8.5" fill="#444">'
    f'Risk allele (rs10191329) shortens median time to walking aid by 3.7 years in homozygous carriers.</text>'
)

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="700" fill="#222">
    Two Genetic Programs in MS: Susceptibility vs. Severity
  </text>
  <text x="{FW2//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Immune genetics (HLA) determines who gets MS; CNS genetics (DYSF-ZNF638) determines how fast it progresses
  </text>
  <text x="340" y="64" text-anchor="middle" font-size="10" fill="{IMM_COL}" font-weight="700">
    Susceptibility track (immune)
  </text>
  <text x="340" y="194" text-anchor="middle" font-size="10" fill="{CNS_COL}" font-weight="700">
    Severity track (CNS resilience)
  </text>
  {svg2_body}
  <text x="{FW2//2}" y="{FH2-4}" text-anchor="middle" font-size="8.5" fill="#aaa">
    International Multiple Sclerosis Genetics Consortium (2023) Nature 619:323–331 · Oksenberg lab, UCSF
  </text>
</svg>"""

with open(os.path.join(OUT, "ms_severity_model.svg"), "w") as f:
    f.write(svg2)
print("Wrote ms_severity_model.svg")
