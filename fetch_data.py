"""
fetch_data.py — Oksenberg Lab Demo (Jorge Oksenberg, UCSF MS Genetics)
Downloads Allen Mouse Brain Atlas ISH expression data for myelin and
neuroprotection genes across MS-relevant brain regions.

Key question: which brain regions have the highest myelin gene expression
relative to neuroprotective capacity (BDNF)?

Context: The Oksenberg lab's 2023 Nature paper (IMSGC) identified the
DYSF-ZNF638 locus as the first MS severity locus, pointing to CNS resilience
— not immune factors — as the key determinant of disease progression.
Here we ask whether regional differences in myelin-neuroprotection balance
predict which regions accumulate the most irreversible MS damage.
"""

import urllib.request, json, csv, os

OUT = os.path.dirname(os.path.abspath(__file__))

STRUCTS = {
    "Ctx":   453,   # Cortex — gray matter, increasing recognition of cortical MS
    "BS":    771,   # Brainstem — high lesion impact, major disability driver
    "Cereb": 512,   # Cerebellum — high disability impact, common in progressive MS
    "Thal":  549,   # Thalamus — significant atrophy in MS
    "Hipp":  382,   # Hippocampus — cognitive MS symptoms
    "CP":    247,   # Striatum — relatively spared in MS
    "CC":    484,   # Corpus callosum — primary MS lesion site (white matter)
}

EXPERIMENTS = {
    "Nefl":  73512198,   # Neurofilament light — axonal density / damage marker
    "Plp1":  75496529,   # Proteolipid protein 1 — major CNS myelin gene
    "Mag":   69736418,   # Myelin-associated glycoprotein
    "Gfap":  100098224,  # GFAP — astrocyte reactive potential
    "Bdnf":  77003265,   # BDNF — neuroprotective / remyelination support
}

print("Querying Allen Mouse Brain Atlas ISH API ...")
results = {name: {} for name in STRUCTS}
for gene, exp_id in EXPERIMENTS.items():
    print(f"  {gene} (exp {exp_id}) ...")
    for name, struct_id in STRUCTS.items():
        url = (f"https://api.brain-map.org/api/v2/data/query.json?"
               f"criteria=model::StructureUnionize,"
               f"rma::criteria,section_data_set[id$eq{exp_id}],"
               f"structure[id$eq{struct_id}]&num_rows=1")
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                d = json.loads(r.read())
            msg = d.get("msg", [])
            val = float(msg[0].get("expression_density", 0) or 0) if msg else 0.0
            results[name][gene] = round(val, 6)
        except Exception as e:
            print(f"    WARNING {name}: {e}")
            results[name][gene] = 0.0

with open(os.path.join(OUT, "ms_genes.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["region"] + list(EXPERIMENTS.keys()))
    for name in STRUCTS:
        w.writerow([name] + [results[name][g] for g in EXPERIMENTS])

print("\nWrote ms_genes.tsv")
print("\nMyelin-to-neuroprotection ratio (Plp1 / Bdnf):")
for name in STRUCTS:
    r = results[name]
    ratio = r["Plp1"] / max(r["Bdnf"], 1e-6)
    print(f"  {name:<6}: {ratio:.1f}x  (Plp1={r['Plp1']:.4f}, Bdnf={r['Bdnf']:.4f})")
