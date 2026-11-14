# Materials Project binary-lithium-alloy screening
import csv
import os
from mp_api.client import MPRester

API_KEY = os.environ.get("MP_API_KEY")
if not API_KEY:
    raise RuntimeError("Set MP_API_KEY in the environment; do not hard-code API keys.")

exclude_elements = {
    "Au", "Ag", "Pt", "Pd", "Rh", "Ir", "Ru", "Os", "Re",
    "Hg", "Pb", "Cd", "Tl",
    "U", "Th", "Ra", "Po", "Tc", "Np", "Pu",
    "Na", "K", "Rb", "Cs", "Fr",
    "Be", "Mg", "Ca", "Sr", "Ba",
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I", "At", "Te", "As",
    "He", "Ne", "Ar", "Kr", "Xe", "Rn",
}

with MPRester(API_KEY) as mpr:
    docs = mpr.materials.summary.search(
        elements=["Li"],
        num_elements=2,
        fields=["material_id", "formula_pretty", "structure",
                "energy_above_hull", "band_gap", "elements"],
    )

filtered = []
for doc in docs:
    other_elements = [str(e) for e in doc.elements if str(e) != "Li"]
    if any(el in exclude_elements for el in other_elements):
        continue
    if doc.energy_above_hull >= 0.05:
        continue
    filtered.append(doc)

os.makedirs("cif_files", exist_ok=True)
os.makedirs("summary", exist_ok=True)
with open("summary/alloy_list.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["mp_id", "formula", "ehull", "band_gap"])
    for doc in filtered:
        cif_path = f"cif_files/{doc.material_id}_{doc.formula_pretty}.cif"
        doc.structure.to(filename=cif_path)
        writer.writerow([doc.material_id, doc.formula_pretty,
                         doc.energy_above_hull, doc.band_gap])

print(f"Completed: {len(filtered)} candidates")
