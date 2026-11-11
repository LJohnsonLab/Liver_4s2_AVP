# Liver4s2 bulk RNA-seq

Liver-specific APOE4 → APOE2 switch (AAV-TBG-iCre vs AAV-TBG-GFP), liver and brain bulk RNA-seq.

- **Read the report:** download the repo (Code → Download ZIP), unzip, open `report/_book/index.html`.
- **Analysis:** Quarto book chapters in `report/` (render with `quarto render report/` from the repo root).
- **Processing:** MCC/SLURM scripts in `mcc_pipeline/scripts/`.

Large intermediates (`objects/`, `csv/`, `figures/`) are not tracked; rendering the book regenerates them
from the count matrices.
