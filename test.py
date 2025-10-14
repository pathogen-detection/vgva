import os
import pandas as pd
from src.vgva import GenomeVariation
from src.vgva import set_text_editable
set_text_editable(False)

if __name__ == '__main__':
    genome_structure_table = pd.read_excel(R"src/test_data/PRRSV_structure.xlsx", sheet_name="PRRSV-1")
    gv = GenomeVariation(R"src/test_data/EU02_mafft.fasta")
    gv.get_genome_variation(
        window=100,
        step=25,
        outliers_type="lower",
        gap_symbol="-",
        is_include_gaps=True,
        threads=os.cpu_count() * 2
    )
    gv.plot_genome_variation(
        None,
        genome_structure_table=genome_structure_table,
        target_x_start=95,
        target_x_end=650,
        xaxis_step=200,
        yaxis_limit=(-0.05, 1.05),
        figures_height=(0.65, 0.2),
        figures_bottom=(0.31, 0.02),
        outdir=R"src/test_data/",
        is_show=True,
    )
