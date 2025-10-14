#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@filename: genome_variation.py
@author: Bangjun Gong
@contact: bangjungong@foxmail.com
@license: GPL v3
@created: 2025-10-13
@updated: 2025-10-13
@version: 1.0.0
@description:
    Viral genome variation analysis,
    The virus whole-genome alignment matrix was segmented into consecutive sub-matrices (window=100bp, step=25bp)
    by sliding from left to right. For each sub-matrix, the python scripts will calculate the 25th percentile (Q1),
    median, and 75th percentile (Q3), upper bound, lower bound, and outliers of pairwise sequence similarities.
"""

import os
import queue
from io import StringIO
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment
from matplotlib import pyplot as plt
from ..simi.simi import Similarity
from ..stat.statistical_analysis import calculate_boxplot_parameters


class GenomeVariation(Similarity):
    '''
    Viral genome variation analysis,
    The virus whole-genome alignment matrix was segmented into consecutive sub-matrices (window=100bp, step=25bp)
    by sliding from left to right. For each sub-matrix, the python scripts will calculate the 25th percentile (Q1),
    median, and 75th percentile (Q3), upper bound, lower bound, and outliers of pairwise sequence similarities.
    '''
    attribute_buffer = queue.Queue()

    def __init__(self, genome_alignment_input, _format: str = "fasta"):
        '''

        :param genome_alignment_input: str like path, StringIO, MultipleSeqAlignment
        :param _format: the format of sequences file
        '''
        self._out_prefix = os.path.splitext(os.path.basename(genome_alignment_input))[0]
        super(GenomeVariation, self).__init__()
        self.genome_alignment: MultipleSeqAlignment = \
            self.__read_alignment(genome_alignment_input, _format)

        self._genome_alignment_input = genome_alignment_input
        self.subalignments_similarities = None
        self.genome_variation_window = None
        self.genome_variation_step = None

    @staticmethod
    def __read_alignment(seq_input, _format: str = "fasta"):
        if isinstance(seq_input, str) and os.path.isfile(seq_input):
            alignment = AlignIO.read(seq_input, _format)
            return alignment
        elif isinstance(seq_input, StringIO):
            alignment = AlignIO.read(seq_input, _format)
            return alignment
        elif isinstance(seq_input, MultipleSeqAlignment):
            return seq_input
        else:
            raise TypeError(
                "Input must be either a path like string, a StringIO object or a MultipleSeqAlignment object")

    def __callback(self, res):
        GenomeVariation.attribute_buffer.put(res.result())

    def subaln_similarities_attribute(
            self,
            start_pos: int,
            sub_alignment: MultipleSeqAlignment,
            gap_symbol: str,
            stop_symbol: str,
            to_stop: bool,
            codon_table_name: int or str,
            aa_substitution_matrices_name: int or str,
            is_include_gaps: bool,
            outliers_type: str
    ):
        similarities = self.similarities_collection(
            alignment=sub_alignment,
            seqtype="nt",
            gap_symbol=gap_symbol,
            stop_symbol=stop_symbol,
            to_stop=to_stop,
            codon_table_name=codon_table_name,
            aa_substitution_matrices_name=aa_substitution_matrices_name,
            is_include_gaps=is_include_gaps,
            return_type="collection",
        )
        res = calculate_boxplot_parameters(list(similarities), outliers_type)
        outliers = np.NaN if len(res.outliers) <= 0 else "|".join(map(str, res.outliers))
        res_dict = {
            "start_pos": [start_pos],
            "q1": [res.q1],
            "median": [res.median],
            "q3": [res.q3],
            "iqr": [res.iqr],
            "real_lower_bound": [res.real_lower_bound],
            "real_upper_bound": [res.real_upper_bound],
            "outliers": [outliers],
            "min_value": [res.min_value],
            "max_value": [res.max_value],
            "mean": [res.mean],
            "std": [res.std],
        }
        return pd.DataFrame(data=res_dict)

    def get_genome_variation(
            self,
            window: int = 100,
            step: int = 25,
            outliers_type="lower",
            gap_symbol: str = "-",
            is_include_gaps: bool = True,
            threads: int = 30,
    ):
        '''
        The virus whole-genome alignment matrix was segmented into consecutive sub-matrices (window=100bp, step=25bp)
        by sliding from left to right. For each sub-matrix, the python scripts will calculate the 25th percentile (Q1),
        median, and 75th percentile (Q3), upper bound, lower bound, and outliers of pairwise sequence similarities.
        :param window: The columns(length) of sub-matrices.
        :param step: The step size (in bases) by which the window advances when generating each sub-alignment matrix.
        :param outliers_type: The outlier type, both| upper| lower.
        :param gap_symbol: The gap symbol in the alignment matrix.
        :param is_include_gaps: Should gap positions be included when assessing pairwise sequence similarity within each submatrix?
        :param threads: The number of cpu core to use.
        :return: subalignments similarities table
        '''

        self.genome_variation_window = window
        self.genome_variation_step = step

        start_array = np.arange(0, self.genome_alignment.get_alignment_length(), step)
        end_array = start_array + window

        pool = ProcessPoolExecutor(threads)
        for start_pos, end_pos in zip(start_array, end_array):
            sub_alignment = self.genome_alignment[:, start_pos: end_pos]
            pool.submit(
                self.subaln_similarities_attribute,
                start_pos=start_pos,
                sub_alignment=sub_alignment,
                gap_symbol=gap_symbol,
                stop_symbol="@",
                to_stop=False,
                codon_table_name=1,
                aa_substitution_matrices_name="BLOSUM62",
                is_include_gaps=is_include_gaps,
                outliers_type=outliers_type,
            ).add_done_callback(self.__callback)
        pool.shutdown()

        res_dfs = []
        while not GenomeVariation.attribute_buffer.empty():
            res_dfs.append(GenomeVariation.attribute_buffer.get())
        merged_df = pd.concat(res_dfs, ignore_index=True)
        merged_df.sort_values(by=['start_pos'], inplace=True)
        self.subalignments_similarities = merged_df
        return merged_df

    def plot_genome_variation(
            self,
            subaln_simi_table: pd.DataFrame or str=None,
            genome_structure_table: pd.DataFrame=None,
            target_x_start: int = None,
            target_x_end: int = None,
            xaxis_step: int = 200,
            yaxis_limit: tuple = (-0.05, 1.05),
            figures_height: tuple = (0.65, 0.2),
            figures_bottom: tuple = (0.31, 0.02),
            outdir: str = None,
            is_show: bool = True
    ):
        '''
        Viral genome variation map.
        :param subaln_simi_table : Subalignments similarities table
        :param genome_structure_table : Genome structure table
        :param target_x_start : The start position of target genes
        :param target_x_end : The end position of target genes
        :param xaxis_step : Interval between adjacent X-axis ticks
        :param yaxis_limit : The limit of yaxis (min, max)
        :param figures_height : The height of the similarity plot and the genome structure diagram.
        :param figures_bottom : The vertical offset from the bottom of the canvas for the similarity plot and the genome-structure diagram.
        :param outdir : Output Directory
        :param is_show: Display the whole-genome variation analysis plot?
        :return: plt, ax1, (ax2)
        '''
        if subaln_simi_table is not None\
                and isinstance(subaln_simi_table, pd.DataFrame):
            df = subaln_simi_table
        elif subaln_simi_table is not None\
                and isinstance(subaln_simi_table, str) \
                and os.path.isfile(subaln_simi_table):
            df = pd.read_csv(subaln_simi_table)
        elif not subaln_simi_table:
            df = self.subalignments_similarities
        else:
            raise TypeError("subalignments_similarities = None, path like str, dataframe")

        start_pos = df["start_pos"]
        q1 = df["q1"]
        median = df["median"]
        q3 = df["q3"]
        real_min = df["real_lower_bound"]
        real_max = df["real_upper_bound"]

        # plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['axes.linewidth'] = 1.5
        # plt.rcParams['font.family'] = 'Times New Roman'
        fig = plt.figure("similarity attribute", facecolor=None, dpi=100, figsize=(16, 8))
        # 1. plot first figure-similarity===============================================================================
        plt.axes([0.045, figures_bottom[0], 0.92, figures_height[0]])
        ax1 = plt.gca()
        ax1.spines['top'].set_color('none')
        ax1.spines['right'].set_color('none')
        ax1.spines["left"].set_linewidth(1.6)
        ax1.spines["bottom"].set_linewidth(2)

        # set xy limit
        plt.xlim(1, self.genome_alignment.get_alignment_length())
        plt.ylim(*yaxis_limit)

        # set base position (x axis)
        fix_position = np.arange(1, self.genome_alignment.get_alignment_length(), xaxis_step)
        plt.xticks(fix_position, rotation=90, ha="center", size=15)

        # set y axis
        plt.yticks(size=15)
        plt.ylabel("Similarity", size=16)

        # set y grid
        plt.grid(True, which='major', axis='y', linestyle='--', color='gray', linewidth=0.5)

        # plot q1 median q3
        plt.plot(start_pos, q1, label="25th percentile(Q1)",color="#0072B2")  #
        plt.plot(start_pos, median, label="Median Values", color="#FF8C00")  #
        plt.plot(start_pos, q3, label="75th percentile(Q3)", color="#009E73")  #
        plt.fill_between(start_pos, real_min, real_max, alpha=0.3)
        plt.legend()

        # plot outliers
        df_outliers = df[["start_pos", "outliers"]].copy(deep=True)
        df_outliers.dropna(how="any", inplace=True)
        dot_legend = False
        for idx, row in df_outliers.iterrows():
            outliers_start_pos = row["start_pos"]
            outliers = row["outliers"]
            outliers_s = pd.Series(map(lambda x: float(f"{float(x):.2f}"), outliers.split("|"))).value_counts()
            outliers_dict = outliers_s.to_dict()
            for similarity, num in outliers_dict.items():
                plt.scatter(outliers_start_pos, similarity, s=np.log2(num) + 20, label="Outliers",)
                if not dot_legend:
                    plt.legend(bbox_to_anchor=(0.01, 1.05), ncol=4, frameon=True, loc='upper left')
                    dot_legend = True

        target_x = [target_x_start, target_x_start, target_x_end, target_x_end]  # eg [95, 95, 601, 601]
        target_y = [-0.05, 1, 1, -0.05]
        if all(target_x):
            ax1.plot(target_x, target_y, linewidth=3, color="#0000ff")

        # 2. show window & step ========================================================================================
        subax1 = fig.add_axes([0.915, 0.94, 0.05, 0.025])
        # set axes[spines & ticks]
        subax1.spines['top'].set_color('none')
        subax1.spines['right'].set_color('none')
        subax1.spines['bottom'].set_color('none')
        subax1.spines['left'].set_color('none')
        subax1.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        # add window & step parameters
        subax1.text(
            0, 0.5,
            f'Window={self.genome_variation_window} Step={self.genome_variation_step}',
            ha='center',
            va='center',
            fontsize=11,
            # color="#0000ff",
            # fontweight='bold',
        )

        # 3. plot genome structure======================================================================================
        if genome_structure_table is None:
            plt.show()
            return plt, ax1
        # set axes[spines & ticks & limit]
        plt.axes([0.045, figures_bottom[1], 0.92, figures_height[1]])  # plt.axes([0.039, 0.33, 0.94, 0.66])
        ax2 = plt.gca()
        ax2.spines['top'].set_color('none')
        ax2.spines['right'].set_color('none')
        ax2.spines['bottom'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.spines["bottom"].set_position(("data", 0))
        plt.xlim(1, self.genome_alignment.get_alignment_length())
        plt.xticks([])
        plt.yticks([])
        plt.ylim(-1, 1)

        # genome structure
        for idx, row_orf in genome_structure_table.iterrows():
            # row["gene"]:开放阅读框名称, row["start_pos"]开放阅读框起始点,
            # row["end_pos"]开放阅读框终止点, row["orientation"]开放阅读框绘图时在X轴的位置（上/下）
            gene = row_orf["gene"]
            start_pos = row_orf["start_pos"]
            end_pos = row_orf["end_pos"]
            orientation = row_orf["orientation"]
            y_offset = row_orf["y_offset"]
            x_offset = row_orf["x_offset"]

            y = 0 if orientation > 0 else -0.5

            rect = plt.Rectangle(
                (start_pos, y),  # Position in the bottom left corner of the open reading box
                end_pos - start_pos + 1,  # Open reading frame length
                0.5,  # Open reading frame width
                fill=False,
                linewidth=2,
                edgecolor='black'
            )
            plt.gca().add_patch(rect)

            text_y = 1 if orientation > 0 else -1
            plt.text(
                start_pos + x_offset * (end_pos - start_pos),
                y_offset * text_y,
                gene,
                size=12
            )

        if outdir and os.path.isdir(outdir):
            plt.savefig(os.path.join(outdir, f"{self._out_prefix}.svg"))
        if is_show:
            plt.show()
        return plt, ax1, ax2


if __name__ == '__main__':
    pass
