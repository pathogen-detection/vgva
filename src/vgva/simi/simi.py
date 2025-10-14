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
"""

import gc
import logging
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Align import substitution_matrices as aa_substitution_matrices_database
from Bio.AlignIO import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

test_logger = logging.getLogger("pim module")


class PairRecords(object):

    def __init__(self):
        self.query_record_from_PairRecords: SeqRecord = None
        self.ref_record_from_PairRecords: SeqRecord = None

    def read2records(
            self,
            query_ref_path,
            query_identifier: str = "query@",
            ref_identifier: str = "ref@",
            _format: str = "fasta"
    ):
        pair_records = list(SeqIO.parse(query_ref_path,_format))
        assert len(pair_records) == 2, f"class PairRecords 只读取两条序列, 当前序列数量:{len(pair_records)}"
        for p_record in pair_records:
            if query_identifier in p_record.id:
                self.query_record_from_PairRecords = p_record
            elif ref_identifier in p_record.id:
                self.ref_record_from_PairRecords = p_record

    def make2record(
            self,
            query_record: SeqRecord,
            ref_record: SeqRecord,
            query_identifier: str = "query@",
            ref_identifier: str = "ref@",
            know_has_identifier: bool = False
    ):
        if not know_has_identifier:
            if query_identifier not in query_record.id:
                query_record.id = query_identifier + query_record.id
            if ref_identifier not in ref_record.id:
                ref_record.id = ref_identifier + ref_record.id

        self.query_record_from_PairRecords: SeqRecord = query_record
        self.ref_record_from_PairRecords: SeqRecord = ref_record

    def to_disk(self, file_path,format):
        assert self.query_record_from_PairRecords and self.ref_record_from_PairRecords, \
            "query_record_from_PairRecords 和 ref_record_from_PairRecords 的值为None"
        SeqIO.write(
            [self.query_record_from_PairRecords, self.ref_record_from_PairRecords],
            file_path,
            format
        )
        return file_path


class Similarity(object):
    def __init__(self):
        self.translate_codon_table_by_id = None
        self.translate_codon_table_by_name = None
        self.start_codons = None
        self.stop_codons = None
        self.aa_characters = None

    @property
    def substitution_matrices_numVSname(self):
        # example = {1: 'BENNER22', 2: 'BENNER6', 3: 'BENNER74', 4: 'BLASTN', 5: 'BLASTP', 6: 'BLOSUM45', 7: 'BLOSUM50',
        #  8: 'BLOSUM62', 9: 'BLOSUM80', 10: 'BLOSUM90', 11: 'DAYHOFF', 12: 'FENG', 13: 'GENETIC', 14: 'GONNET1992',
        #  15: 'HOXD70', 16: 'JOHNSON', 17: 'JONES', 18: 'LEVIN', 19: 'MCLACHLAN', 20: 'MDM78', 21: 'MEGABLAST',
        #  22: 'NUC.4.4', 23: 'PAM250', 24: 'PAM30', 25: 'PAM70', 26: 'RAO', 27: 'RISLER', 28: 'SCHNEIDER', 29: 'STR',
        #  30: 'TRANS'}
        return {index: sub for index, sub in
                zip(range(1, len(aa_substitution_matrices_database.load()) + 1),
                    aa_substitution_matrices_database.load())}

    def simple_pair_simiarity(
            self,
            input_pair_alignment: PairRecords = None,
            seqtype: str = "nt",
            gap_symbol: str = "-",
            stop_symbol: str = "@",
            to_stop: bool = False,
            codon_table_name: int or str = 1,
            aa_substitution_matrices_name: int or str = "BLOSUM62",
            is_include_gaps: bool = True
    ):
        '''
        pair similarity of two records
        :param input_pair_alignment: alignment of two records
        :param seqtype: nt | nt2aa | protein
        :param gap_symbol:  gap symbol
        :param stop_symbol: stop symbol
        :param to_stop : Translation stops when a stop codon is encountered.
        :param codon_table_name: the name of codon table
        :param aa_substitution_matrices_name: the name of aa substitution matrices
        :param is_include_gaps: Include gaps in similarity calculation?
        :return: similarities of nt, nt2aa or protein
        '''

        assert seqtype in ["nt", "nt2aa", "protein"], "seqtype=nt, protein or nt2aa"
        global translate_codon_table
        global aa_substitution_matrix

        # get alignment of two records
        assert isinstance(input_pair_alignment, PairRecords), "not a PairRecords object"
        query_seq = input_pair_alignment.query_record_from_PairRecords.seq.upper()
        ref_seq = input_pair_alignment.ref_record_from_PairRecords.seq.upper()
        assert len(query_seq) == len(ref_seq), "not an alignment"

        # 1. calculate nt similarities==================================================================================
        if seqtype == "nt":
            identity_base_num = 0
            different_base_num = 0
            for query_base, ref_base in zip(query_seq, ref_seq):
                # do not include gaps
                if (not is_include_gaps) and (query_base == gap_symbol or ref_base == gap_symbol):
                    continue
                # When both bases at the same position are gaps, they are not counted
                if query_base == gap_symbol and ref_base == gap_symbol:
                    continue
                if query_base != ref_base:
                    different_base_num += 1
                    continue
                identity_base_num += 1
            return identity_base_num / (identity_base_num + different_base_num)

        # 1. calculate aa similarities==================================================================================
        # get codon table
        if isinstance(codon_table_name, int):
            translate_codon_table = CodonTable.unambiguous_dna_by_id[codon_table_name]
            self.translate_codon_table_by_id = codon_table_name
        elif isinstance(codon_table_name, str):
            translate_codon_table = CodonTable.unambiguous_dna_by_name[codon_table_name]
            self.translate_codon_table_by_name = codon_table_name

        # get codon table data: start codons, stop codons, and corresponding amino acids.
        self.start_codons = translate_codon_table.start_codons
        self.stop_codons = translate_codon_table.stop_codons
        self.aa_characters = translate_codon_table.forward_table

        # translate nucleotide sequences to amino acid sequences?
        query_protein = query_seq
        ref_protein = ref_seq
        if seqtype == "nt2aa":
            query_protein = query_seq.translate(table=translate_codon_table, stop_symbol=stop_symbol, to_stop=to_stop).upper()
            ref_protein = ref_seq.translate(table=translate_codon_table, stop_symbol=stop_symbol, to_stop=to_stop).upper()

        # remove stop symbols in the terminal protein sequences
        query_protein = query_protein.rstrip(chars=stop_symbol)
        ref_protein = ref_protein.rstrip(chars=stop_symbol)

        # get aa substitution matrices
        if isinstance(aa_substitution_matrices_name, int):
            aa_substitution_matrix = aa_substitution_matrices_database.load(
                self.substitution_matrices_numVSname.get(aa_substitution_matrices_name)
            )
        elif isinstance(aa_substitution_matrices_name, str):
            aa_substitution_matrix = aa_substitution_matrices_database.load(
                aa_substitution_matrices_name
            )

        identity_aa_num: int = 0
        different_aa_num: int = 0
        for query_aa, ref_aa in zip(query_protein, ref_protein):
            # do not include gaps，and query_aa or ref_aa is gaps
            if (not is_include_gaps) and (query_aa == gap_symbol or ref_aa == gap_symbol):
                continue
            # query_aa and ref_aa is gap_symbol or stop_symbol
            # - * - *
            # * - - *

            if (query_aa in [gap_symbol, stop_symbol]) and (ref_aa in [gap_symbol, stop_symbol]):
                continue
            # query_aa or ref_aa is not gap_symbol and stop_symbol
            # M M
            # * -

            if ((query_aa in [gap_symbol, stop_symbol]) and (ref_aa not in [gap_symbol, stop_symbol])) \
                    or ((query_aa not in [gap_symbol, stop_symbol]) and (ref_aa in [gap_symbol, stop_symbol])):
                different_aa_num += 1
                continue

            # Check if amino acid is in library
            upper_aa_pool = list(map(lambda x: x.upper(), self.aa_characters.values()))
            test_logger.debug(f"Query_aa: {query_aa.upper()} | Refer_aa: {ref_aa.upper()} | AA_pool: {upper_aa_pool}")
            if ((query_aa.upper() not in upper_aa_pool) or (ref_aa.upper() not in upper_aa_pool)) \
                    and (query_aa.upper() != ref_aa.upper()):
                test_logger.debug(
                    f"Query_aa:[{query_aa.upper()}] or Refer_aa:[{ref_aa.upper()}] is not in AA_pool"
                    f"and they are different"
                )
                different_aa_num += 1
                continue

            if ((query_aa.upper() not in upper_aa_pool) or (ref_aa.upper() not in upper_aa_pool)) \
                    and (query_aa.upper() == ref_aa.upper()):
                test_logger.debug(
                    f"Query_aa:[{query_aa.upper()}] or Refer_aa:[{ref_aa.upper()}] is not in AA_pool"
                    f"and they are different"
                )
                identity_aa_num += 1
                continue

            # For unequal amino acids, assess their similarity using an amino acid matrix.
            if aa_substitution_matrix[query_aa.upper(), ref_aa.upper()] < 0:
                different_aa_num += 1
                continue
            identity_aa_num += 1

        test_logger.debug(
            f"AA similarity {identity_aa_num}/({identity_aa_num}+{different_aa_num})={identity_aa_num/(identity_aa_num+different_aa_num)}")
        return identity_aa_num / (identity_aa_num + different_aa_num)

    def similarities_collection(
            self,
            alignment: MultipleSeqAlignment,
            seqtype: str = "nt",
            gap_symbol: str = "-",
            stop_symbol: str = "@",
            to_stop: bool = False,
            codon_table_name: int or str = 1,
            aa_substitution_matrices_name: int or str = "BLOSUM62",
            is_include_gaps: bool = True,
            return_type: str = "collection",
    ):
        '''
        This method is exclusively designed for global alignment and computes the pairwise similarity between sequences.
        The input must be a MultipleSeqAlignment object, i.e., the return value of AlignIO.read().
        It calculates the set of pairwise similarity scores among all sequences in the sequence alignment matrix.
        :param alignment: MultipleSeqAlignment object
        :param seqtype: nt | nt2aa | protein
        :param gap_symbol:  gap symbol
        :param stop_symbol: stop symbol
        :param to_stop : Translation stops when a stop codon is encountered.
        :param codon_table_name: the name of codon table
        :param aa_substitution_matrices_name: the name of aa substitution matrices
        :param is_include_gaps: Include gaps in similarity calculation?
        :param return_type: the returned data type， "collection" or "table"
            "collection" a series contained similarities ，"table" a dataframe the header is current_id|next_id|similarities
        :return: similarities(Series||similarities_table["similarities"])
            or table(dataframe, current_id|next_id|similarities, || similarities_table)
        '''
        # calculate the number of records in MultipleSeqAlignment object
        records_num = len(alignment)
        from Bio.SeqRecord import SeqRecord
        current_record: SeqRecord
        similarities_series = []
        for idx, current_record in enumerate(alignment):
            # break, when encounter the last record
            if idx == records_num - 1:
                break
            for next_record in alignment[idx + 1:]:

                pair_similarity_series = pd.Series()
                pair_similarity_series["current_id"] = current_record.id
                pair_similarity_series["next_id"] = next_record.id

                # When both aligned residues are gaps, similarity = np.NaN, continue
                current_seq_unique_base = set(str(current_record.seq))
                next_seq_unique_base = set(str(next_record.seq))
                if len(current_seq_unique_base) == len(next_seq_unique_base) \
                        and len(current_seq_unique_base) == 1 \
                        and list(current_seq_unique_base)[0] == gap_symbol:
                    pair_similarity_series["similarities"] = np.NaN
                    continue

                # calculate the pair similarity
                pair_records = PairRecords()
                pair_records.make2record(
                    query_record=current_record,
                    ref_record=next_record,
                    query_identifier="",
                    ref_identifier="",
                )
                bisimilarity = \
                    self.simple_pair_simiarity(
                        input_pair_alignment=pair_records,
                        seqtype=seqtype,
                        gap_symbol=gap_symbol,
                        stop_symbol=stop_symbol,
                        to_stop=to_stop,
                        codon_table_name=codon_table_name,
                        aa_substitution_matrices_name=aa_substitution_matrices_name,
                        is_include_gaps=is_include_gaps
                    )
                pair_similarity_series["similarities"] = bisimilarity
                similarities_series.append(
                    pair_similarity_series
                )
        similarities_table = pd.concat(similarities_series, axis=1, ignore_index=True).T

        del similarities_series
        gc.collect()

        if return_type == "collection":
            return similarities_table["similarities"]
        elif return_type == "table":
            return similarities_table
        else:
            raise ValueError(
                "return_type must be 'collection' or 'table'"
            )


if __name__ == "__main__":
    # test
    pass


