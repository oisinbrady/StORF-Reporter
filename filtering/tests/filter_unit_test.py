import sys; sys.path.insert(1, '../')
import os
import argparse
from unittest import mock, TestCase
from storf_filter import filter_by_overlap, filter_by_size_range, filter_by_gc_range, filter_by_stop_codons, get_con_groups, read_fasta, FILTERS
from filter_constants import STORF_CAP_PERCENTAGE, GC_LB, GC_UB, OLAP_LB, SIZE_LB, ARB_MAX_STORF_SIZE, ARB_MAX_OLAP_SIZE


class TestFilters(TestCase):
    @mock.patch.dict(FILTERS, {"min_olap": 0, "max_olap": 50, 'test_disable_o_count': True}, clear=True)
    def test_filter_by_overlap_size_range(self):
        contig = read_fasta("mock_contig.fasta")
        contig_values=[[s[0], 1] for s in contig]
        selection = filter_by_overlap(storf_group_values=contig_values, storf_group=contig)
        expected = [1, 0, 1, 0, 1, 0, 1]  # selected StORFs based of overlap size filter
        self.assertEqual(expected, [v[1] for v in selection])

    @mock.patch.dict(FILTERS, {"min_olap": 0, "max_olap": 50, "max_olap_count": 44, 'test_disable_o_count': False}, clear=True)
    def test_filter_by_overlap_size_and_count(self):
        contig = read_fasta("mock_contig.fasta")
        contig_values=[[s[0], 1] for s in contig]
        selection = filter_by_overlap(storf_group_values=contig_values, storf_group=contig)
        expected = [1, 0, 1, 0, 1, 0, 1]
        self.assertEqual(expected, [v[1] for v in selection])

    @mock.patch.dict(FILTERS, {"min_gc": 0.3, "max_gc": 0.5}, clear=True)
    def test_filter_by_gc_range(self):
        contig = read_fasta("mock_contig.fasta")
        contig_values=[[s[0], 1] for s in contig]
        selection = filter_by_gc_range(storf_group_values=contig_values, storf_group=contig)
        expected = [1, 1, 0, 0, 0, 0, 1]
        self.assertEqual(expected, [v[1] for v in selection])

    @mock.patch.dict(FILTERS, {"stop_codons": ["TGA"]}, clear=True)
    def test_filter_by_stop_codons(self):
        contig = read_fasta("mock_contig.fasta")
        contig_values=[[s[0], 1] for s in contig]
        selection = filter_by_stop_codons(storf_group_values=contig_values, storf_group=contig)
        expected = [1, 0, 0, 0, 0, 0, 0]
        self.assertEqual(expected, [v[1] for v in selection])

    def test_contig_getter(self):
        storfs = read_fasta("mock_storfs.fasta")
        contigs = get_con_groups(storfs)
        expected_contig_number = 6
        expected_ids = [
            ["26306-26441", "26351-26468", "26404-26521", "26513-26636", "26518-26656", "26556-26685", "26739-26859"],
            ["30081-30258", "30104-30236"],
            ["31807-31942"],
            ["35309-35456", "35356-35566", "35384-35495"],
            ["35779-35884"],
            ["45134-45251"]
        ]
        # test the total number of contigs
        self.assertEqual(len(contigs), expected_contig_number)
        # test the exact sequence order for each contig
        for i, contig in enumerate(contigs):
            for j, storf in enumerate(contig):
                id_pretext = ">Chromosome-UR:"
                storf_id = storf[0][len(id_pretext):len(id_pretext)+len(expected_ids[i][j])]
                self.assertEqual(storf_id, expected_ids[i][j])
                

if __name__ == '__main__':
    unittest.main()