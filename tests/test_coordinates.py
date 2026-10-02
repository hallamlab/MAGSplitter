import tempfile
import subprocess
import sys
import unittest
from pathlib import Path

import pandas as pd

from magsplitter.coordinates import restore_coordinates
from magsplitter.dataframe_manipulation import undo_orf_removal_test, combine_rxn_contig_map
from magsplitter.ptools_writing import pf_file_creator


class CoordinatesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.annotation = self.root / 'sample.ORF_annotation_table.txt'
        self.source = self.root / 'sample.ptinput.tsv'
        self.table = pd.DataFrame([
            ['a', 'sample-C1', '100', '200', '+', '500'],
            ['b', 'sample-C2', '7', '30', '-', '40'],
        ], columns=['orf_id', 'seqname', 'start', 'end', 'strand', 'contig_length'])
        self.table.to_csv(self.source, sep='\t', index=False)
        self.template = pd.DataFrame([dict(
            ORF_ID='ID\ta', NAME='NAME\ta', STARTBASE='STARTBASE\t100',
            ENDBASE='ENDBASE\t200', FUNCTION='FUNCTION\tenzyme',
            PRODUCT_TYPE='PRODUCT-TYPE\tP', metacyc=['METACYC\tRXN-1'], ec=['EC\t1.2.3.4'])])
        self.expanded = undo_orf_removal_test(self.template, [['ID\ta', 'ID\tb']])
        self.expanded['CONTIG_ID'] = ['sample-C1', 'sample-C2']

    def test_expansion_restores_member_coordinates_and_preserves_annotation(self):
        fixed = restore_coordinates(self.expanded, self.annotation)
        self.assertEqual(fixed.STARTBASE.tolist(), ['STARTBASE\t100', 'STARTBASE\t30'])
        self.assertEqual(fixed.ENDBASE.tolist(), ['ENDBASE\t200', 'ENDBASE\t7'])
        pf_file_creator(fixed, self.root)
        text = (self.root / '0.pf').read_text()
        self.assertIn('ID\tb\nNAME\tb\nSTARTBASE\t30\nENDBASE\t7\n', text)
        self.assertEqual(text.count('METACYC\tRXN-1'), 2)
        self.assertEqual(text.count('EC\t1.2.3.4'), 2)
        self.assertEqual(self.expanded.iloc[1].STARTBASE, 'STARTBASE\t100')

    def test_cli_filters_unsupported_rna_and_validates_before_replacing_outputs(self):
        pf_file_creator(self.template, self.root)
        with (self.root / '0.pf').open('a') as handle:
            handle.write('ID\trna1\nNAME\trna1\nSTARTBASE\t1\nENDBASE\t3\n'
                         'FUNCTION\tRNA\nPRODUCT-TYPE\ttRNA\n//\n')
        self.annotation.write_text('# ORF_ID\tCONTIG_ID\na\tsample-C1\nb\tsample-C2\n')
        (self.root / 'orfs.tsv').write_text('a\tb\n')
        (self.root / 'contigs.tsv').write_text('sample-C1\toriginal1\t500\nsample-C2\toriginal2\t40\n')
        (self.root / 'bins.tsv').write_text('original1\tmag1\noriginal2\tmag2\n')
        command = [sys.executable, '-m', 'magsplitter.main',
                   '-p', str(self.root / '0.pf'), '-c', str(self.annotation),
                   '-r', str(self.root / 'orfs.tsv'), '-i', str(self.root / 'contigs.tsv'),
                   '-m', str(self.root / 'bins.tsv'), '-o', str(self.root)]
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        output = self.root / 'results/mag2/0.pf'
        original = output.read_text()
        self.assertIn('STARTBASE\t30\nENDBASE\t7', original)
        self.source.unlink()
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('Missing coordinate table', run.stderr)
        self.assertEqual(output.read_text(), original)

    def test_explicit_table(self):
        self.assertEqual(len(restore_coordinates(self.expanded, 'custom.tsv', self.source)), 2)

    def test_no_duplicates_to_expand(self):
        pd.testing.assert_frame_equal(undo_orf_removal_test(self.template, []), self.template)
        pd.testing.assert_frame_equal(undo_orf_removal_test(self.template, [['ID\tmissing', 'ID\tb']]), self.template)

    def test_duplicate_expansion_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            undo_orf_removal_test(self.template, [['ID\ta', 'ID\tb', 'ID\tb']])

    def test_missing_table_rejected(self):
        self.source.unlink()
        with self.assertRaisesRegex(ValueError, '--feature_table'):
            restore_coordinates(self.expanded, self.annotation)

    def test_invalid_metadata_rejected(self):
        for column, value, message in [
            ('start', '0', 'outside'), ('end', '41', 'outside'),
            ('start', '31', 'outside'), ('start', '1.5', 'integers'),
            ('strand', '?', 'strand'), ('seqname', 'wrong', 'contig'),
            ('orf_id', 'absent', 'missing coordinates'), ('orf_id', 'a', 'duplicate'),
        ]:
            with self.subTest(column=column, value=value):
                table = self.table.copy()
                table.loc[1, column] = value
                table.to_csv(self.source, sep='\t', index=False)
                with self.assertRaisesRegex(ValueError, message):
                    restore_coordinates(self.expanded, self.annotation)

    def test_missing_or_ambiguous_contig_assignment_rejected(self):
        mapping = pd.DataFrame({'ORF_ID': ['ID\ta'], 'CONTIG_ID': ['sample-C1']})
        with self.assertRaisesRegex(ValueError, 'Missing ORF-to-contig'):
            combine_rxn_contig_map(self.expanded.drop(columns='CONTIG_ID'), mapping)
        with self.assertRaises(pd.errors.MergeError):
            combine_rxn_contig_map(self.template, pd.concat([mapping, mapping]))


if __name__ == '__main__':
    unittest.main()
