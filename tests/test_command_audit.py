import json
from pathlib import Path
import shlex
import tempfile
import unittest
from command_audit import configure_command_log, record_command
from haplotype_split import BACKENDS, split_peaks

class CommandAuditTests(unittest.TestCase):
    def tearDown(self):
        configure_command_log(None)

    def test_exact_arguments_routing_and_append(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'commands.jsonl'
            configure_command_log(path)
            argv = ['cutadapt', '-g', 'ACGT...TGCA;rightmost', '-o', "a path/it's.fastq"]
            for _ in range(2):
                record_command(argv, stdin='pipe:upstream', stdout='captured', stderr='captured')
            records = [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0]['argv'], argv)
            self.assertEqual(shlex.split(records[0]['command']), argv)
            self.assertEqual(records[0]['stdin'], 'pipe:upstream')
            self.assertFalse(records[0]['shell'])

    def test_only_dada2_or_disabled(self):
        self.assertEqual(set(BACKENDS), {'dada2', 'none'})
        for method in ['isonclust', 'amplici', 'gmm2pass', 'gmm2pass+dada2']:
            with self.assertRaises(ValueError):
                split_peaks(method, comps=None, reads_df=None, peak_subsets={}, wt_amplicon_length=336)
