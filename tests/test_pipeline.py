import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import main_pipeline
from config import PipelineConfig
from examples.create_demo_corpus import create_demo_database
from local_pubmed_backend import LocalPubMedBackend


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = patch.object(PipelineConfig, 'CACHE_DIR', str(self.root / 'cache'))
        self.cache.start()
        self.addCleanup(self.cache.stop)
        self.pairs = self.root / 'pairs.csv'
        self.pairs.write_text('Gene,miRNA\nDEMO1,hsa-miR-999-3p\n', encoding='utf-8')

    def cli(self, *extra):
        argv = ['main_pipeline.py', '--mode', 'direct', '--pairs', str(self.pairs),
                '--output', str(self.root / 'runs'), *extra]
        with patch('sys.argv', argv), contextlib.redirect_stdout(io.StringIO()):
            return main_pipeline.main()

    def report(self):
        return next((self.root / 'runs').glob('run_*/pipeline_report.txt')).read_text(encoding='utf-8')

    def test_candidate_only_run_assigns_no_scores(self):
        self.assertEqual(self.cli('--skip-literature', '--skip-bert'), 0)
        report = self.report()
        self.assertIn('Status: completed', report)
        self.assertIn('skipped; no scores assigned', report)
        self.assertNotIn('Pairs with Overall_Score', report)
        self.assertFalse(list((self.root / 'runs').rglob('mti_validation_results.csv')))

    def test_failed_literature_retrieval_stops_before_scoring(self):
        with patch.object(main_pipeline.MiRNAPipeline, 'step2_literature_mining', return_value=False), \
             patch.object(main_pipeline.MiRNAPipeline, 'step3_bert_validation') as scoring:
            self.assertEqual(self.cli(), 1)
            scoring.assert_not_called()
        self.assertIn('Status: failed', self.report())
        self.assertIn('Literature retrieval: failed', self.report())

    def test_failed_bert_stops_before_network(self):
        with patch.object(main_pipeline.MiRNAPipeline, 'step2_literature_mining', return_value=True), \
             patch.object(main_pipeline.MiRNAPipeline, 'step3_bert_validation', return_value=False), \
             patch.object(main_pipeline.MiRNAPipeline, 'step5_network_generation') as network:
            self.assertEqual(self.cli('--generate-network', '--allow-placeholder-expression'), 1)
            network.assert_not_called()
        self.assertIn('BERT scoring: failed', self.report())

    def test_invalid_step_combinations_are_rejected(self):
        for extra in [('--skip-literature',), ('--skip-bert', '--generate-network'),
                      ('--generate-network',), ('--max-articles', '0')]:
            with self.subTest(extra=extra), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    self.cli(*extra)
                self.assertEqual(caught.exception.code, 2)

    def test_direct_formats_and_duplicate_removal(self):
        with contextlib.redirect_stdout(io.StringIO()):
            pipeline = main_pipeline.MiRNAPipeline(str(self.root / 'formats'))
        data = pd.DataFrame({'gene_name': ['DEMO1', 'DEMO1', None],
                             'miRNA_name': ['hsa-miR-999-3p'] * 3})
        for suffix in ('.csv', '.tsv', '.xlsx'):
            path = self.root / ('pairs' + suffix)
            if suffix == '.xlsx':
                data.to_excel(path, index=False)
            else:
                data.to_csv(path, sep='\t' if suffix == '.tsv' else ',', index=False)
            with self.subTest(suffix=suffix), contextlib.redirect_stdout(io.StringIO()):
                results = pipeline._load_direct_pairs(str(path))
                self.assertEqual(len(results), 1)
                self.assertEqual(results[0]['Gene'], 'DEMO1')
                self.assertEqual(results[0]['Databases'], 0)
                self.assertNotIn('Overall_Score', results[0])

    def test_empty_or_invalid_pairs_fail(self):
        for content in ('Gene,miRNA\n', 'wrong,columns\nx,y\n'):
            self.pairs.write_text(content, encoding='utf-8')
            self.assertEqual(self.cli('--skip-literature', '--skip-bert'), 1)

    def test_cartesian_product_uses_unique_nonempty_entries(self):
        genes = self.root / 'genes.txt'
        mirnas = self.root / 'mirnas.txt'
        genes.write_text('DEMO1\nDEMO2\nDEMO1\n\n', encoding='utf-8')
        mirnas.write_text('hsa-miR-999-3p\nhsa-miR-998-3p\n', encoding='utf-8')
        with contextlib.redirect_stdout(io.StringIO()):
            pipeline = main_pipeline.MiRNAPipeline(str(self.root / 'cartesian'))
            pairs = pipeline._build_direct_pairs(str(genes), str(mirnas))
        self.assertEqual(len(pairs), 4)

    def test_direct_uses_local_retrieval_even_if_online_is_configured(self):
        database = self.root / 'synthetic.sqlite'
        create_demo_database(database)
        with patch.object(PipelineConfig, 'PUBMED_LOCAL_DB', str(database)), \
             patch.object(PipelineConfig, 'PUBMED_BACKEND', 'eutils'), \
             patch('requests.get', side_effect=AssertionError('network access is forbidden')):
            self.assertEqual(self.cli('--skip-bert'), 0)
        summary = pd.read_csv(next((self.root / 'runs').rglob('mining_summary.csv')))
        self.assertEqual(summary.loc[0, 'Backend'], 'local')
        self.assertEqual(summary.loc[0, 'Articles_Found'], 1)

    def test_model_arguments_are_recorded_and_forwarded(self):
        observed = []
        def scoring(instance, functions, threshold):
            observed.append((instance.bert_model, instance.llm_model, threshold))
            return True
        with patch.object(main_pipeline.MiRNAPipeline, 'step2_literature_mining', return_value=True), \
             patch.object(main_pipeline.MiRNAPipeline, 'step3_bert_validation', scoring):
            self.assertEqual(self.cli('--bert-model', 'local-model', '--llm-model', 'test-model'), 0)
        self.assertEqual(observed, [('local-model', 'test-model', 30)])
        parameters = json.loads(next((self.root / 'runs').rglob('run_parameters.json')).read_text())
        self.assertEqual(parameters['arguments']['bert_model'], 'local-model')

    def test_ambiguous_mature_arms_fail_before_overwriting_evidence(self):
        database = self.root / 'synthetic.sqlite'
        create_demo_database(database)
        self.pairs.write_text('Gene,miRNA\nDEMO1,hsa-miR-999-3p\nDEMO1,hsa-miR-999-5p\n', encoding='utf-8')
        with patch.object(PipelineConfig, 'PUBMED_LOCAL_DB', str(database)):
            self.assertEqual(self.cli('--skip-bert'), 1)
        self.assertIn('Literature retrieval: failed', self.report())
        self.assertFalse(list((self.root / 'runs').rglob('DEMO1_miR-999.txt')))

    def test_backend_query_error_is_not_reported_as_absent_literature(self):
        database = self.root / 'synthetic.sqlite'
        create_demo_database(database)
        with patch.object(PipelineConfig, 'PUBMED_LOCAL_DB', str(database)), \
             patch.object(LocalPubMedBackend, 'search_pair', side_effect=RuntimeError('fixture query failure')):
            self.assertEqual(self.cli('--skip-bert'), 1)
        self.assertIn('Literature retrieval: failed', self.report())
        summary = pd.read_csv(next((self.root / 'runs').rglob('mining_summary.csv')))
        self.assertIn('fixture query failure', summary.loc[0, 'Error'])

    def test_bert_subprocess_receives_selected_model(self):
        with contextlib.redirect_stdout(io.StringIO()):
            pipeline = main_pipeline.MiRNAPipeline(str(self.root / 'model'), bert_model='local-model')
        pipeline.mining_results = [{'Gene': 'DEMO1', 'miRNA': 'hsa-miR-999-3p',
            'miRNA_Search_Term': 'miR-999', 'Articles_Found': 1, 'File': 'DEMO1_miR-999.txt',
            'Direction': 'DirectInput', 'Databases': 0, 'Database_Sources': 'DirectInput'}]
        def inference(command, check):
            self.assertEqual(command[command.index('--model') + 1], 'local-model')
            pd.DataFrame([{'miRNA': 'miR-999', 'Target_Gene': 'DEMO1', 'Overall_Score': 31}]).to_csv(
                command[command.index('-o') + 1], index=False)
        with patch('main_pipeline.subprocess.run', side_effect=inference), contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(pipeline.step3_bert_validation(['stress'], 30))
        self.assertEqual(pipeline.validation_results.loc[0, 'miRNA'], 'hsa-miR-999-3p')

    def test_nested_llm_error_is_not_success(self):
        with contextlib.redirect_stdout(io.StringIO()):
            pipeline = main_pipeline.MiRNAPipeline(str(self.root / 'llm'))
        pipeline.validation_results = pd.DataFrame([{'miRNA': 'hsa-miR-999-3p', 'Target_Gene': 'DEMO1'}])
        with patch('main_pipeline.MTILLMSummarizer') as summarizer, \
             patch('main_pipeline.check_ollama_connection', return_value=(True, ['test-model'])), \
             contextlib.redirect_stdout(io.StringIO()):
            summarizer.return_value.llm = object()
            summarizer.return_value.batch_summarize.return_value = {'pair': {'standard_summary': {'error': 'fixture'}}}
            self.assertFalse(pipeline.step4_llm_analysis(['stress'], 'standard'))

    def test_synthetic_database_refuses_overwrite(self):
        database = self.root / 'synthetic.sqlite'
        create_demo_database(database)
        self.assertEqual(len(LocalPubMedBackend(database).search_pair('DEMO1', 'hsa-miR-999-3p')), 1)
        with self.assertRaises(FileExistsError):
            create_demo_database(database)

    def test_extractor_configuration_is_separate_and_environment_based(self):
        path = Path(__file__).resolve().parents[1] / 'extractor' / 'config.py'
        spec = importlib.util.spec_from_file_location('extractor_settings', path)
        module = importlib.util.module_from_spec(spec)
        with patch.dict('os.environ', {'NCBI_EMAIL': 'researcher@example.org', 'NCBI_API_KEY': ''}):
            spec.loader.exec_module(module)
        self.assertEqual(module.ENTREZ_EMAIL, 'researcher@example.org')
        self.assertTrue(hasattr(PipelineConfig, 'PUBMED_LOCAL_DB'))


if __name__ == '__main__':
    unittest.main()
