import gzip
import json
import sqlite3
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

from pubmed_year_export.build_local_index import build
from pubmed_year_export.download_pubmed_segments import month_segments


class LocalIndexTest(unittest.TestCase):
    def test_month_segments_cover_requested_year(self):
        segments = list(month_segments(2024, 2024))
        self.assertEqual(len(segments), 12)
        self.assertEqual(segments[1], ("2024-02", "2024/02/01:2024/02/29[dp]"))

    def test_builds_json_manifest_and_searchable_sqlite(self):
        record = {
            "pmid": "123",
            "title": "miR-27a directly targets TXNIP",
            "abstract": "A luciferase assay showed that miR-27a targets TXNIP.",
            "abstract_sections": [{"text": "A luciferase assay showed that miR-27a targets TXNIP."}],
            "other_abstracts": [],
            "authors": [],
            "languages": ["eng"],
            "journal": {"title": "Test", "publication_date": {"year": "2025"}},
            "publication_types": [],
            "mesh_headings": [{"descriptor": "MicroRNAs", "qualifiers": []}],
            "keywords": [{"text": "TXNIP"}],
            "identifiers": {"doi": "10.1/test"},
            "date_revised": {"year": "2025"},
            "source_file": "batch.xml.gz",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "parsed").mkdir()
            (root / "raw").mkdir()
            (root / "state.json").write_text(
                json.dumps(
                    {
                        "query": "test query",
                        "next_retstart": 1,
                        "target_records": 1,
                        "pubmed_match_count": 1,
                    }
                ),
                encoding="utf-8",
            )
            with gzip.open(root / "parsed" / "batch_000.jsonl.gz", "wt", encoding="utf-8") as handle:
                handle.write(json.dumps(record) + "\n")

            build(Namespace(export_dir=root, corpus_version="test-v1"))

            connection = sqlite3.connect(root / "pubmed_mirna.sqlite")
            hits = connection.execute(
                "SELECT pmid FROM articles_fts WHERE articles_fts MATCH ?",
                ('TXNIP AND "miR"',),
            ).fetchall()
            connection.close()
            self.assertEqual(hits, [("123",)])
            self.assertTrue((root / "articles.jsonl.gz").exists())
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["record_count"], 1)

    def test_accepts_completed_baseline_state_and_records_provenance(self):
        record = {
            "pmid": "456",
            "title": "let-7 regulates a target",
            "abstract": "MicroRNA evidence.",
            "journal": {
                "title": "Test",
                "publication_date": {"medline_date": "2024 Winter"},
            },
            "source_file": "pubmed26n0001.xml.gz",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "parsed").mkdir()
            with gzip.open(
                root / "parsed" / "pubmed26n0001.jsonl.gz",
                "wt",
                encoding="utf-8",
            ) as handle:
                handle.write(json.dumps(record) + "\n")
            (root / "baseline_state.json").write_text(
                json.dumps(
                    {
                        "complete": True,
                        "baseline_year": 2026,
                        "base_url": "https://example.test/baseline/",
                        "filter_version": "test-filter-v1",
                        "filter_description": "test filter",
                        "start_year": 2020,
                        "end_year": 2026,
                        "records_scanned": 10,
                        "files": {
                            "pubmed26n0001.xml.gz": {
                                "completed": True,
                                "official_md5": "abc",
                                "records_scanned": 10,
                                "records_selected": 1,
                                "filtered_sha256": "def",
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )

            build(Namespace(export_dir=root, corpus_version="baseline-test-v1"))

            connection = sqlite3.connect(root / "pubmed_mirna.sqlite")
            year = connection.execute(
                "SELECT publication_year FROM articles WHERE pmid = '456'"
            ).fetchone()[0]
            connection.close()
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))

        self.assertEqual(year, 2024)
        self.assertEqual(manifest["baseline"]["year"], 2026)
        self.assertEqual(
            manifest["baseline"]["source_files"][0]["official_md5"], "abc"
        )


if __name__ == "__main__":
    unittest.main()
