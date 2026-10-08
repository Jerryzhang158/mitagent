import json
import csv
import tempfile
import unittest
from pathlib import Path

from literature_mining import LiteratureMiner
from local_pubmed_backend import LocalPubMedBackend, normalize_mirna
from pubmed_year_export.build_local_index import (
    article_row,
    fts_row,
    init_database,
)


def record(pmid: str, title: str, abstract: str, doi: str) -> dict:
    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "abstract_sections": [{"label": "RESULTS", "text": abstract}],
        "other_abstracts": [
            {"language": "eng", "type": "plain-language-summary", "text": "Plain summary."}
        ],
        "authors": [],
        "languages": ["eng"],
        "journal": {"title": "Test Journal", "publication_date": {"year": "2025"}},
        "article_dates": [],
        "publication_types": [],
        "mesh_headings": [{"descriptor": "MicroRNAs", "qualifiers": []}],
        "keywords": [{"text": "TXNIP"}],
        "identifiers": {"doi": doi, "pmc": f"PMC{pmid}"},
        "date_revised": {"year": "2025"},
        "source_file": "pubmed26n1000.xml.gz",
    }


def build_test_database(path: Path) -> None:
    records = [
        record(
            "1",
            "hsa-miR-27a-3p directly targets TXNIP",
            "The complete primary abstract demonstrates direct repression of TXNIP.",
            "10.1/exact",
        ),
        record(
            "2",
            "miR-27a regulates TXNIP",
            "A family-level experiment without a mature arm designation.",
            "10.1/core",
        ),
        record(
            "3",
            "miR-27a-5p and TXNIP in stress",
            "The opposite mature arm was measured in a stress model.",
            "10.1/other-arm",
        ),
        record(
            "4",
            "Unrelated miRNA study",
            "This article does not mention the requested gene.",
            "10.1/unrelated",
        ),
    ]
    connection = init_database(path)
    with connection:
        connection.executemany(
            "INSERT INTO articles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [article_row(item) for item in records],
        )
        connection.executemany(
            "INSERT INTO articles_fts VALUES (?,?,?,?,?)",
            [fts_row(item) for item in records],
        )
    connection.close()


class LocalPubMedBackendTest(unittest.TestCase):
    def test_normalization_preserves_mature_arm_and_core_family(self):
        terms = normalize_mirna("hsa-miR-27a-3p")
        self.assertEqual(terms.canonical, "miR-27a-3p")
        self.assertEqual(terms.core_family, "miR-27a")
        self.assertTrue(terms.has_mature_arm)

        unprefixed = normalize_mirna("miR-16-5p")
        self.assertEqual(unprefixed.canonical, "miR-16-5p")
        self.assertEqual(unprefixed.core_family, "miR-16")

        let_terms = normalize_mirna("hsa-let-7a-5p")
        self.assertEqual(let_terms.canonical, "let-7a-5p")
        self.assertEqual(let_terms.core_family, "let-7a")

    def test_exact_mature_search_precedes_deduplicated_family_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "pubmed.sqlite"
            build_test_database(database)
            hits = LocalPubMedBackend(database).search_pair(
                "TXNIP", "hsa-miR-27a-3p", limit=3
            )

        self.assertEqual([hit["pmid"] for hit in hits], ["1", "2", "3"])
        self.assertEqual(hits[0]["match_type"], "exact_mature")
        self.assertEqual(
            [hit["match_type"] for hit in hits[1:]],
            ["core_family", "core_family"],
        )
        self.assertEqual(len({hit["pmid"] for hit in hits}), 3)
        self.assertEqual(hits[0]["doi"], "10.1/exact")
        self.assertIn("Plain summary.", hits[0]["complete_abstract"])

    def test_exact_only_mode_excludes_family_fallback_records(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "pubmed.sqlite"
            build_test_database(database)
            hits = LocalPubMedBackend(database).search_pair(
                "TXNIP",
                "miR-27a-3p",
                limit=10,
                match_mode="exact_only",
            )

        self.assertEqual([hit["pmid"] for hit in hits], ["1"])
        self.assertEqual(hits[0]["match_mode"], "exact_only")

    def test_literature_miner_writes_compatible_txt_and_provenance_jsonl(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "pubmed.sqlite"
            output = root / "articles"
            build_test_database(database)

            miner = LiteratureMiner(
                str(output), backend="local", database_path=str(database)
            )
            results = miner.mine_literature(
                [
                    {
                        "Gene": "TXNIP",
                        "miRNA": "hsa-miR-27a-3p",
                        "Direction": "miRNA→Gene",
                    }
                ],
                max_articles=3,
            )

            txt_path = output / "TXNIP_miR-27a.txt"
            jsonl_path = output / "TXNIP_miR-27a.jsonl"
            txt = txt_path.read_text(encoding="utf-8")
            jsonl = [
                json.loads(line)
                for line in jsonl_path.read_text(encoding="utf-8").splitlines()
            ]
            with open(
                output / "mining_summary.csv", encoding="utf-8", newline=""
            ) as handle:
                summary = list(csv.DictReader(handle))

        self.assertEqual(results[0]["Articles_Found"], 3)
        self.assertEqual(results[0]["Backend"], "local")
        self.assertEqual(results[0]["JSONL_File"], "TXNIP_miR-27a.jsonl")
        self.assertIn("Abstract\nThe complete primary abstract", txt)
        self.assertIn("Plain summary.", txt)
        self.assertEqual(len(jsonl), 3)
        self.assertEqual(jsonl[0]["match_type"], "exact_mature")
        self.assertIn("pmid", jsonl[0])
        self.assertIn("bm25_score", jsonl[0])
        self.assertEqual(summary[0]["Backend"], "local")


if __name__ == "__main__":
    unittest.main()
