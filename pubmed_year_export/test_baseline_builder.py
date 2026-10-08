import gzip
import hashlib
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from pubmed_year_export.build_from_baseline import (
    discover_source_files,
    download_with_resume,
    parse_official_md5,
    process_source_file,
    publication_years,
    record_matches,
    select_source_files,
)
from pubmed_year_export.test_parser import SAMPLE
from pubmed_year_export.supervise_baseline_build import (
    artifacts_complete,
    baseline_complete,
    child_command,
)


def matching_record() -> dict:
    return {
        "title": "hsa-miR-27a-3p regulates TXNIP",
        "abstract": "A focused microRNA experiment.",
        "other_abstracts": [],
        "journal": {
            "publication_date": {"medline_date": "2024 Jan-Feb"},
        },
        "article_dates": [],
        "mesh_headings": [],
        "keywords": [],
    }


class BaselineBuilderTest(unittest.TestCase):
    def test_discovers_only_requested_baseline_year(self):
        listing = b"""
        <a href="pubmed26n0002.xml.gz">two</a>
        <a href="pubmed26n0001.xml.gz">one</a>
        <a href="pubmed26n0001.xml.gz.md5">checksum</a>
        <a href="pubmed25n1274.xml.gz">old</a>
        """
        with patch(
            "pubmed_year_export.build_from_baseline.request_bytes",
            return_value=listing,
        ):
            files = discover_source_files("https://example.test/", 2026)
        self.assertEqual(files, ["pubmed26n0001.xml.gz", "pubmed26n0002.xml.gz"])

    def test_selects_a_bounded_range_from_an_exact_source_file(self):
        files = [f"pubmed26n{number:04d}.xml.gz" for number in range(1, 6)]
        self.assertEqual(
            select_source_files(files, "pubmed26n0004.xml.gz", 2),
            ["pubmed26n0004.xml.gz", "pubmed26n0005.xml.gz"],
        )
        with self.assertRaisesRegex(ValueError, "was not discovered"):
            select_source_files(files, "pubmed26n9999.xml.gz", 2)

    def test_parses_common_official_md5_format(self):
        checksum = "MD5 (pubmed26n0001.xml.gz) = 0123456789abcdef0123456789ABCDEF"
        self.assertEqual(
            parse_official_md5(checksum, "pubmed26n0001.xml.gz"),
            "0123456789abcdef0123456789abcdef",
        )

    def test_resumes_partial_download_with_http_range(self):
        content = b"hello world"

        class Response:
            status = 206

            def __init__(self):
                self.chunks = iter([b"world", b""])

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def read(self, _size):
                return next(self.chunks)

        with tempfile.TemporaryDirectory() as directory:
            final = Path(directory) / "source.xml.gz"
            final.with_name(final.name + ".part").write_bytes(b"hello ")

            def open_request(request, timeout):
                self.assertEqual(request.get_header("Range"), "bytes=6-")
                self.assertEqual(timeout, 10)
                return Response()

            with patch(
                "pubmed_year_export.build_from_baseline.urllib.request.urlopen",
                side_effect=open_request,
            ):
                size = download_with_resume(
                    "https://example.test/source.xml.gz",
                    final,
                    hashlib.md5(content).hexdigest(),
                    timeout=10,
                    retries=1,
                )

            self.assertEqual(size, len(content))
            self.assertEqual(final.read_bytes(), content)
            self.assertFalse(final.with_name(final.name + ".part").exists())

    def test_filter_uses_publication_dates_and_precise_mirna_terms(self):
        record = matching_record()
        self.assertEqual(publication_years(record), {2024})
        self.assertEqual(record_matches(record), (True, [2024]))

        record["title"] = "A mirror-image assay for TXNIP"
        record["abstract"] = "No small RNA terminology is present."
        self.assertEqual(record_matches(record), (False, [2024]))

        record["title"] = "miRNA study"
        record["abstract"] = ""
        record["other_abstracts"] = [{"text": "The let-7 family was tested."}]
        self.assertEqual(record_matches(record), (True, [2024]))

    def test_completion_and_revision_dates_do_not_count_as_publication_years(self):
        record = matching_record()
        record["journal"]["publication_date"] = {"year": "2019"}
        record["date_completed"] = {"year": "2025"}
        record["date_revised"] = {"year": "2026"}
        self.assertEqual(publication_years(record), {2019})
        self.assertEqual(record_matches(record), (False, []))

    def test_streams_sample_to_deterministic_filtered_shard(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "pubmed26n0001.xml.gz"
            shard = root / "pubmed26n0001.jsonl.gz"
            with gzip.open(raw, "wb") as handle:
                handle.write(SAMPLE)

            first = process_source_file(raw, shard, 2020, 2026)
            second = process_source_file(raw, shard, 2020, 2026)
            with gzip.open(shard, "rt", encoding="utf-8") as handle:
                records = [json.loads(line) for line in handle if line.strip()]

        self.assertEqual(first[:2], (1, 1))
        self.assertEqual(second[:2], (1, 1))
        self.assertEqual(first[2], second[2])
        self.assertEqual(records[0]["pmid"], "12345678")
        self.assertEqual(
            records[0]["baseline_filter"]["matched_publication_years"], [2025]
        )

    def test_supervisor_completion_checks_and_child_command(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "baseline_state.json").write_text(
                json.dumps({"complete": True}), encoding="utf-8"
            )
            self.assertTrue(baseline_complete(root))
            self.assertFalse(artifacts_complete(root))
            for name in (
                "articles.jsonl.gz",
                "pmids.txt.gz",
                "pubmed_mirna.sqlite",
                "manifest.json",
            ):
                (root / name).touch()
            self.assertTrue(artifacts_complete(root))

        command = child_command(
            Namespace(
                output_dir=Path("output"),
                workers=2,
                timeout=120,
                retries=10,
                corpus_version="test-v1",
            )
        )
        self.assertIn("pubmed_year_export.build_from_baseline", command)
        self.assertIn("--build-local-index", command)


if __name__ == "__main__":
    unittest.main()
