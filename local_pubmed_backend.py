"""Read-only local PubMed FTS5 backend for gene-miRNA pair searches."""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


DEFAULT_DATABASE = (
    Path(__file__).resolve().parent
    / "pubmed_local_backend"
    / "data"
    / "pubmed_mirna_baseline_2026"
    / "pubmed_mirna.sqlite"
)


@dataclass(frozen=True)
class MirnaSearchTerms:
    original: str
    canonical: str
    core_family: str
    body: str
    core_body: str
    is_let_family: bool
    has_mature_arm: bool


def normalize_mirna(mirna: str) -> MirnaSearchTerms:
    """Normalize species-prefixed miRNA names without discarding 3p/5p."""
    original = str(mirna or "").strip()
    if not original:
        raise ValueError("miRNA must not be empty")

    value = original.replace("–", "-").replace("—", "-").strip()
    value = re.sub(r"\.\d+$", "", value)
    value = re.sub(r"(?i)^(?!mir-|let-)[a-z]{3}-", "", value)

    let_match = re.match(r"(?i)^let[-_\s]*(.+)$", value)
    if let_match:
        body = re.sub(r"[-_\s]+", "-", let_match.group(1)).strip("-")
        if not body:
            raise ValueError(f"Cannot normalize miRNA name: {original}")
        core_body = re.sub(r"-(?:3p|5p)$", "", body, flags=re.IGNORECASE)
        return MirnaSearchTerms(
            original=original,
            canonical=f"let-{body}",
            core_family=f"let-{core_body}",
            body=body,
            core_body=core_body,
            is_let_family=True,
            has_mature_arm=core_body != body,
        )

    mir_match = re.match(r"(?i)^(?:micro[-_\s]*rna|mirna|mir)[-_\s]*(.+)$", value)
    if mir_match is None:
        raise ValueError(f"Cannot normalize miRNA name: {original}")
    body = re.sub(r"[-_\s]+", "-", mir_match.group(1)).strip("-")
    if not body or not re.search(r"\d", body):
        raise ValueError(f"Cannot normalize miRNA name: {original}")
    core_body = re.sub(r"-(?:3p|5p)$", "", body, flags=re.IGNORECASE)
    return MirnaSearchTerms(
        original=original,
        canonical=f"miR-{body}",
        core_family=f"miR-{core_body}",
        body=body,
        core_body=core_body,
        is_let_family=False,
        has_mature_arm=core_body != body,
    )


def _tokens(value: str) -> list[str]:
    return re.findall(r"[A-Za-z0-9]+", value)


def _fts_phrase(tokens: Iterable[str]) -> str:
    cleaned = [token for token in tokens if token]
    if not cleaned:
        raise ValueError("FTS phrase must contain at least one alphanumeric token")
    return '"' + " ".join(cleaned) + '"'


def _mirna_alias_phrases(terms: MirnaSearchTerms, exact: bool) -> list[str]:
    body = terms.body if exact else terms.core_body
    body_tokens = _tokens(body)
    if terms.is_let_family:
        aliases = [["let", *body_tokens]]
    else:
        aliases = [
            ["miR", *body_tokens],
            ["miRNA", *body_tokens],
            ["microRNA", *body_tokens],
        ]
    return list(dict.fromkeys(_fts_phrase(alias) for alias in aliases))


def build_fts_query(gene: str, terms: MirnaSearchTerms, exact: bool) -> str:
    gene_phrase = _fts_phrase(_tokens(gene))
    alternatives = " OR ".join(_mirna_alias_phrases(terms, exact))
    return f"{gene_phrase} AND ({alternatives})"


def _json_value(value: str | None, default: Any) -> Any:
    if not value:
        return default
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return default


def _complete_abstract(primary: str, other_abstracts: list[dict[str, Any]]) -> str:
    chunks = [primary.strip()] if primary and primary.strip() else []
    chunks.extend(
        str(item.get("text", "")).strip()
        for item in other_abstracts
        if str(item.get("text", "")).strip()
    )
    return "\n".join(dict.fromkeys(chunks))


class LocalPubMedBackend:
    """Search the reproducible local PubMed corpus without network access."""

    _RANK_EXPRESSION = "bm25(articles_fts, 0.0, 5.0, 2.0, 1.0, 1.0)"

    def __init__(self, database_path: str | Path = DEFAULT_DATABASE):
        self.database_path = Path(database_path).expanduser().resolve()
        if not self.database_path.is_file():
            raise FileNotFoundError(f"Local PubMed database not found: {self.database_path}")
        self._validate_database()

    def _connect(self) -> sqlite3.Connection:
        uri = self.database_path.as_uri() + "?mode=ro"
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only=ON")
        return connection

    def _validate_database(self) -> None:
        connection = self._connect()
        try:
            names = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type IN ('table', 'view')"
                )
            }
        finally:
            connection.close()
        missing = {"articles", "articles_fts"} - names
        if missing:
            raise ValueError(
                f"Invalid local PubMed database; missing tables: {sorted(missing)}"
            )

    def _query(self, fts_query: str, limit: int, match_type: str) -> list[dict[str, Any]]:
        sql = f"""
            SELECT
                articles.pmid,
                articles.title,
                articles.abstract,
                articles.publication_year,
                articles.publication_date_json,
                articles.journal,
                articles.doi,
                articles.pmcid,
                articles.languages_json,
                articles.abstract_sections_json,
                articles.other_abstracts_json,
                articles.authors_json,
                articles.publication_types_json,
                articles.mesh_json,
                articles.keywords_json,
                articles.source_file,
                {self._RANK_EXPRESSION} AS bm25_score
            FROM articles_fts
            JOIN articles ON articles.pmid = articles_fts.pmid
            WHERE articles_fts MATCH ?
            ORDER BY bm25_score, COALESCE(articles.publication_year, 0) DESC
            LIMIT ?
        """
        connection = self._connect()
        try:
            rows = connection.execute(sql, (fts_query, limit)).fetchall()
        finally:
            connection.close()

        hits: list[dict[str, Any]] = []
        for row in rows:
            other_abstracts = _json_value(row["other_abstracts_json"], [])
            score = float(row["bm25_score"])
            hits.append(
                {
                    "pmid": row["pmid"],
                    "title": row["title"],
                    "abstract": row["abstract"],
                    "complete_abstract": _complete_abstract(
                        row["abstract"], other_abstracts
                    ),
                    "abstract_sections": _json_value(
                        row["abstract_sections_json"], []
                    ),
                    "other_abstracts": other_abstracts,
                    "publication_year": row["publication_year"],
                    "publication_date": _json_value(
                        row["publication_date_json"], {}
                    ),
                    "journal": row["journal"],
                    "doi": row["doi"] or "",
                    "pmcid": row["pmcid"] or "",
                    "languages": _json_value(row["languages_json"], []),
                    "authors": _json_value(row["authors_json"], []),
                    "publication_types": _json_value(
                        row["publication_types_json"], []
                    ),
                    "mesh_headings": _json_value(row["mesh_json"], []),
                    "keywords": _json_value(row["keywords_json"], []),
                    "source_file": row["source_file"],
                    "match_type": match_type,
                    "fts_query": fts_query,
                    "bm25_score": score,
                    "relevance_score": -score,
                }
            )
        return hits

    def search_pair(
        self,
        gene: str,
        mirna: str,
        limit: int = 50,
        match_mode: str = "exact_then_family",
    ) -> list[dict[str, Any]]:
        gene = str(gene or "").strip()
        if not gene:
            raise ValueError("gene must not be empty")
        if limit < 1:
            raise ValueError("limit must be >= 1")
        if match_mode not in {"exact_only", "exact_then_family", "core_only"}:
            raise ValueError(
                "match_mode must be exact_only, exact_then_family, or core_only"
            )

        terms = normalize_mirna(mirna)
        results: list[dict[str, Any]] = []
        seen_pmids: set[str] = set()

        if match_mode in {"exact_only", "exact_then_family"}:
            exact_query = build_fts_query(gene, terms, exact=True)
            exact_match_type = (
                "exact_mature" if terms.has_mature_arm else "exact_name"
            )
            for hit in self._query(exact_query, limit, exact_match_type):
                results.append(hit)
                seen_pmids.add(hit["pmid"])

        remaining = limit - len(results)
        if remaining > 0 and match_mode in {"exact_then_family", "core_only"}:
            core_query = build_fts_query(gene, terms, exact=False)
            core_candidates = self._query(
                core_query,
                remaining + len(seen_pmids),
                "core_family",
            )
            for hit in core_candidates:
                if hit["pmid"] in seen_pmids:
                    continue
                results.append(hit)
                seen_pmids.add(hit["pmid"])
                if len(results) >= limit:
                    break

        for rank, hit in enumerate(results, 1):
            hit["rank"] = rank
            hit["query_gene"] = gene
            hit["query_mirna"] = mirna
            hit["normalized_mirna"] = terms.canonical
            hit["core_family"] = terms.core_family
            hit["match_mode"] = match_mode
        return results
