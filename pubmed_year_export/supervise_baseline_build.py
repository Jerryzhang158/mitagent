#!/usr/bin/env python3
"""Restart a resumable full Baseline build after transient process failures."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def baseline_complete(output_dir: Path) -> bool:
    state_path = output_dir / "baseline_state.json"
    if not state_path.exists():
        return False
    try:
        return bool(json.loads(state_path.read_text(encoding="utf-8")).get("complete"))
    except (OSError, json.JSONDecodeError):
        return False


def artifacts_complete(output_dir: Path) -> bool:
    return all(
        (output_dir / name).exists()
        for name in (
            "articles.jsonl.gz",
            "pmids.txt.gz",
            "pubmed_mirna.sqlite",
            "manifest.json",
        )
    )


def child_command(args: argparse.Namespace) -> list[str]:
    return [
        sys.executable,
        "-u",
        "-m",
        "pubmed_year_export.build_from_baseline",
        "--output-dir",
        str(args.output_dir),
        "--workers",
        str(args.workers),
        "--timeout",
        str(args.timeout),
        "--retries",
        str(args.retries),
        "--build-local-index",
        "--corpus-version",
        args.corpus_version,
    ]


def run(args: argparse.Namespace) -> int:
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    restarts = 0
    while True:
        if baseline_complete(output_dir) and artifacts_complete(output_dir):
            print(f"[{utc_now()}] Baseline and local index are already complete", flush=True)
            return 0

        print(
            f"[{utc_now()}] Starting resumable build (restart={restarts})",
            flush=True,
        )
        result = subprocess.run(child_command(args), cwd=args.workspace.resolve())
        if result.returncode == 0 and baseline_complete(output_dir) and artifacts_complete(output_dir):
            print(f"[{utc_now()}] Full build completed successfully", flush=True)
            return 0

        restarts += 1
        print(
            f"[{utc_now()}] Child exited with code {result.returncode}; "
            f"retrying in {args.restart_delay:g} seconds",
            file=sys.stderr,
            flush=True,
        )
        if args.max_restarts and restarts >= args.max_restarts:
            print(
                f"[{utc_now()}] Maximum restarts reached ({args.max_restarts})",
                file=sys.stderr,
                flush=True,
            )
            return result.returncode or 1
        time.sleep(args.restart_delay)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("pubmed_local_backend/data/pubmed_mirna_baseline_2026"),
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--retries", type=int, default=10)
    parser.add_argument("--restart-delay", type=float, default=60)
    parser.add_argument(
        "--max-restarts",
        type=int,
        default=100,
        help="maximum child-process restarts; 0 means unlimited",
    )
    parser.add_argument(
        "--corpus-version",
        default="pubmed-mirna-baseline-2026-v1",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.workers < 1 or args.workers > 4:
        print("--workers must be between 1 and 4", file=sys.stderr)
        return 2
    if args.retries < 1 or args.restart_delay < 0 or args.max_restarts < 0:
        print("invalid retry configuration", file=sys.stderr)
        return 2
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
