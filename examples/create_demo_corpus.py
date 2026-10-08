"""Create a synthetic FTS5 fixture, not biological observations or PubMed data."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pubmed_year_export.build_local_index import article_row, fts_row, init_database


def create_demo_database(path: Path) -> None:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite existing database: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    article = {
        'pmid': 'synthetic-1',
        'title': 'Synthetic fixture: miR-999-3p and DEMO1',
        'abstract': (
            'This is synthetic test text, not an experimental observation. '
            'In this fictional example, miR-999-3p targets DEMO1 and modulates '
            'the oxidative stress response. A simulated reporter assay supports '
            'binding to the DEMO1 3\'UTR. These statements are for software testing only.'
        ),
        'abstract_sections': [], 'other_abstracts': [], 'authors': [],
        'languages': ['eng'],
        'journal': {'title': 'Synthetic fixture', 'publication_date': {'year': '2025'}},
        'article_dates': [], 'publication_types': [], 'mesh_headings': [],
        'keywords': [{'text': 'miRNA'}], 'identifiers': {},
        'date_revised': {}, 'source_file': 'synthetic-fixture',
    }
    connection = init_database(path)
    try:
        with connection:
            connection.execute('INSERT INTO articles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', article_row(article))
            connection.execute('INSERT INTO articles_fts VALUES (?,?,?,?,?)', fts_row(article))
    finally:
        connection.close()
    path.with_suffix('.manifest.json').write_text(json.dumps({
        'corpus_version': 'mitagent-synthetic-demo-1', 'record_count': 1,
        'synthetic': True, 'description': 'Software fixture; no biological claims.'
    }, indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('outputs/demo/pubmed.sqlite'))
    args = parser.parse_args()
    create_demo_database(args.output)
    print(f'Synthetic database written to {args.output}')


if __name__ == '__main__':
    main()
