"""One-off generator for the tiny OT-shaped parquet fixtures under tests/fixtures/opentargets/.

Not run by the test suite itself (the parquet files it writes are committed
alongside the hand-written Reactome text fixtures) — re-run this manually if
the fixture gene universe needs to change. Column layouts mirror the real OT
26.06 `target`, `disease`, and `association_by_datasource_indirect` parquet
schemas exactly (verified against the live files; see ingest.py's module
docstring), just with a handful of rows/genes instead of the real ~80k/3.3GB.

Gene universe: TP53 (clean approved-symbol match), PTGS2/COX2 (approved
symbol PTGS2, deprecated/alias symbol COX2 in symbolSynonyms — the
synonym-mapping case), BRCA1 (a second clean gene to round out interactions/
pathways). FAKEGENE123 deliberately has no corresponding row anywhere in this
target table — it's the unmappable-symbol case.
"""

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

FIXTURES_DIR = Path(__file__).parent / "opentargets"

TARGET_SCHEMA = pa.schema(
    [
        ("id", pa.string()),
        ("approvedSymbol", pa.string()),
        ("biotype", pa.string()),
        ("symbolSynonyms", pa.list_(pa.struct([("label", pa.string()), ("source", pa.string())]))),
        ("obsoleteSymbols", pa.list_(pa.struct([("label", pa.string()), ("source", pa.string())]))),
        ("proteinIds", pa.list_(pa.struct([("id", pa.string()), ("source", pa.string())]))),
        ("transcriptIds", pa.list_(pa.string())),
    ]
)

TARGET_ROWS = [
    {
        "id": "ENSG00000141510",
        "approvedSymbol": "TP53",
        "biotype": "protein_coding",
        "symbolSynonyms": [{"label": "P53", "source": "HGNC"}],
        "obsoleteSymbols": [],
        "proteinIds": [
            {"id": "ENSP00000269305", "source": "ensembl_PRO"},
            {"id": "P04637", "source": "uniprot_swissprot"},
        ],
        "transcriptIds": ["ENST00000269305"],
    },
    {
        "id": "ENSG00000073756",
        "approvedSymbol": "PTGS2",
        "biotype": "protein_coding",
        "symbolSynonyms": [{"label": "COX2", "source": "HGNC"}, {"label": "COX-2", "source": "uniprot"}],
        "obsoleteSymbols": [],
        "proteinIds": [{"id": "ENSP00000263429", "source": "ensembl_PRO"}],
        "transcriptIds": ["ENST00000263429"],
    },
    {
        "id": "ENSG00000012048",
        "approvedSymbol": "BRCA1",
        "biotype": "protein_coding",
        "symbolSynonyms": [{"label": "BRCC1", "source": "HGNC"}],
        "obsoleteSymbols": [],
        "proteinIds": [{"id": "ENSP00000350283", "source": "ensembl_PRO"}],
        "transcriptIds": ["ENST00000357654"],
    },
]

DISEASE_SCHEMA = pa.schema(
    [
        ("id", pa.string()),
        ("name", pa.string()),
        ("dbXRefs", pa.list_(pa.string())),
    ]
)

DISEASE_ROWS = [
    {"id": "EFO_9000001", "name": "hackathon test disease", "dbXRefs": ["MONDO:9000001"]},
]

ASSOCIATION_SCHEMA = pa.schema(
    [
        ("targetId", pa.string()),
        ("diseaseId", pa.string()),
        ("aggregationType", pa.string()),
        ("aggregationValue", pa.string()),
        ("associationScore", pa.float64()),
    ]
)

ASSOCIATION_ROWS = [
    # genetic_association evidence for the test disease (3 genes — deliberately
    # below MIN_GENETIC_ASSOCIATION_GENES=5; tests pass --allow-low-coverage).
    {"targetId": "ENSG00000141510", "diseaseId": "EFO_9000001", "aggregationType": "datasourceId", "aggregationValue": "eva", "associationScore": 0.8},
    {"targetId": "ENSG00000073756", "diseaseId": "EFO_9000001", "aggregationType": "datasourceId", "aggregationValue": "eva", "associationScore": 0.6},
    {"targetId": "ENSG00000012048", "diseaseId": "EFO_9000001", "aggregationType": "datasourceId", "aggregationValue": "gene2phenotype", "associationScore": 0.9},
    # affected_pathway and known_drug evidence, to exercise datatype mapping breadth.
    {"targetId": "ENSG00000141510", "diseaseId": "EFO_9000001", "aggregationType": "datasourceId", "aggregationValue": "reactome", "associationScore": 0.5},
    {"targetId": "ENSG00000012048", "diseaseId": "EFO_9000001", "aggregationType": "datasourceId", "aggregationValue": "clinical_precedence", "associationScore": 0.7},
    # a different disease, to prove predicate-pushdown filtering actually excludes it.
    {"targetId": "ENSG00000141510", "diseaseId": "EFO_OTHERDISEASE", "aggregationType": "datasourceId", "aggregationValue": "eva", "associationScore": 0.99},
]


def main() -> None:
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(TARGET_ROWS, schema=TARGET_SCHEMA), FIXTURES_DIR / "target.parquet")
    pq.write_table(pa.Table.from_pylist(DISEASE_ROWS, schema=DISEASE_SCHEMA), FIXTURES_DIR / "disease.parquet")
    pq.write_table(pa.Table.from_pylist(ASSOCIATION_ROWS, schema=ASSOCIATION_SCHEMA), FIXTURES_DIR / "association.parquet")


if __name__ == "__main__":
    main()
