# AlphaFold Server custom-template preparation

Standalone Python CLI, requiring Python 3.10+ and Biopython (`pip install biopython`).

```
python af_server_templates.py config.json --out prepared_templates
```

Provide one query sequence and 1–4 template specifications. Paths resolve relative to the configuration file. The query can be an amino-acid string or `{"fasta":"query.fasta"}` (one FASTA record).

## Automatic mode

```json
{
  "sequence": "ACDEFGHIKLMNPQRSTVWY",
  "templates": [
    {"file": "structure.cif", "chain": "A", "mode": "auto", "model": 1}
  ]
}
```

`chain` is the mmCIF **label_asym_id**, which may differ from the displayed author chain ID. Omit it only when exactly one coordinate-containing polymer chain is present. One selected chain is exported per template; other chains and nonpolymer molecules are excluded. Model numbers refer to `_atom_site.pdbx_PDB_model_num` (usually 1).

Automatic mode finds a maximum-score semiglobal alignment using BLOSUM62, gap-open −10, gap-extension −0.5, and free terminal gaps. It maps aligned residue pairs, including substitutions and unresolved residues represented in the full mmCIF polymer sequence. It does not invent coordinates. This maximizes alignment score, not raw pair count: maximizing pair count alone would force unrelated residues together. Tied optimal alignments are reported; the first is selected deterministically by Biopython. Review coverage, identity, and residue pairs before upload, particularly for low similarity, repeats, or reordered domains. This initial version does not automatically select a chain or combine manual anchors with automatic alignment.

## Manual mode

All manual inputs use **1-based** query positions. Default template numbering is 1-based mmCIF label_seq_id.

```json
{
  "sequence": "ACDEFGHIKLMNPQRSTVWY",
  "templates": [{
    "file": "structure.cif", "chain": "A", "mode": "manual",
    "numbering": "label",
    "blocks": [
      {"query": [1, 9], "template": [5, 13]},
      {"query": [12, 20], "template": [14, 22]}
    ]
  }]
}
```

Ranges are inclusive, must have equal lengths, and imply a direct one-to-one correspondence. Split them around insertions/deletions. Alternatively, provide `"pairs": [[1,5],[2,6],[12,14]]`; pairs and blocks can be combined if they do not overlap. Use `"numbering":"auth"` with pairs such as `[[1,"42"],[2,"42A"],[3,"43"]]` to address original author residue numbers and insertion codes. Author blocks refer to consecutive polymer positions between two resolved endpoints, not consecutive numerical author labels. For unresolved residues use label numbering. Duplicate, reversed, out-of-range, empty, and entirely unresolved mappings are rejected. Manual mappings may deliberately map substitutions.

## Outputs

- `query_to_hit.json`: Server-compatible list containing only `name`, `queryIndices`, and `templateIndices`.
- `template_0.cif`, etc.: single-chain coordinate files with full polymer sequence tables.
- `mapping_report.json`: source provenance, sequence, coverage, identity, ambiguity and resolved-coordinate status, plus every mapped pair in human-readable 1-based numbering.

JSON output indices are **0-based**. Missing-coordinate residues still occupy indices in the template sequence. Source mmCIF label numbering is retained. Existing output files are protected against overwrite; use a fresh directory. Alternate atom locations are retained.

PDB inputs are supported through Biopython conversion, using the coordinate-derived sequence. PDB SEQRES and missing-residue records are not reconstructed in this version. Missing residue numbers are not treated as evidence of sequence: indices refer to the exported coordinate-only sequence. Prefer mmCIF when unresolved residues matter. Standard and recognized modified amino acids are converted for alignment using Biopython's extended residue dictionary; unrecognized polymer residues become X. Microheterogeneous sequence positions are rejected.

## Server use

In the protein chain's custom-template setting, upload the generated mapping JSON and its corresponding CIF files. For homomers, the same template set can be used for each identical query chain; these are single-chain templates and do not encode a trimeric arrangement. Do not upload mapping_report.json as the mapping file.

## Example validation

`example_tnf.json` reproduces the supplied TNF mappings automatically. Copy the four supplied CIFs into `example_inputs/`, then run:

```
python af_server_templates.py example_tnf.json --out tnf_prepared
```

The original downloaded mapping lists are unsorted; the generated pairs are sorted by query position and are semantically identical. Three examples map 157 residues directly; template 2 maps query indices 0–156 to template indices 1–157 because of an extra N-terminal serine. Tests reproduced every pair in all four supplied mappings and confirmed full-sequence preservation and single-chain output. No live server submission was performed.

Official format references:
- https://alphafoldserver.com/ (Help: custom templates and mapping files)
- https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md (template indexing; local AF3 input schema differs from Server upload schema)
