## Plan: Add `retrieve structure` with format standardization

Implement a new `protein-quest retrieve structure` subcommand that reads `search structure` CSV output, downloads each `model_url`, and writes provider-prefixed files. Add a `--format` flag to standardize output format to one of `StructureFileExtensions`, with default `.cif.gz`.

### Steps
1. Add CLI parser `retrieve structure` in `src/protein_quest/cli.py`.
2. Positional args: `structures_csv`, `output_dir`.
3. Optional args: `--format` (choices from `StructureFileExtensions`, default `.cif.gz`), `--max-parallel-downloads`, and existing cacher args.
4. Register parser in `_add_retrieve_subcommands()`.
5. Add handler `_handle_retrieve_structure()` with provenance decorator `@prov(input_files=["structures_csv"], output_dirs=["output_dir"])`.
6. Read rows from CSV (`provider`, `model_identifier`, `model_url`, `model_format`) and validate required fields.
7. Build output filenames as `{provider}_{model_identifier}{native_ext}` to avoid collisions.
8. Infer native extension from `model_format`: `PDB -> .pdb`, `MMCIF -> .cif`, `BCIF -> .bcif`; preserve `.gz` if URL suffix is gzipped.
9. Download with existing async helper (`utils.retrieve_files`) and configured cacher.
10. Post-process each downloaded file: if extension differs from requested `--format`, convert to target extension using IO helpers (`read_structure`, `write_structure`, `gunzip_file`, `split_name_and_extension`).
11. Skip conversion when already in requested format.
12. Report summary counts (downloaded, converted, final files written).
13. Register command dispatch in `HANDLERS` as `("retrieve", "structure")`.

### Tests
1. Add unit tests in `tests/pdbe_3dbeacons/test_fetch.py` for:
   - format-to-extension mapping
   - provider-prefixed filename generation
   - gz inference behavior
2. Add CLI tests in `tests/test_cli.py` for:
   - default `--format .cif.gz`
   - custom `--format` (example: `.pdb`)
   - naming collision prevention via provider prefix
3. Run regression tests for IO conversion (`tests/test_io.py`).

### Documentation
1. Update `README.md` with new command:
   - `protein-quest retrieve structure structures.csv downloads-structure/`
   - default standardization to `.cif.gz`
   - custom example with `--format`.

### Verification
1. `pytest tests/pdbe_3dbeacons/test_fetch.py -k "retrieve or format or naming"`
2. `pytest tests/test_cli.py -k "retrieve_structure or search_structure"`
3. `pytest tests/test_io.py`
4. Manual smoke test with generated `structures.csv`.

### Confirmed decisions
1. Use `--format` for standardization.
2. Default `--format` is `.cif.gz`.
3. Use provider-prefixed filenames (`{provider}_{model_identifier}{ext}`).

TODO remove once implemented
