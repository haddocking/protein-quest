# wasm version of protein-quest

TODO remove once PR is ready to merge

```shell
uv run marimo convert docs/notebooks/uniprot.ipynb -o docs/notebooks/uniprot.py
# Added inline dep to protein-quest
marimo export html-wasm docs/notebooks/uniprot.py  -o docs/notebookswasm --mode edit
caddy file-server -r docs/notebookswasm -l 0.0.0.0:8080
```

```python
import micropip

gemmi_whl = "http://localhost:8080/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(gemmi_whl)
mmcif_whl = "http://localhost:8080/mmcif-1.1.1-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(mmcif_whl)
pq_whl = "http://localhost:8080/protein_quest-1.6.0-py3-none-any.whl"
await micropip.install(pq_whl)
```


## psutil

Make non-dep in wasm and added workaround

```shell
uv build --wheel
cp dist/protein_quest-1.6.0-py3-none-any.whl docs/notebookswasm/
```

## gemmi

Make wasm using cibuildwheel for python 3.14.2

```shell
gh extension install https://github.com/nektos/gh-act
cp .github/workflows/wheels2.yml .github/workflows/wheels_wasm.yml
# pick medium sized runner image
gh act workflow_dispatch -j build_wasm_wheels --artifact-server-path $PWD/.artifacts
unzip -d ../protein-quest/docs/notebookswasm/ .artifacts/1/cibw-wheels--/cibw-wheels--.zip
```

`.github/workflows/wheels2wasm.yml`:

```yaml
# Use cibuildwheel v2 to build wheels for Python 3.
# Based on https://cibuildwheel.readthedocs.io/en/stable/setup/

name: Wheels2 wasm only

on:
  workflow_dispatch:

jobs:
  build_wasm_wheels:
    name: Wheels wasm
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v5

    - name: Build wheels
      uses: pypa/cibuildwheel@v4.2.0
      env:
        CIBW_BUILD: "cp314-pyodide_wasm32"
        CIBW_PLATFORM: pyodide
        CIBW_ARCHS: wasm32
        # Skip tests, as some others fail for pyodide
        CIBW_TEST_COMMAND_PYODIDE: "true"

    - uses: actions/upload-artifact@v4
      with:
        name: cibw-wheels-wasm
        path: ./wheelhouse/*.whl
```

test
```python
gemmi_whl = "http://localhost:8080/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(gemmi_whl)
import urllib.request
import gemmi

url = "https://www.ebi.ac.uk/pdbe/entry-files/download/3jrs_updated.cif"
structure = gemmi.cif.read_string(urllib.request.urlopen(url).read())
structure
```

## mmcif

```shell
git clone  --recurse-submodules https://github.com/rcsb/py-mmcif
cd py-mmcif
mkdir -p .github/workflows
cp ../gemmi/.github/workflows/wheels2wasm.yml .github/workflows/wheels_wasm.yml
rm  ~/.config/act/actrc 
# use large image with cmake
gh act workflow_dispatch -j build_wasm_wheels --artifact-server-path $PWD/.artifacts
unzip -d ../protein-quest/docs/notebookswasm/ .artifacts/1/cibw-wheels-wasm/cibw-wheels-wasm.zip
```

Needed changes to code see https://github.com/i-VRESSE/py-mmcif fork.

test
```
from mmcif.api.DictionaryApi import DictionaryApi
from mmcif.io.BinaryCifReader import BinaryCifReader
from mmcif.io.BinaryCifWriter import BinaryCifWriter
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter
```

## distributed

Made map_with_progress work when distributed is available or not.
