
This is SYOTools, available from PyPI:

> pip install syotools

This is SYOTools, available from TestPyPI:

> python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps syotools==1.0.6

this reconfiguration of the SYOTools repo is set up to be packaged so that the end user
can "pip install syotools" instead of cloning the repo. This should make things easier
for the google colab platform

## Development Setup

You need `PYSYN_CDBS` and `SCI_ENG_DIR` set in your environment. I use `venv` and add this to the activate script:

```bash
export PYSYN_CDBS=$(python -c "import sys; print([p for p in sys.path if 'site-packages' in p][0])")/syotools/reference_data/pysynphot_data/
export SCI_ENG_DIR=$(python -c "import sys; print([p for p in sys.path if 'site-packages' in p][0])")/hwo_sci_eng
```