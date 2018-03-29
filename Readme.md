# Computational Results To Paper

## Usage

Use [pcn_existance_checker.py](./ff_pcn/pcn_existance_checker.py) to test specific ranges.

For factorizations yafu is used. [yafu.py](./ff_pcn/yafu.py) provides an batchprocessing interface to yafu.


## Installation
1. Clone [yafu-setup-package](https://github.com/KingBowser/yafu-setup-package).
1. Edit Makefile, set `YAFU_OPTS=NFS=1 USE_SSE41=1`.
1. `make all`.
1. Setup paths in `yafu.ini`.

### Note
`yafu.ini` must always placed at current working directory.
