# Packaging and source layout

MDscope uses a conventional `src/mdscope/` package layout.
The public imports remain `mdscope.*` and the CLI remains `mdscope`.

Install from a checkout with `python -m pip install -e ".[dev]"`.
For HPC where dependencies already exist, `python -m pip install --no-deps -e .`
registers the checkout in the active environment; do not rely on a global
`PYTHONPATH` pointing to a named parent directory. Run `mdscope --help`,
`python -m mdscope.cli validate <config.yaml>`, or the existing
SLURM scripts after installing.

Build with `python -m build` and validate using `python -m twine check dist/*`.
Check both wheel and source distributions in clean virtual environments from
outside the repository. All seven functional subpackages are included through
setuptools package discovery, not a manually curated list.

The package relocation did not change scientific kernels, YAML schemas, output
paths, units, or trajectory passes. The supported imports are `mdscope.*`;
`mdscope.analysis.*` is not an alias. HPC scripts use the installed lowercase
package consistently.
