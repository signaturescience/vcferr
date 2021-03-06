## Local development

1. Check out a new branch. 
1. Make changes. 
1. Build with `pip install .`. 
1. Test.
1. Push developmental branch to GitHub.
1. Open a PR and request review from the maintainers.

## PyPI

1. Bump version in `setup.py` and `__init__.py`
1. Build the release files by running `pyton -m build` at repo root
1. Optionally remove previous `.whl` and `.tar.gz` versions of package from `dist/`.
1. Upload to PyPI with `twine upload dist/*` (NOTE: This will ask you for username/password. Make sure to set up a token with PyPi and use `__token__` as username and the token value with the pypi- prefix for password.)

## Bioconda

1. Upload new package version to PyPI first. Bioconda's build relies up the updated package at PyPI.
1. Read [documentation on contributing to Bioconda](https://bioconda.github.io/contributor/index.html).
1. Fork bioconda/bioconda-recipes. Add bioconda/bioconda-recipes as upstream and `git pull origin upstream` to sync changes.
1. Create a new branch.
1. Update [vcferr's meta.yaml recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/vcferr/meta.yaml).
   1. Update the version number on line 2. 
   2. Update the sha256 has on line 10. You can get this hash on the PyPI page (e.g., https://pypi.org/project/vcferr/1.0.0/#files), or via `curl -sL https://pypi.io/packages/source/v/vcferr/vcferr-1.0.0.tar.gz | sha256sum` (substitute package version).
1. Locally test updates as described in the [testing recipes locally docs](https://bioconda.github.io/contributor/building-locally.html).
1. Push your changes and open a PR to bioconda/bioconda-recipes.
