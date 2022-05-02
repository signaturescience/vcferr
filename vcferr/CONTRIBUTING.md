## Local development

Check out a new branch. 
Make changes. 
Build with `pip install .`. 
Test.
Push developmental branch to GitHub.
Open a PR and request review from the maintainers.

## PyPI

FIXME

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