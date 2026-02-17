# Continuous Integration (CI) tests

The following files are part of the SCARAP test suite, written with the [:bat: Bats framework](https://github.com/bats-core/bats-core).
The test suite is run when there are new pushes to the repository, using a [github-actions workflow]().
It collects the input testdata needed for the scripts from [a repository containing the testdata](https://github.com/LebeerLab/SCARAP-testdata).

## Running the test suite locally

Some preparation is needed to run the suite locally. In short, bats needs to be installed, the test data needs to be downloaded and the test environment needs to be installed. 

First, install bats by adding a few git submodules: 

    mkdir test/test_helper
    git submodule add https://github.com/bats-core/bats-core.git test/bats
    git submodule add https://github.com/bats-core/bats-support.git test/test_helper/bats-support
    git submodule add https://github.com/bats-core/bats-assert.git test/test_helper/bats-assert

(Note: these instructions were taken from [the bats manual](https://bats-core.readthedocs.io/en/stable/tutorial.html#quick-installation).)

Next, download the test data: 

    git clone https://github.com/LebeerLab/SCARAP-testdata.git
    mv SCARAP-testdata testdata

Then, install and activate the test environment: 

    conda env create --file test/environment_max.yml --prefix ./env
    conda activate ./env

Finally, run one of the test scripts as follows (from the project root): 

    ./test/bats/bin/bats test/01_pan.bats
