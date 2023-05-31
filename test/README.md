# CI tests written with Bats

The following files are part of the SCARAP test suite, written with the [:bat: Bats framework](https://github.com/bats-core/bats-core).
The test suite is run when there are new pushes to the repository, using a [github-actions workflow]().
It collects the input testdata needed for the scripts from [a repository containing the testdata](https://github.com/LebeerLab/SCARAP-testdata).

## Running the test suite locally

Some preparation is needed to run the suite locally however. In short, bats-core needs to be installed and the testdata needs to be cloned into the right folder.

- Installing bats-core: the easiest way to install bats localy is by adding it as a git submodule. [Follow the instructions in the manual](https://bats-core.readthedocs.io/en/stable/tutorial.html#quick-installation) to do this.
- Fetching the testdata: for the scripts to work, the [data folder from the external repo](https://github.com/LebeerLab/SCARAP-testdata) needs to be moved into a new directory in the root called "testdata". The easiest way to do this would be  executing the following code in the root of this project:
  
  ```
  git clone https://github.com/LebeerLab/SCARAP-testdata.git
  mv SCARAP-testdata testdata
  ```
