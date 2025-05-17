# Workflow notes

These notes are here to remind me how to perform basic maintenance
and testing procedures.  They are specific to my computer and
Python/Stata/R/git setup.

Workflow:

- Before posting commits to the devel branch:
  - Run the full set of local tests.
  - Run the full set of local lints.
- Before accepting pull requests from the devel branch
  to the master branch:
  - Run the full set of remote tests.
  - Run the full set of remote lints.
- Before releasing a distribution:
  - Install a test distribution.
  - Run the full set of install tests on the test distribution.
- Immediately after releasing a distribution:
  - Install the released distribution.
  - Run the full set of install tests on the released distribution.

## Testing the Python version

### Local tests

"Local" tests are conducted on my computer on the version of
rcrbounds.py in this folder (or more generally, on the version in
the current working directory).

To perform an interactive test of basic functionality:

1. Use the "Anaconda Prompt" or "Anaconda Powershell Prompt" to open
   a command window that has the correct environment variables to
   run the Anaconda distribution of Python.
2. Execute the shell command:  
   `cd "C:\Users\brian\Desktop\GitHub repos\rcr\python"`  
   (or similar) to move to the folder containing the
   rcrbounds.py version you want to test.
   - The folder should also include the rcr_example.py script.
   - You will also need access to the online data
   (at <http://www.sfu.ca/~bkrauth/code/rcr_example.dta>.
3. Execute the shell command:  
   `ipython`  
   to open an interactive Python session.
4. Execute the ipython command:  
   `%run rcr_example`  
   to run the basic example. You will see two tables of
   results (one OLS and one RCR) for the same regression
   model.
5. Assuming everything has worked out so far, execute
   the ipython command:  
   `%run rcr_example all`  
   to run and display some additional examples including graphs.
6. Execute the ipython command  
   `exit`  
   to exit ipython.

Once basic functionality has been established, you can
perform a full set of unit tests by following these steps:

1. Use the "Anaconda Prompt" or "Anaconda Powershell Prompt" to open
   a command window that has the correct environment variables to
   run the Anaconda distribution of Python.
2. Execute the shell command:  
   `cd "C:\Users\brian\Desktop\GitHub repos\rcr\python"`  
   (or similar) to change to this folder.
   - The folder should also include the testing subfolder with
     all of the tests.
   - You will also need access to the online data
    (at <http://www.sfu.ca/~bkrauth/code/rcr_example.dta>).
3. Execute the shell command:  
   `python -m pytest`  
   to run the full set of unit tests.

### Install tests

"Install" tests are conducted on my computer on an installed
distribution of the RCRBOUNDS pacakge.

1. Use the "Anaconda Prompt" or "Anaconda Powershell Prompt" to open
   a command window that has the correct environment variables to
   run the Anaconda distribution of Python.
2. (optional) Create and activate a virtual environment using
   venv.
3. Execute the shell command:  
   `cd "C:\Users\brian\Desktop\GitHub repos\rcr\python"`  
   (or similar) to change to this folder.
   - The folder should also include the testing subfolder with
     all of the tests.
   - You will also need access to the online data
    (at <http://www.sfu.ca/~bkrauth/code/rcr_example.dta>).
4. Install the version of RCRBOUNDS you wish to test:
   - To install the public distribution posted at PyPi:  
     `pip install rcrbounds`
   - To install the test distribution posted at PyPi:  
     `pip install -i https://test.pypi.org/simple/ rcrbounds`
   - To install the local version of the test distribution in the dist folder:  
     `pip install --no-index --find-links="C:\Users\brian\Desktop\GitHub repos\rcr\python\dist" rcrbounds`  
   - To install the local rcrbounds.py file (including all dependencies)
     and allow the installation to be editable:  
     `pip install --editable .`  
     Note that "editable" means changes to rcrbounds.py
     will be immediately reflected in the installed version.
5. Execute the shell command:  
   `pytest`  
   to run the full set of unit tests.
6. To uninstall the version of RCRBOUNDS you have installed,
   execute the shell command:  
   `python -m pip uninstall rcrbounds`  
   Note that this will not uninstall the dependencies, so
   using a virtual environment would probably be better.

### Remote tests

Remote tests are automatically run by GitHub whenever
a pull request is generated. The tests are run on virtual machines, and include:

1. Running the unit tests on virtual Windows, MacOSX,
   and Linux machines running at least two versions
   of Python.
   - The virtual environments start out "clean" and
     install all required dependencies, so they are
     cleaner tests than the local ones.
2. Running the flake8 and pylint linters. The
   configuration for these linters is slightly more strict
   than for the local versions of the linters, so they
   will sometimes flag things that are not flagged
   locall.

## Linting

### Python linting

flake8 is a basic linter for Python.  To use locally:

1. Use the "Anaconda Prompt" or "Anaconda Powershell Prompt" to open
   a command window that has the correct environment variables to
   run the Anaconda distribution of Python.
2. Execute the shell command:  
   `cd "C:\Users\brian\Desktop\GitHub repos\rcr\python"`  
   (or similar) to change to the folder containing the code
   to lint.
3. Execute the shell command:  
   `flake8`  
   to lint all Python files in the current folder and its
   subfolders.  To lint a particular file, execute:  
   `flake8 rcrbounds.py`  
   for example.

pylint is a more extensive linter. To use locally:

1. Use the "Anaconda Prompt" or "Anaconda Powershell Prompt" to open
   a command window that has the correct environment variables to
   run the Anaconda distribution of Python.
2. Execute the shell command:  
   `cd "C:\Users\brian\Desktop\GitHub repos\rcr\python"`  
   (or similar) to change to the folder containing the code
   to lint.
3. Execute the shell command:  
   `pylint rcrbounds.py rcr_example.py`  
   to lint the Python files in this folder.  For complex reasons,
   pylint needs to be called from the same folder that contains
   the file(s) to be linted, and the filenames have to be provided
   explicitly.

Finally, both flake8 and pylint are run remotely.

## Creating distributions

### Python

To create a distribution for the Python module:

1. Edit the version number in `setup.cfg`.
2. Update build and twine:  
   `python -m pip install --upgrade build`  
   `python -m pip install --upgrade twine`
3. Build the distribution:  
   `python -m build`  
   This will create a distribution (in both `.whl` and
   `.tar.gz` formats) in the dist folder
4. Delete any old distributions in the dist folder.
5. Log into testpypi and go to <https://test.pypi.org/manage/account/#api-tokens>
   to get an API token.
6. Upload the packages to testpypi:  
   `python -m twine upload --repository testpypi dist/*`  
   You will be asked for the API token.
7. If necessary, uninstall the current rcrbounds:  
   `pip uninstall rcrbounds`
8. Wait a minute or two and install the test rcrbounds:  
   `pip install -i https://test.pypi.org/simple/ --no-deps rcrbounds`  
   The `--no-deps` prevents pip from installing any dependencies
   from testpypi since they might be nonstandard.
9. Use pytest to "install test" the package as described above.
10. Log into pypi and get an API token.
11. Upload the packages to pypi:  
    `python -m twine upload dist/*`  
    You will be asked for the API token.
12. Wait a minute or two and install the production rcrbounds:
13. Use pytest to "install test" the package as described above.

Additional details are at <https://packaging.python.org/en/latest/tutorials/packaging-projects/>

### Stata

To update the Stata program:

1. Update the version number and date in setup.cfg.
2. Push the changes to the master branch.
3. Open Stata.
4. Run the Stata command:  
   `ado uninstall rcr`  
   to clear out any current installation of the package.
5. Run the Stata commands:
   `net install rcr, from("https://raw.githubusercontent.com/bvkrauth/rcr/master/stata/")`  
   and (optionally)  
   `net get rcr, from("https://raw.githubusercontent.com/bvkrauth/rcr/master/stata/")`  
6. Run the Stata command `ado describe rcr` to check the version number and date.
7. Run the Stata command:  
   `do rcr_example.do`  
   to run a test case.
8. Transfer the full package of Stata files to the
   SFU web server, etc.
