ChemicalComposition
====================


Originally, the ChemicalComposition was part of different Python packages from the fufezan-lab, but we decided to create its own package for the sake of modularity.

Contributors
------------

- J.Leufken
- S.Schulze
- M.Koesters
- C.Fufezan


Installation from  source
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download chemical-composition using `GitHub`_ **or** the zip file:

    GitHub version: Starting from your command line, the easiest way is to clone the GitHub repo.::

    user@localhost:~$ git clone git@github.com:computational-ms/chemical-composition.git

2. Next, navigate into the Ursgal folder and install the requirements::

    user@localhost:~$ cd chemical-composition
    user@localhost:~/ursgal$ pip install -r requirements.txt

3. Finally, use setup.py to install chemical-composition into the Python site-packages::

    user@localhost:~/ursgal$ python setup.py install

(You might need administrator privileges to write in the Python site-package folder.
On Linux or OS X, use ```sudo python setup.py install``` or write into a user folder
by using this command ```python setup.py install --user```. On Windows, you have to
start the command line with administrator privileges.)
