ChemicalComposition
=====================

<!-- ![RTD](https://readthedocs.org/projects/unimod-mapper/badge/?version=latest)
![Travis-CI](https://travis-ci.org/computational-ms/unimod-mapper.svg?branch=master) -->

Originally, the chemical_composition was part of different Python packages from the fufezan-lab, but we decided to create its own package for the sake of modularity.


Installation from  source
~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download Ursgal using `GitHub`_ **or** the zip file:

* GitHub version: Starting from your command line, the easiest way is to clone the GitHub repo.::

   user@localhost:~$ git clone https://github.com/ursgal/ursgal.git

* ZIP version: Alternatively, download and extract the `ursgal zip file`_

.. _GitHub:
   https://github.com/ursgal/ursgal

.. _ursgal zip file:
   https://github.com/ursgal/ursgal/archive/master.zip

2. Next, navigate into the Ursgal folder and install the requirements::

    user@localhost:~$ cd ursgal
    user@localhost:~/ursgal$ pip install -r requirements.txt

.. note::

    Pip is included in Python 3.4 and higher. However, it might not be
    included in in your system's PATH environment variable.
    If this is the case, you can either add the Python scripts directory to your
    PATH env variable or use the path to the pip.exe directly for the
    installation, e.g.: ~/Python34/Scripts/pip.exe install -r requirements.txt

.. note::

    On Mac it may be neccesary to use Python3.6, since it comes with its
    own OpenSSL now. This may avoid problems when using pip.

3. Finally, use setup.py to download third-party engines (those that we are allowed to distribute)
and to install Ursgal into the Python site-packages::

    user@localhost:~/ursgal$ python setup.py install

If you want to install the third-party engines without installing Ursgal
into the Python site-packages you can use::

    user@localhost:~/ursgal$ python setup.py install_resources

.. note::

    Since we are not allowed to distribute all third party engines, you might need to
    download and install them on your own. See FAQ (`How to install third party engines`_) and
    the respective engine documentation for more information.

.. note::

    Under Linux, it may be required to change the permission in the
    python site-package folder so that all files are executable

(You might need administrator privileges to write in the Python site-package folder.
On Linux or OS X, use ```sudo python setup.py install``` or write into a user folder
by using this command ```python setup.py install --user```. On Windows, you have to
start the command line with administrator privileges.)



Testing
---------


nosetests

Contributors
------------

* J. Leufken
* S. Schulze
* M. Koesters
* C. Fufezan
