*****************
Unix Installation
*****************

Installing Java JDK 8+
======================

#.  Check if Java is installed, and if so the version, using::

        $ java -version

#.  If it is not installed, install the latest version of java using the following in a console::

        $ sudo apt install default-jre

#.  After installing the latest version, or if it was already installed, check to insure the version is JDK8+.

    It should output something similar::

        openjdk version "11.0.11" 2021-04-20
        OpenJDK Runtime Environment (build 11.0.11+9-Ubuntu-0ubuntu2)
        OpenJDK 64-Bit Server VM (build 11.0.11+9-Ubuntu-0ubuntu2, mixed mode)

Installing Gradle
=================

#. `Download <https://gradle.org/releases/>`_ the latest gradle distribution
   which comes in two flavors:

     * Binary-only (bin)

     * Complete (all), with documents and sources

   Our recommendation is to download the binary-only version as the documents
   and sources can be viewed online.

#. Unpack the distribution; unzip the distribution zip file in the directory of
   your choosing, e.g.: ::

    $ mkdir /opt/gradle
    $ unzip -d /opt/gradle gradle-7.1-bin.zip
    $ ls /opt/gradle/gradle-7.1
    bin  init.d  lib  LICENSE  NOTICE  README

#. Configure the system environment; configure your ``PATH`` environment
   variable to include the ``bin`` directory of the unzipped distribution by
   editing your bashrc or bash_profile file, e.g.: ::

    $ gedit ~/.bashrc

  and at the bottom line of this bashrc (or bash_profile) file insert the following: ``export PATH=$PATH:/opt/gradle/gradle-7.1/bin``

#. Verify your installation; open a console and run ``gradle -v`` to run gradle and display the version, e.g.: ::

    $ gradle -v
    ------------------------------------------------------------
    Gradle 7.1
    ------------------------------------------------------------

Installing Libraries Dependencies
==================================

#. Installing CFITSIO::

    $ tar xf cfitsio_latest.tar.gz
    $ cd cfitsio-X.XX                   # Replace X.XX with version
    $ ./configure --prefix=/usr/local --enable-sse2 --enable-reentrant
    $ make
    $ make utils
    $ sudo make install
    $ ./testprog > testprog.lis
    $ diff testprog.lis testprog.out    # Should have no output
    $ cmp testprog.fit testprog.std     # Should have no output
    $ rm cookbook fitscopy imcopy smem speed testprog

The CFITSIO library should be tested by building and running
the testprog.c program that is included with the release.
On Unix systems, type::
 
    % make testprog
    % testprog > testprog.lis
    % diff testprog.lis testprog.out
    % cmp testprog.fit testprog.std

On VMS systems, (assuming cc is the name of the C compiler command), type::
 
    $ cc testprog.c
    $ link testprog, cfitsio/lib
    $ run testprog
 
The testprog program should produce a FITS file called `testprog.fit`
that is identical to the testprog.std FITS file included in this
release.  The diagnostic messages (which were piped to the file
testprog.lis in the Unix example) should be identical to the listing
contained in the file testprog.out.  The `diff` and `cmp` commands
shown above should not report any differences in the files.

Try typing `echo $LD_LIBRARY_PATH`

- If it returns nothing, then type

`export LD_LIBRARY_PATH=/usr/local/lib`

- If it returns something, then type

`export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib`

ldconfig

Installing SuperPhot
====================

Make sure the following libraries are installed:
    * cfitsio

Navigate to the path where the SuperPhot repository was downloaded and run the
following commands::

    $ ./gradlew clean
    $ ./gradlew build

********************
Windows Installation
********************

Installing Java JDK 8+
======================

#use MSVC 2017
#TODO test gradle wrapper installation
#TODO edit build.gradle,
#TODO libz still found?
#TODO look up windows ldd equivalent (dependency walker may work check, try cygwin newest), if not try to see what works with the current dll's on Ashkans and hope to god it works :)
#TODO
#TODO See if cfitsio is removable at all, look where its included and if feasible to remove (might be good to do surgery and get rid of it at some point)
#boost requires program options, which arent a header only option
#TODO figure out dll shared library with Ashkan
