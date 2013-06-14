#!/bin/sh
#
# Customized for PETSc tarball for PyLith by Brad Aagaard.
# Original is $PETSC_DIR/bin/maint/builddist.
#
# This script builds the PETSc tar file distribution
#
# Usage: builddist petscrepo branch
# example usage: builddist /sandbox/petsc/petsc.clone maint
# Notes: version info is automatically obtained from include/petscversion.h
#
#echo "------- Have you done ALL of the following?? -----------------------------"
#echo "(0) Excluded any directories not for public release"
#echo "(1) Set the version number in manual.tex, and intro.tex, website/index.html?"
#echo "(2) Set the version number and release date in petscversion.h?"
#echo "(3) Made sure Fortran include files match C ones?"
#echo "(4) Latest version of win32fe.exe is copied over into /sandbox/petsc/petsc-dist on login.mcs"
#echo "(5) Update update-docs.py with additional docs on the petsc web page (for eg new changes file)"
#echo "(6) Update version info on all petsc web pages"
#echo "(7) tag the new release with 'hg' and make a new clone for the release"

# If version specified on the commandline, set the version
echo "Starting date: `date`"

if [ $# = 2 ]; then
  petscrepo=$1
  branch=$2
elif [ $# = 1 ]; then
  petscrepo=$1
  branch='master'
else
  echo 'Error: petscrepo not specified. Usge: builddist petscrepo'
  exit
fi

# check petscrepo to be valid
if [ ! -d $petscrepo ]; then
  echo 'Error: dir $petscrepo does not exist'
  exit
fi

if [ ! -d $petscrepo/.git ]; then
  echo 'Error: dir $petscrepo is not a git repo'
  exit
fi

# Initialize vars
PETSC_ARCH=arch-build; export PETSC_ARCH
PETSC_DIR=`cd $petscrepo;pwd -P`; export PETSC_DIR

# Clean and Update the git repository and check for branch
cd $PETSC_DIR
#git clean -q -f -d -x
#git fetch -q origin
#git checkout -f origin/$branch
if [ "$?" != "0" ]; then
  echo 'Error: branch: $branch does not exist in $PETSC_DIR'
  exit
fi

pdir=`basename $PETSC_DIR`

# Create a tmp dir
if [ ! -d /tmp ]; then
  echo 'Error: /tmp does not exist!'
  exit
fi
tmpdir=~/tmp/petsc-dist-tmp.$USER.$$
if [ -d $tmpdir ]; then
  /bin/rm -rf $tmpdir
fi
mkdir -p $tmpdir
if [ ! -d $tmpdir ]; then
  echo 'Error: Cannot create $tmpdir'
  exit
fi

# check petscversion.h and set the version string
version_release=`grep '^#define PETSC_VERSION_RELEASE ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`
version_major=`grep '^#define PETSC_VERSION_MAJOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`
version_minor=`grep '^#define PETSC_VERSION_MINOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`
version_subminor=`grep '^#define PETSC_VERSION_SUBMINOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`

#generated a couple of more values
version_date=`date +"%b, %d, %Y"`
version_git=`git log -1 --pretty=format:%H`
version_date_git=`git log -1 --pretty=format:%ci`
if  [ ${version_release} = 0 ]; then
  if [ "$branch" = "master" ]; then
    version=-dev
  elif [ "$branch" = "knepley/pylith" ]; then
    version=-pylith
  else
    version=-`echo $branch | sed -e s#/#-#g`
  fi
elif [ ${version_release} = 1 ]; then
  version=-${version_major}.${version_minor}.${version_subminor}
else
  echo "Unknown PETSC_VERSION_RELEASE: ${version_release}"
  exit
fi
echo "Building ~/petsc$version.tar.gz and ~/petsc-lite$version.tar.gz"

# create fortran stubdocs and cleanup
cd $PETSC_DIR
# we should have an option of running configure - without compilers/mpi/blas etc..
./config/configure.py --with-mpi=0 --with-matlab=0 --with-c2html=0
make allfortranstubs
#make alldoc LOC=${PETSC_DIR}
#
make ACTION=clean tree_basic

# now tar up the PETSc tree
cd $PETSC_DIR/..
cat ${PETSC_DIR}/bin/maint/xclude | sed -e s/petsc-dist/$pdir/ > $tmpdir/xclude
/bin/tar --create --file $tmpdir/petsc.tar --exclude-from  $tmpdir/xclude $pdir
cd $tmpdir
/bin/tar xf $tmpdir/petsc.tar
# just make sure we are not doing 'mv petsc petsc'
if [ ! -d petsc$version ]; then
  /bin/mv $pdir petsc$version
fi
#
# Before Creating the final tarfile, make changes to the bmakefiles/makefiles,
# create tagfiles etc. permissions etc.
#
# Eliminate chmods in the makefile & conf/rules
cd $tmpdir/petsc$version
/bin/mv makefile makefile.bak
/bin/grep -v 'chmod' makefile.bak > makefile
/bin/rm -f makefile.bak

/bin/mv conf/rules conf/rules.bak
/bin/grep -v 'chmod' conf/rules.bak > conf/rules
/bin/rm -f  conf/rules.bak

#add in PETSC_VERSION_DATE, PETSC_VERSION_GIT, PETSC_VERSION_DATE_GIT
echo Using PETSC_VERSION_DATE: ${version_date}
echo Using PETSC_VERSION_GIT: ${version_git}
echo Using PETSC_VERSION_DATE_GIT: ${version_date_git}
/bin/mv include/petscversion.h include/petscversion.h.bak
cat include/petscversion.h.bak | \
  sed -e "s/#define PETSC_VERSION_DATE\ .*/#define PETSC_VERSION_DATE       \"${version_date}\"/" | \
  sed -e "s/#define PETSC_VERSION_GIT\ .*/#define PETSC_VERSION_GIT        \"${version_git}\"/" | \
  sed -e "s/#define PETSC_VERSION_DATE_GIT\ .*/#define PETSC_VERSION_DATE_GIT   \"${version_date_git}\"/" \
  > include/petscversion.h
/bin/rm -f include/petscversion.h.bak

# just to be sure
/bin/rm -rf configure.log* RDict.* lib RESYNC PENDING
# eliminate .html files in src/contrib
/usr/bin/find src/contrib -type f -name "*.html" -exec /bin/rm -f {} \;
# eliminate .pyc files if any
/usr/bin/find . -type f -name "*.pyc" -exec /bin/rm -f {} \;
# eliminate misc files mercurial leaves
/usr/bin/find . -type f -name "*.orig" -exec /bin/rm -f {} \;
# Create EMACS-TAGS
cd $tmpdir/petsc$version; ${PETSC_DIR}/bin/maint/generateetags.py

# Set the correct file permissions.
cd $tmpdir
chmod -R a-w petsc$version
chmod -R u+w petsc$version
chmod -R a+r petsc$version
find petsc$version -type d -name "*" -exec chmod a+x {} \;

# Now create the tar files
cd $tmpdir
/bin/tar -czf ~/petsc$version.tar.gz petsc$version

# create lite version
/bin/rm -rf $tmpdir/petsc$version/docs $tmpdir/petsc$version/zope $tmpdir/petsc$version/config/BuildSystem/docs $tmpdir/petsc$version/config/examples/old
find $tmpdir/petsc$version -type f -name "*.html" -exec rm {} \;
# recreate EMACS-TAGS [after deletion]
cd $tmpdir/petsc$version; ${PETSC_DIR}/bin/maint/generateetags.py

cd $tmpdir
/bin/tar -czf ~/petsc-lite$version.tar.gz petsc$version

#cleanup
/bin/rm -rf $tmpdir

echo "Ending date: `date`"
