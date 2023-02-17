#!/bin/sh

# ADFRsuite installation script
pythonargs=" "
pyoptimize=0
TarDir=`pwd`
export ADS_ROOT=""
usage() {
    echo "Usage: ./install.sh [-d InstDir] [-c optimization]"
    exit 
}
# Parse the command-line arguments
opts=`getopt "hlc:d:" "$@"`
if [ "$?" != 0 ]
then
    usage
fi

set -- $opts

while true; do
    case "$1" in 

    -c) pythonargs="$pythonargs -c"; pyoptimize="$2"; shift; shift ;;
    -d) export ADS_ROOT="$2"; shift ; shift ;;
    -l) pythonargs="$pythonargs -l"; shift ;;
     
    -h) echo "Optional parameters:"
    echo "[-h]  help message;"
    echo "[ -d  InstDir] specifies installation directory (default-current directory)"
    echo "[ -c optimization] compile Python code with or without optimization:"
    echo "    0 - no optimization (generates .pyc files)"
    echo "    1 - with optimization (generates .pyo files);"
    exit ;;
    --) shift;;
    *)  if [ -z "$1" ] ; then break ; else echo "$1 is not a valid option" ; usage; fi ;; 
    esac
done


#echo "script options" python args "'$pythonargs'"  dest "'$ADS_ROOT'" pyoptimize $pyoptimize

if [ "$ADS_ROOT" != "" ]; then
    # check if the user has write access to the installation directory
    if [ -e "$ADS_ROOT" ]; then
	if [ -d "$ADS_ROOT" ]; then
	    if [ ! -w  "$ADS_ROOT" ]; then 
		echo "Can not complete installation - specified directory $ADS_ROOT does not have write access."
		exit 1

	    fi
	else 
	    echo "$ADS_ROOT" is not a directory
	    exit 1
	fi
    else 
	echo Creating directory "$ADS_ROOT"
	mkdir  -p "$ADS_ROOT"
    fi

else
    export ADS_ROOT="$(pwd)"
fi

echo "Installing ADFRsuite to $ADS_ROOT"

cd "$ADS_ROOT"
echo "Installing Python Interpreter to $ADS_ROOT"
tar xzvf $ADS_ROOT/Python*.tar.gz

if [ "$?" != 0 ]; then
    echo "Error in Python installation"
    exit 1
fi
echo Python installed, please wait for the rest of ADFRsuite to be installed 

cd "$ADS_ROOT"

## plaform we run on

export ADS_ARCHOSV=`$ADS_ROOT/Tools/archosv`

## add the path to the directory holding the python interpreter to your path

export PATH="$ADS_ROOT/bin:"$PATH

## use Python interpreter locally installed

PYTHON="$ADS_ROOT/bin/python2.7"
export PYTHONHOME="$ADS_ROOT"
if [ "`uname -s`" = "Linux" ] ; then
    export LD_LIBRARY_PATH="$ADS_ROOT/lib"
fi

## run python script - install.py - to install the packages and create scripts

if [ "$pyoptimize" -eq 1 ]; then
    echo "Running $PYTHON -O Tools/install.py $pythonargs"
    $PYTHON -O Tools/install.py $pythonargs
else
    echo "Running $PYTHON Tools/install.py $pythonargs"
    $PYTHON Tools/install.py  $pythonargs
fi

unset PYTHONHOME
