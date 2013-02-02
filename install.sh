##
## Installation script relevant only for users who want
## By Marcus Breese
##
#!/usr/bin/env bash
VIRTUALENV=`which virtualenv-2.7`
if [ "$VIRTUALENV" == "" ]; then
    VIRTUALENV=`which virtualenv-2.6`
fi
if [ "$VIRTUALENV" == "" ]; then
    VIRTUALENV=`which virtualenv`
fi
if [ "$VIRTUALENV" == "" ]; then
    echo "Missing virtualenv!"
    exit 1
fi

if [ ! -e env ]; then 
    echo "Initializing virtualenv folder (env)"
    $VIRTUALENV env
fi
. env/bin/activate
python setup.py install --prefix env
pip install numpy
pip install -r requirements.txt

#
# If you're on a Mac, and matplotlib won't install, try this:
# LDFLAGS="-L/usr/X11/lib" CFLAGS="-I/usr/X11/include -I/usr/X11/include/freetype2 -I/usr/X11/include/libpng15" pip install matplotlib