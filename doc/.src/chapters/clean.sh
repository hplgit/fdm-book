#!/bin/sh -x
doconce clean

rm -rf .ptex2tex.cfg* newcommands*.tex automake_sphinx.py *-solarized*.* *~ slides-decay/*~ .*.exerinfo tmp* _doconce_deb* *-beamer*.pdf *-4print.pdf *-4screen.pdf *.html *-plain.* \#* libpeerconnection.log *-2up.pdf *-A4*.pdf deck.js reveal.js latex_figs html_images slides_*.tex *-sol.*

if [ -d src-$1 ]; then
  cd src-$1
  if [ -f clean.sh ]; then
    sh clean.sh
  fi
  rm -f *.pyc
fi
cd ..
