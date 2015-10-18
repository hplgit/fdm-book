#!/bin/sh
# Compile DocOnce errata document $name.do.txt to PDF and HTML formats
name=erratalist

doconce format pdflatex $name
doconce ptex2tex $name
pdflatex $name

doconce format html $name --html_style=bloodish

# Publish
dest=../../pub  # local dir (gh-pages on github)
dest=/some/repo/for/publishing/html   # external repo
if [ ! -d $dest ]; then
exit 0
fi
cp $name.pdf $name.html $dest

