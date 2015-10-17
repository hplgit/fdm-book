#!/bin/bash
# Generic make.sh
# Usage: make.sh name [sphinx publish src]
#
# if publish: -DNOTREAD preprocess option is added so that all
# sections are included. Use #ifdef NOTREAD to mark sections
# that should not yet be published because they are not sufficiently
# proof read.

set -x

nickname=$1
mainname=main_$nickname

function system {
  "$@"
  if [ $? -ne 0 ]; then
    echo "make.sh: unsuccessful command $@"
    echo "abort!"
    exit 1
  fi
}


if [ $# -eq 0 ]; then
  print 'name of document missing!'
  exit 1
fi

function oliver_art {
if [ -f fig-${nickname}/oliver_sin.jpg ]; then
doconce subst "<h2>Table of contents</h2>" "<h2>Table of contents</h2>\n\n<a href=\"https://www.facebook.com/oliversinofficial\"><img src=\"fig-$name/oliver_sin.jpg\" width=300 align=\"right\"></a>" ${mainname}.html
cp fig-${nickname}/oliver_sin.jpg sphinx-rootdir/_build/html/_images
doconce subst '<div class="toctree-wrapper compound">' '<div class="toctree-wrapper compound">\n\n<a href="https://www.facebook.com/oliversinofficial"><img src="_images/oliver_sin.jpg" width=300 align="right"></a>\n' sphinx-rootdir/_build/html/index.html
fi
}

rm -f tmp_*

doconce spellcheck -d .dict4spell.txt *.do.txt
if [ $? -ne 0 ]; then
  echo "make.sh aborts due to misspellings"
  exit 1
fi
rm -rf tmp_stripped*

egrep "[^\\]thinspace" *.do.txt
if [ $? -eq 0 ]; then echo "wrong thinspace commands - abort"; exit; fi

#comments="--skip_inline_comments"
comments=""
doc=document
appendix=document

preprocessor_opt="DOCUMENT=$doc APPENDIX=$appendix BOOK=standalone"
no_solutions='--without_solutions --without_answers'
sphinx=0
publish=0

if [ $# -ge 2 ]; then
  sphinx=1
fi
if [ $# -ge 3 ]; then
  publish=1
else
  publish=0
  # allow everything to be compiled...
  preprocessor_opt="${preprocessor_opt} -DNOTREAD"
fi

rm -f *.aux
preprocess -DFORMAT=pdflatex ../newcommands_keep.p.tex > newcommands_keep.tex

# PDF with solutions (start with this and let .aux be the most
# relevant version for xr references
system doconce format pdflatex ${mainname} $preprocessor_opt $comments --latex_table_format=center --device=screen "--latex_code_style=default:vrb-blue1@sys:vrb[frame=lines,label=\\fbox{{\tiny Terminal}},framesep=2.5mm,framerule=0.7pt,fontsize=\fontsize{9pt}{9pt}]"
doconce latex_exercise_toc ${mainname}
doconce subst 'frametitlebackgroundcolor=.*?,' 'frametitlebackgroundcolor=blue!5,' ${mainname}.tex
system pdflatex ${mainname}
makeindex ${mainname}
bibtex ${mainname}
pdflatex ${mainname}
pdflatex ${mainname}
cp ${mainname}.pdf ${nickname}-sol.pdf
rm -f ${mainname}.pdf

# See http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/firth/software/pdfjam for examples on pdfnup and pdfjam

# PDF for printing, standard paper size
system doconce format pdflatex ${mainname} $preprocessor_opt $solutions $comments --latex_table_format=center --device=paper "--latex_code_style=default:vrb-blue1@sys:vrb[frame=lines,label=\\fbox{{\tiny Terminal}},framesep=2.5mm,framerule=0.7pt,fontsize=\fontsize{9pt}{9pt}]"
doconce subst 'frametitlebackgroundcolor=.*?,' 'frametitlebackgroundcolor=blue!5,' ${mainname}.tex
pdflatex ${mainname}
makeindex ${mainname}
bibtex ${mainname}
pdflatex ${mainname}
pdflatex ${mainname}
cp ${mainname}.pdf ${nickname}-4print.pdf

# 2 pages per sheet to save trees
#pdfnup --frame true --outfile ${nickname}-4print-2up.pdf ${nickname}-4print.pdf

# PDF for electronic viewing
system doconce format pdflatex ${mainname} $preprocessor_opt $no_solutions $comments --latex_table_format=center --device=screen "--latex_code_style=default:vrb-blue1@sys:vrb[frame=lines,label=\\fbox{{\tiny Terminal}},framesep=2.5mm,framerule=0.7pt,fontsize=\fontsize{9pt}{9pt}]"
doconce latex_exercise_toc ${mainname}
doconce subst 'frametitlebackgroundcolor=.*?,' 'frametitlebackgroundcolor=blue!5,' ${mainname}.tex
system pdflatex ${mainname}
makeindex ${mainname}
bibtex ${mainname}
pdflatex ${mainname}
pdflatex ${mainname}
cp ${mainname}.pdf ${nickname}-4screen.pdf
rm -f ${mainname}.pdf

# HTML
system preprocess -DFORMAT=html ../newcommands_keep.p.tex > newcommands_keep.tex

style=solarized
# Note that we use perldoc and not native solarized3 code
system doconce format html ${mainname} $preprocessor_opt $no_solutions --html_style=solarized3 --html_output=${nickname}-$style --pygments_html_style=perldoc $comments
doconce replace 'P$d$' 'Pd' ${nickname}-${style}.html
system doconce split_html ${nickname}-${style}.html --nav_button=text

# default nickname.html is bootswatch journal style
system doconce format html ${mainname} $preprocessor_opt $no_solutions $comments --html_style=bootswatch_journal --html_code_style=inherit --html_output=${nickname}
doconce replace 'P$d$' 'Pd' ${nickname}.html
system doconce split_html ${nickname}.html

# nickname-sol.html is bootswatch journal with solutions
system doconce format html ${mainname} $preprocessor_opt $comments --html_style=bootswatch_journal --html_code_style=inherit --html_output=${nickname}-sol
doconce replace 'P$d$' 'Pd' ${nickname}-sol.html
system doconce split_html ${nickname}-sol.html

if [ $sphinx -eq 1 ]; then
system doconce format sphinx ${mainname} $preprocessor_opt $no_solutions $comments --sphinx_keep_splits
doconce replace 'P$d$' 'Pd' ${mainname}.rst
system doconce split_rst ${mainname}
system doconce sphinx_dir theme=cbc copyright="H. P. Langtangen" ${mainname}
system python automake_sphinx.py
fi

# Add Oliver Sin's image for this chapter (if it exists), html and sphinx
oliver_art


# Note: this script is run from a subdirectory
if [ -f slides_${nickname}.do.txt ]; then
bash -x ../make_slides.sh slides_${nickname}.do.txt
if [ $? -ne 0 ]; then echo "doconce could not compile slides - abort"; exit; fi
fi
rsync="rsync -rtDvz -u -e ssh -b --exclude-from=$HOME/1/.rsyncexclude --delete --force "

#repo=~/vc/fdm-book/doc
repo=../../..
if [ $publish -eq 1 ]; then
# Copy compiled documents to destination repo
dest=$repo/pub
dest_sol=$repo/Trash
if [ ! -d $dest/${nickname} ]; then
echo "making directory ${dest}/${nickname}"
mkdir $dest/${nickname}
mkdir $dest/${nickname}/pdf
mkdir $dest/${nickname}/html
fi
if [ ! -d $dest_sol ]; then
    mkdir $dest_sol
fi
if [ ! -d $dest_sol/${nickname} ]; then
echo "making directory ${dest_sol}/${nickname}"
mkdir $dest_sol/${nickname}
mkdir $dest_sol/${nickname}/pdf
mkdir $dest_sol/${nickname}/html
fi

echo "copying compiled documents to $dest"
# Should not copy chapters since that requires ref[][][] and
# publishing of separate chapters...or should chapters exist in html?
cp ${nickame}*.pdf slides*.pdf ${dest}/${nickname}/pdf
cp ${nickname}*.html .${nickname}* ._${nickname}*.html slides*.html ._slides* ${dest}/${nickname}/html
# Remove documents with solutions
rm -f ${dest}/${nickname}/html/*-sol.html ${nickname}/html/._*-sol.html
rm -f ${dest}/${nickname}/pdf/*-sol.pdf

# Replace the sphinx doc entirely
rm -rf ${dest}/${nickname}/sphinx
cp -r sphinx-rootdir/_build/html ${dest}/${nickname}/sphinx

rm -f fig-${nickname}/*~ fig-${nickname}/tmp*
# drop copying deck.js
dirs="fig-${nickname} mov-${nickname} reveal.js"
for dir in $dirs; do
  if [ -d $dir ]; then
    $rsync $dir ${dest}/${nickname}/html
  fi
done

echo "copying compiled documents to $dest_sol"

cp ${nickname}-sol.pdf ${dest_sol}/${nickname}/pdf
cp ${nickname}-sol.html ._${nickname}-sol*.html ${dest_sol}/${nickname}/html
dirs="fig-${nickname} mov-${nickname}"
for dir in $dirs; do
  if [ -d $dir ]; then
    $rsync $dir ${dest_sol}/${nickname}/html
  fi
done

# Copy slides source
#dest=$repo/doc/slides/src/$name
#if [ ! -d $dest ]; then
#  mkdir $dest
#fi
#$rsync lecture_${nickname}.do.txt slides-$nickname mov-$name fig-$namename src-$nickname $dest
#cp ../newcommands_keep.p.tex ../make_lectures.sh ../mako_code.txt ../papers.pub ../venues.list $dest/..

cd ..
cp index_files.do.txt index.do.txt
system doconce format html index --html_style=bootstrap --html_links_in_new_window --html_bootstrap_navbar=off
cd -
cp ../index.html $dest
rm -f ../index.*

# Copy src
if [ $# -ge 3 ]; then
dest=$repo/src
if [ ! -d $dest/$nickname ]; then
  mkdir $dest/$nickname
fi
if [ -d src-$nickname ]; then
cd src-$nickname
if [ -f clean.sh ]; then
sh clean.sh
fi
cd ..
files=`find -L src-$nickname \( -name '*.py' -o -name '*.pyx' -o -name '*.f' -o -name '*.c' -o -name '*.cpp' -o -name 'make*.sh' \) -print`
if [ -f src-${nickname}/README ]; then
  files="$files src-$name/README"
fi
tar cfz tmp.tar.gz $files
mv -f tmp.tar.gz  $dest
cd $dest
rm -rf src-$nickname
tar xfz tmp.tar.gz
rm -f tmp.tar.gz
python ~/1/programs/rsync_git.py src-$nickname $nickname
rm -rf src-$nickname
cd -
fi
fi
fi

cd $dest
git add .
