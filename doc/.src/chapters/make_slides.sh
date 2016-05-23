#!/bin/bash
#
# bash ../make_slides.sh slides_nickname.do.txt
#
# But this script is normally run from make.sh to make both
# chapter and slides.

dofile=$1
if [ ! -f $dofile ]; then
  echo "No such file: $dofile"
  exit 1
fi
nickname=$2

filename=`echo $dofile | sed 's/\.do\.txt//'`

function system {
  output=`"$@"`
  if [ $? -ne 0 ]; then
    echo "make.sh: unsuccessful command $@"
    echo "abort!"
    exit 1
  fi
  echo "$output" > tmp.log
}

# LaTeX PDF for printing
# (Smart to make this first to detect latex errors - HTML/MathJax
# gives far less errors and warnings)
rm -f *.aux
preprocess -DFORMAT=pdflatex ../newcommands_keep.p.tex > newcommands_keep.tex
system doconce format pdflatex $filename --latex_admon=paragraph --latex_code_style=lst-yellow2
doconce replace 'section{' 'section*{' ${filename}.tex
system pdflatex $filename
system pdflatex $filename
mv -f $filename.pdf ${filename}-plain.pdf

# HTML
preprocess -DFORMAT=html ../newcommands_keep.p.tex > newcommands_keep.tex

# reveal.js HTML5 slides
html=${filename}-reveal
system doconce format html $filename --pygments_html_style=native --keep_pygments_html_bg --html_output=$html
doconce replace 'P$d$' 'Pd' ${html}.html
system doconce slides_html ${html}.html reveal --html_slide_theme=darkgray

html=${filename}-reveal-beige
system doconce format html $filename --pygments_html_style=perldoc --keep_pygments_html_bg --html_output=$html
doconce replace 'P$d$' 'Pd' ${html}.html
system doconce slides_html ${html}.html reveal --html_slide_theme=beige

# deck.js HTML5 slides
#html=${filename}-deck
#system doconce format html $filename --pygments_html_style=perldoc --keep_pygme#nts_html_bg --html_output=$html
#doconce replace 'P$d$' 'Pd' ${html}.html
#system doconce slides_html ${html}.html deck --html_slide_theme=sandstone.default

# Plain HTML with everything in one file
html=${filename}-1
system doconce format html $filename --html_style=bloodish --html_output=$html -DWITH_TOC
doconce replace 'P$d$' 'Pd' ${html}.html
doconce replace "<li>" "<p><li>" ${html}.html
doconce split_html ${html}.html --method=space8

# HTML with solarized style and one big file
html=${filename}-solarized
system doconce format html $filename --html_style=solarized3 --html_output=$html --pygments_html_style=perldoc --pygments_html_linenos  -DWITH_TOC
doconce replace 'P$d$' 'Pd' ${html}.html
doconce replace "<li>" "<p><li>" ${html}.html
doconce split_html ${html}.html --method=space8

# Drop split HTML files - too much hassle with cranking up the font
# for each slide...

# LaTeX Beamer
rm -f *.aux
preprocess -DFORMAT=pdflatex ../newcommands_keep.p.tex > newcommands_keep.tex
system doconce format pdflatex $filename --latex_title_layout=beamer --latex_table_format=footnotesize --latex_admon_title_no_period --latex_code_style=pyg --movie_prefix=https://raw.githubusercontent.com/hplgit/fdm-book/master/doc/.src/chapters/${nickname}/[[[
system doconce slides_beamer $filename --beamer_slide_theme=red_shadow
system pdflatex -shell-escape $filename
pdflatex -shell-escape $filename
pdflatex -shell-escape $filename
cp ${filename}.pdf ${filename}-beamer.pdf
rm -f ${filename}.pdf
cp ${filename}.tex ${filename}-beamer.tex  # sometimes nice to look at

# Handouts
system doconce format pdflatex $filename --latex_title_layout=beamer --latex_table_format=footnotesize --latex_admon_title_no_period --latex_code_style=pyg
system doconce slides_beamer $filename --beamer_slide_theme=red_shadow --handout
system pdflatex -shell-escape $filename
pdflatex -shell-escape $filename
pdflatex -shell-escape $filename
pdfnup --nup 2x3 --frame true --delta "1cm 1cm" --scale 0.9 --outfile ${filename}-beamer-handouts2x3.pdf ${filename}.pdf
rm -f ${filename}.pdf

# Ordinary plain LaTeX
rm -f *.aux  # important
system doconce format pdflatex $filename --latex_admon=paragraph --latex_code_style=lst-yellow2
doconce replace 'section{' 'section*{' ${filename}.tex
system pdflatex $filename
system pdflatex $filename
mv -f $filename.pdf ${filename}-plain.pdf
