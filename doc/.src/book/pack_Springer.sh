#!/bin/bash
# Pack a book project for Springer
set -x

author_name=langtangen
# Name of main text file in this directory
name=book
# Name of main tex file for the book as Springer will see
# it in subdir $author_name
book=book

# Put all files in directory $author_name
rm -rf $author_name
mkdir $author_name
cd $author_name
# Copy the single tex file for the entire book
cp ../${name}.tex $book.tex

# Copy all figures to one directory
mkdir figs
for dir in ../fig-*; do
  cp $dir/* figs
done
doconce subst '\{fig-.+?/' '{figs/' $book.tex

# Copy my hacked style files and give them new name with my initials
cp ~/texmf/tex/latex/misc/t2do.sty t2_hpl.sty
cp ~/texmf/tex/latex/misc/svmonodo.cls svmono_hpl.cls
doconce replace '{t2do}' '{t2_hpl}' $book.tex
doconce replace '{svmonodo}' '{svmono_hpl}' $book.tex
doconce subst '% Use .+ with doconce modifications.*' '' book.tex

# Copy ready-made discription of how this dir is organized
cp ../README_Springer_dir.txt 00README.txt

# Copy .bib file and newcommands
cp ../../chapters/papers.bib .
doconce replace '{../chapters/papers}' '{papers}' $book.tex
cp ../newcommands_keep.tex .

# Test that the book can be compiled in this subdir
rm -rf tmp.txt
pdflatex book | tee tmp.txt   # output of command in tmp.txt
rm -rf *.dvi *.aux *.out *.log *.loe *.toc

# Copy the log file from last run in the parent directory
# and analyze potential problems (too long lines, etc.) with the script
# doconce latex_problems
cp ../${name}.log book_last_run.log
doconce latex_problems book_last_run.log > book_last_run.log.summary

# Copy all style files needed for the book to a subdir stylefiles.
# Make use of 1) doconce's output of all needed style files (found
# in tmp.txt from running pdflatex above, but run doonce grab to
# get just the list of files), 2) the script ../grab_stylefiles.py
# to find each style file (the script creates a new script tmpcp.sh
# for copying the style files from their various locations).
doconce grab --from- '\*File List\*' --to- '\*\*\*\*' tmp.txt > tmp2.txt
python ../grab_stylefiles.py tmp2.txt  # make script tmpcp.sh
mkdir stylefiles
cd stylefiles
sh ../tmpcp.sh  # copy all style files
cd ..
rm tmpcp.sh
rm *~ tmp*

# Use most recently compiled PDF in the parent dir as official PDF
cp ../${name}.pdf $book.pdf

echo 'Pause before sending tar file to Springer...'
sleep 2
cd ..
# Springer's FTP site info
user=b3143
passwd=spvb3143
url=213.71.6.142
#--------------------------------------------
# Mark name of tar file with the date
file=TCSE6_Mar4_2031.tar.gz
#--------------------------------------------
tar czf $file $author_name
ftp -n $url <<EOF
quote $user
quote $passwd
cd langtangen
put $file
EOF

# cp a copy of the tex file and the tar file to my Google Drive
# in case Springer wants to get the files from there
#cp langtangen/$book.tex $file "/mnt/hgfs/hpl/Google Drive/Springer-TCSE6"
