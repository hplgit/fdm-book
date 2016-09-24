reports="report1and2"
for report in $reports; do
doconce format pdflatex $report --latex_code_style=vrb --latex_papersize=a4
doconce replace 'report on' 'report on \\ ' $report.tex
pdflatex $report
done
