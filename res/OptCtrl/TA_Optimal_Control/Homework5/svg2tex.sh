#! /bin/bash

#FILES= eval "ls $1 | grep .svg"
#FILENAMES="${FILES%.*}"
for f in $1/*.svg 
do 
#base="${f##*/}"
#base= "${f%.svg}"
 echo "inkscape -z -D --file=./$f --export-pdf=./${f%.svg}.pdf --export-latex"
eval "inkscape -z -D --file=./$f --export-pdf=./${f%.svg}.pdf --export-latex"
done
#echo $FILENAMES
#
#echo FILES