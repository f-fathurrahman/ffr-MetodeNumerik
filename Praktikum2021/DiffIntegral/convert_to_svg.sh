#!/bin/bash
BASNAM=`basename $1 .tex`
pdflatex $1
pdfcrop $BASNAM.pdf $BASNAM.pdf
pdf2svg $BASNAM.pdf $BASNAM.svg