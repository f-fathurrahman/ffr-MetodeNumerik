#!/bin/bash

#./do_pweave.sh Bisection.texw
#./do_pweave.sh RegulaFalsi.texw
./do_pweave.sh FixedPoint.texw

lualatex --shell-escape BUKU_METODE_NUMERIK.tex
lualatex --shell-escape BUKU_METODE_NUMERIK.tex