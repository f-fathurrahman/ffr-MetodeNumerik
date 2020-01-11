#!/bin/bash

./do_pweave.sh Test1.texw
./do_pweave.sh Test2.texw

lualatex --shell-escape BUKU_METODE_NUMERIK.tex
lualatex --shell-escape BUKU_METODE_NUMERIK.tex