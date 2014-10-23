#!/usr/bin/env bash

addDataNew.py  --simulation --time 0 --range 1:180 > pans.txt
plotAddSurfaces.py --data1 ./pans.txt --title1 K14_120gpm_73psi_uniform_coupled_fixed --max 3.5

# addData.py  --testData --test 34 --data ../data/bert-pdfs/ > testData.txt
# plotAddSurfaces.py --data1 ./testData.txt --title1 elo_11psi_2MW_test_34 --max 0.5
# addData.py  --testData --test 35 --data ../data/bert-pdfs/ > testData.txt
# plotAddSurfaces.py --data1 ./testData.txt --title1 elo_11psi_2MW_test_35 --max 0.5
# addData.py  --testData --test 36 --data ../data/bert-pdfs/ > testData.txt
# plotAddSurfaces.py --data1 ./testData.txt --title1 elo_11psi_2MW_test_36 --max 0.5
