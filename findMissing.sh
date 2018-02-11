#!/bin/sh

while grep -q False result/ex_$1.txt; do
    python2 ff_pcn/database.py result/ex_$1.txt;
done
