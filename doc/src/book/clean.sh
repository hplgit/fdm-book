#!/bin/sh
python -c 'import scripts; scripts.clean()'
rm -rf sphinx-* *.pyc automake* decay-book*
