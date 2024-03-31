#!/bin/bash -x
do="python3 ../../src/prepare_annotation.py"

# Should pass
$do vector-ok.bed reference-ok.tsv -o anno-simple-ok.txt && echo "Success"
$do vector-ok.bed empty.txt -o anno-simple-v-empty.txt && echo "Success"
$do vector-etc.bed reference-ok.tsv -o anno-multi-ok.txt && echo "Success"

# Should fail
$do no-vector.bed reference-ok.tsv -o fail.txt || echo "Failed successfully"
$do vector-ok.bed reference-bad-vector.tsv -o fail.txt || echo "Failed successfully"
$do vector-ok.bed reference-bad-dupes.tsv -o fail.txt || echo "Failed successfully"
$do empty.txt reference-ok.tsv -o fail.txt || echo "Failed successfully"

