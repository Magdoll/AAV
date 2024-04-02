#!/bin/bash -x
#
do="python3 ../../src/prepare_annotation.py"

# Should pass:
# - 1-line vector, empty others
$do vector-ok.bed reference-ok.tsv -o anno-simple-ok.txt && echo "Success"
# - 1-line vector, 1-line others with the same vector
$do vector-ok.bed empty.txt -o anno-simple-v-empty.txt && echo "Success"
# - multi-line vector with ITR labels (~140bp), multi-line others with the same vector and multi host
$do vector-etc.bed reference-ok.tsv -o anno-multi-ok.txt && echo "Success"

# Should fail:
# - 1-line vector with no 'vector' label, any others
$do no-vector.bed reference-ok.tsv -o fail.txt || echo "Failed successfully"
# - OK 1-line vector, 1-line others with different seq name for 'vector' label
$do vector-ok.bed reference-bad-vector.tsv -o fail.txt || echo "Failed successfully"
# - multi-line others with conflicting types for same seq name
$do vector-ok.bed reference-bad-dupes.tsv -o fail.txt || echo "Failed successfully"
# - empty vector file
$do empty.txt reference-ok.tsv -o fail.txt || echo "Failed successfully"
