#!/bin/bash

set -eu

mydir=$(dirname "$0")

source_url=https://code.soundsoftware.ac.uk/attachments/download/1698/Zweieck-Duell.ogg

testfile="$mydir/input.ogg"
outfile="$mydir/output.csv"
expfile="$mydir/expected.csv"

if sonic-annotator -v >/dev/null ; then
    :
else
    echo "Failed to find required binary sonic-annotator"
    exit 1
fi

if wget --version >/dev/null ; then
    wget -O "$testfile" "$source_url"
else
    curl -o "$testfile" "$source_url"
fi

VAMP_PATH="$mydir/.." \
	 sonic-annotator \
	 -d vamp:nnls-chroma:chordino:simplechord \
	 -w csv \
	 --csv-omit-filename \
	 --csv-one-file "$outfile" \
	 --csv-force \
	 "$testfile"

if cmp "$outfile" "$expfile" ; then
    echo 
    echo PASS
    exit 0
else
    echo
    echo "*** FAIL: Result does not match expected output. Diff follows:"
    echo
    diff "$outfile" "$expfile"
    exit 1
fi


