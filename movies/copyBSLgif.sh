#!/bin/bash

BSLPATH=michai@bslcenc.nerc-bas.ac.uk:/users/michai/mitgcmPy

scp ${BSLPATH}/$1.gif .

ffmpeg -i $1.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $1.mp4

rm $1.gif
