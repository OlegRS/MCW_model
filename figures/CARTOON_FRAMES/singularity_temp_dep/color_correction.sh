#!/bin/bash
# for file in *.png; do
#     convert file -resize 500 file;
# done

find ./ -name "*png" | while read filename; do convert $filename -fuzz 10%  -fill '#D8E219FF' -opaque '#3B528BFF' $filename; done
