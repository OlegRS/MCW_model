#!/bin/bash

find ./ -name "*png" | while read filename; do convert $filename -resize 700 $filename; done
    
