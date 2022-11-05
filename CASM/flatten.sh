#! /usr/bin/env bash

for d in */ ; do

    for FILE in $d*; do
        echo $FILE
        mv $FILE .
done
# Remove dir

    rm -r $d

done

