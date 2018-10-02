#!/bin/bash

for value in {1..5}
do
  octave --eval "SparkFrequency($1, $2, $value, $3, $4)"
done

