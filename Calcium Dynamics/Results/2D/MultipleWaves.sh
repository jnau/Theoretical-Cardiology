#!/bin/bash

for value in {1..5}
do
  octave --eval "WaveSpeed($1, $2, $value)"
done

