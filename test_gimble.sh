#!/usr/bin/env sh
./gIMble setup -s data/test.samples.csv -v data/test.vcf -b data/test.bed -o test -g data/test.genomefile && \
./gimble blocks -z test.z -l 10 -r 2 -m 10 && \
./gimble blocks -z test.z -l 10 -r 2 -m 10 -f && \
./gimble windows -w 3 -s 1 -z test.z && \
./gimble windows -w 3 -s 1 -z test.z -f
