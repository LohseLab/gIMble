#!/usr/bin/env sh
./gIMble setup -s data/test.samples.csv -v data/test.vcf -b data/test.bed -o test -g data/test.genomefile && \
./gimble blocks -z test.z -l 10 -m 20 && \
./gimble blocks -z test.z -l 10 -m 20 -f && \
./gimble windows -w 25 -s 10 -z test.z && \
./gimble windows -w 25 -s 10 -z test.z -f
./gIMble simulate -z test.z -m data/sim.model.tsv -c data/sim.model.config.yaml
