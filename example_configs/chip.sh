#!/bin/bash

curl -s https://raw.githubusercontent.com/biomedbigdata/TF-Prioritizer/master/docker.py | python3 - -c chipSeq.json -o /nfs/data/COM2POSE/ATAC-seq/out_chip -t 30 -m 90
