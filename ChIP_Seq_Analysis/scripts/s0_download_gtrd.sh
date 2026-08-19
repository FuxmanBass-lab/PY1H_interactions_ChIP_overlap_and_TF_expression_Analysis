#!/bin/sh
# Downloads GTRD ChIP-seq peak files and the sample metadata table.

wget -r -np -nH --cut-dirs=3 -R index.html* http://gtrd.biouml.org:8888/egrid/bigBeds/hg38/ChIP-seq/Peaks/MACS2/
wget http://gtrd.biouml.org:8888/downloads/current/metadata/ChIP-seq.metadata.txt
