#! bin/bash

scons

files="data_sets/snap/as-skitter data_sets/snap/ca-AstroPh data_sets/snap/ca-CondMat data_sets/snap/ca-GrQc data_sets/snap/ca-HepPh data_sets/snap/ca-HepTh data_sets/law/cnr-2000 data_sets/law/dblp-2010 data_sets/law/dblp-2011 delaunay_n15 data_sets/snap/email-Enron data_sets/snap/email-EuAll data_sets/law/eu-2005 data_sets/konect/flickr-growth data_sets/konect/flickr-links data_sets/law/hollywood-2009 data_sets/law/hollywood-2011 data_sets/law/in-2004 data_sets/konect/libimseti data_sets/konect/orkut-links data_sets/snap/p2p-Gnutella04 data_sets/snap/p2p-Gnutella05 data_sets/snap/p2p-Gnutella06 data_sets/snap/p2p-Gnutella08 data_sets/snap/p2p-Gnutella24 data_sets/snap/p2p-Gnutella25 data_sets/snap/p2p-Gnutella30 data_sets/snap/p2p-Gnutella31 data_sets/konect/petster-carnivore data_sets/konect/petster-friendships-cat-uniq data_sets/konect/petster-friendships-dog-uniq rgg_n_2_15_s0 data_sets/snap/roadNet-CA data_sets/snap/roadNet-PA data_sets/snap/roadNet-TX data_sets/snap/soc-Epinions1 data_sets/snap/soc-LiveJournal1 data_sets/snap/soc-pokec-relationships data_sets/snap/soc-Slashdot0811 data_sets/snap/soc-Slashdot0902 data_sets/snap/web-Google data_sets/snap/web-NotreDame data_sets/snap/web-Stanford data_sets/snap/wiki-Talk data_sets/snap/wiki-Vote data_sets/konect/youtube-links data_sets/konect/youtube-u-growth data_sets/konect/zhishi-baidu-internallink data_sets/konect/zhishi-baidu-relatedpages data_sets/konect/zhishi-hudong-internallink"
for file in $files
do
    ./optimized/vcc --preconfiguration=fsocial --k=2 examples/$file.graph
done

echo "Done!"
