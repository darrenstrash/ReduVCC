#! bin/bash 

snap="as-skitter ca-AstroPh ca-CondMat ca-GrQc ca-HepPh ca-HepTh email-Enron email-EuAll p2p-Gnutella04 p2p-Gnutella05 p2p-Gnutella06 p2p-Gnutella08 p2p-Gnutella09 p2p-Gnutella24 p2p-Gnutella25 p2p-Gnutella30 p2p-Gnutella31 roadNet-CA roadNet-PA roadNet-TX soc-Epinions1 soc-LiveJournal1 soc-pokec-relationships soc-Slashdot0811 soc-Slashdot0902 web-BerkStan web-Google web-NotreDame web-Stanford wiki-Talk wiki-Vote" 

for file in $snap
do 
	qsub $1/job_$file.run
done

konect="flickr-growth flickr-links libimseti petster-carnivore petster-friendships-cat-uniq petster-friendships-dog-uniq youtube-links youtube-u-growth zhishi-baidu-internallink zhishi-baidu-relatedpages zhishi-hudong-internallink"

for file in $konect
do 
        qsub $1/job_$file.run
done

law="cnr-2000 dblp-2010 dblp-2011 eu-2005 hollywood-2009 hollywood-2011 in-2004 indochina-2004 uk-2002"
for file in $law
do 
        qsub $1/job_$file.run
done
