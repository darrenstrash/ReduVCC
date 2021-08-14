#! bin/bash 

snap="as-skitter email-Enron roadNet-CA roadNet-PA roadNet-TX soc-Epinions1 soc-LiveJournal1 soc-pokec-relationships soc-Slashdot0902 web-BerkStan web-Google web-NotreDame web-Stanford wiki-Talk" 

for file in $snap
do 
	qsub $1/job_$file.run
done

konect="flickr-growth flickr-links libimseti petster-carnivore petster-friendships-cat-uniq petster-friendships-dog-uniq youtube-links youtube-u-growth zhishi-baidu-internallink zhishi-baidu-relatedpages zhishi-hudong-internallink"

for file in $konect
do 
        qsub $1/job_$file.run
done

law="cnr-2000 dblp-2011 eu-2005 hollywood-2009 hollywood-2011 in-2004 indochina-2004 uk-2002"
for file in $law
do 
        qsub $1/job_$file.run
done
