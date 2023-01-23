#! /bin/bash
biosample=$1
rep1=$2
rep2=$3
exp=$4
cd $biosample
# removed.bedpe are ENCODE bedpe files removing interchromosomal links and links with overlapping anchors (seven column)
awk '{if ($5-$3>8000) print $0 "\t" 0}' $rep1.removed.bedpe > $rep1.removed+.bedpe
awk '{if ($5-$3>8000) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" 0 "\t" $7}' $rep2.removed.bedpe > $rep2.removed+.bedpe
bedtools pairtopair -a $rep1.removed+.bedpe -b $rep2.removed+.bedpe -type both | sort -k1,1 -k2,3 -k5,6n > repLinks.bedpe
cur=0
uniq=1
# merge replicated links between reps
while [ $cur != $uniq ]
do
	 awk  'BEGIN{c1=0;c2=0;}{
		if ( !rep1[$1,$2,$3,$4,$5,$6] && !rep2[$9,$10,$11,$12,$13,$14])
			{
			if (c1+c2>0) print chr "\t" start1 "\t" end1 "\t" chr "\t" start2 "\t" end2 "\t" c1 "\t" c2 "\t";
			rep1[$1,$2,$3,$4,$5,$6]++;
			rep2[$9,$10,$11,$12,$13,$14]++;
			chr=$1
			if ($2<$10) start1=$2;
			else start1=$10;
			if ($3>$11) end1=$3;
			else end1=$11;
			if ($5<$13) start2=$5;
			else start2=$13;
			if ($6>$14) end2=$6;
			else end2=$14;
			if ($1!=$9 || $2!=$10 || $3!=$11 || $4!=$12 || $5!=$13 || $6!=$14 || $7!=$15 || $8!=$16){
				c1=$7+$15;
				c2=$8+$16;
			}
			else {
			        c1=$7;
				c2=$8;
			}
		 }
	         if ( rep1[$1,$2,$3,$4,$5,$6]>0 && !rep2[$9,$10,$11,$12,$13,$14] )
			{
			rep2[$9,$10,$11,$12,$13,$14]++;
			if (start1>$10) start1=$10;
			if (end1<$11) end1=$11;
			if (start2>$13) start2=$13;
			if (end2<$14) end2=$14;
			if ($1!=$9 || $2!=$10 || $3!=$11 || $4!=$12 || $5!=$13 || $6!=$14 || $7!=$15 || $8!=$16){
				c1=c1+$15;
				c2=c2+$16;
			}
		}
		if ( !rep2[$1,$2,$3,$4,$5,$6] && !rep1[$1,$2,$3,$4,$5,$6] && rep2[$9,$10,$11,$12,$13,$14]>0 )
		{
			if (c1+c2>0) print chr "\t" start1 "\t" end1 "\t" chr "\t" start2 "\t" end2 "\t" c1 "\t" c2 "\t";
			rep1[$1,$2,$3,$4,$5,$6]++;
			start1=$2;
			end1=$3;
			start2=$5;
			end2=$6;
			c1=$7
			c2=$8
		}
        } END{print chr "\t" start1 "\t" end1 "\t" chr "\t" start2 "\t" end2 "\t" c1 "\t" c2 "\t";}' repLinks.bedpe > out
# check for links with overlapping anchors
 	echo "number of overlapping links with overlapping anchors"
 	awk '{if ($3>$5) print $0;}' out | wc -l
	echo "removed links with overlapping anchors"
	awk '{if ($3<=$5) print $0;}' out > out2
 	bedtools pairtopair -a out2 -b out2 -type both | sort -k1,1 -k2,3n -k5,6n > repLinks.bedpe
 	cur=$(wc -l out2 | awk '{print $1}')
	uniq=$(wc -l repLinks.bedpe | awk  '{print $1}')
done
echo $cur
rep1=$(awk 'BEGIN{count=0}{count=count+$7;}END{print count;}' ${rep1}.removed+.bedpe | awk '{print $1/1000000}')
rep2=$(awk 'BEGIN{count=0}{count=count+$8;}END{print count;}' ${rep2}.removed+.bedpe | awk '{print $1/1000000}')
awk -v rep1=$rep1 -v rep2=$rep2 '{print $0 $7/rep1 "\t" $8/rep2 "\t" ($7/rep1+$8/rep2)/2}' out2 > $exp.bedpe

