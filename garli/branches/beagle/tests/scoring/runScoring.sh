#!/bin/bash


echo $score

expect=14498.88532
echo $expect

line=1

for i in *.conf
do
	$TESTING_EXEC $i
	score=`tail -1 scr.${i%.conf}.sitelikes.log | awk '{print $2}'`
	expect=`head -n$line data/expected.scr | tail -n1`
	diff=`echo \($score\) - \($expect\) | bc`
	echo ***********TEST**************
	echo ***Score is $score
	echo ***Excpected is $expect
	echo ***SCORE DIFFERENCE IS $diff
	OK=`echo "$diff < 0.01 && $diff > -0.01" | bc`

	if (( $OK ))
	then
		echo ***Score OK***
	else
		echo ***Score Not OK***
		exit 1
	fi
	line=`expr $line + 1`
done

