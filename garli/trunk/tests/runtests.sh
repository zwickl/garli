#!/bin/bash

if [ ! $# -eq 3 ];
then
	echo Usage: pass exactly three arguments:
	echo        '$0 <path to tests directory with data subdirectory> <path of GARLI binary> <path of NEXUSvalidator, part of NCL installation>'
	exit 1
fi

TESTS_DIR=$1
GARLI_BIN=$2
NEXUS_VAL=$3

rm  -f *.log00.log *.screen.log *.best*.tre *.best*.tre.phy *.boot.tre *.boot.phy *treelog00.tre *treelog00.log *problog00.log *fate00.log .*lock* *swaplog* *.check out.* qout.* mpi_m* *SiteLikes.log *sitelikes.log *best.all.phy *best.phy *current.phy *internalstates.log

echo "Linking to data ...."
ln -sf $TESTS_DIR/data | exit 1

echo "**************************"
echo "Running internal tests ..."
echo "**************************"

for i in $TESTS_DIR/internal/*.conf
do
        base=${i/*\/}
        base=${base/.conf/}
        echo "Running test $i"

        $GARLI_BIN -t $i || exit 1
done

echo "**************************"
echo "Running scoring tests ..."
echo "**************************"
line=1
for i in $TESTS_DIR/scoring/a.conf $TESTS_DIR/scoring/a.G3.conf $TESTS_DIR/scoring/a.G4.conf $TESTS_DIR/scoring/c.conf $TESTS_DIR/scoring/c.M3x2.conf $TESTS_DIR/scoring/n.conf $TESTS_DIR/scoring/n.G4.conf $TESTS_DIR/scoring/n.G5.conf 
do
	base=${i/*\/}
	base=${base/.conf/}

	$GARLI_BIN $i || exit 1
#	score=`tail -1 ${i%.conf}.sitelikes.log | awk '{print $2}'`
	score=`tail -1 scr.$base.sitelikes.log | awk '{print $2}'`
	expect=`head -n$line data/expected.scr | tail -n1`
	diff=`echo \($score\) - \($expect\) | bc`
	echo ***********TEST**************
	echo ***Score is $score
	echo ***Excpected is $expect
	echo ***SCORE DIFFERENCE IS $diff
	OK=`echo "$diff < 0.01 && $diff > -0.01" | bc`

	if (( $OK ))
	then
		echo "***Scoring OK for $i ***"
	else
		echo "***Scoring test failed for $i ***"
		exit 1
	fi
	line=`expr $line + 1`
done

echo "**************************"
echo "Running constraint tests ..."
echo "**************************"

for i in $TESTS_DIR/const/*.conf
do
        base=${i/*\/}
        base=${base/.conf/}
        echo "Running test $base"

        $GARLI_BIN $i || exit 1

        #NEXUSvalidator gives a warning every time it reads a tree file
        #without a taxa block.  So, shut it up initially and then if it
        #fails let it output whatever error
        $NEXUS_VAL $base*.tre 2> /dev/null
        if [ $? -eq 0 ]
        then
                echo TREEFILES ARE OK
        else
                $NEXUS_VAL $base*.tre || exit 1
        fi
done

echo "**************************"
echo "Running output tests ..."
echo "**************************"

for i in $TESTS_DIR/output/*.conf
do
	base=${i/*\/}
	base=${base/.conf/}
	echo "Running test $base"

	$GARLI_BIN $i || exit 1

	#NEXUSvalidator gives a warning every time it reads a tree file
        #without a taxa block.  So, shut it up initially and then if it
        #fails let it output whatever error
        $NEXUS_VAL $base*.tre 2> /dev/null
        if [ $? -eq 0 ]
        then
                echo TREEFILES ARE OK
        else
                $NEXUS_VAL $base*.tre || exit 1
        fi
done


echo "**************************"
echo "Running checkpoint tests ..."
echo "**************************"

for i in $TESTS_DIR/check/*.conf
do
	base=${i/*\/}
	base=${base/.conf/}
	echo "Running test $i"
	
	$GARLI_BIN $i || exit 1

	$GARLI_BIN $TESTS_DIR/restart/$base.conf || exit 1

	#NEXUSvalidator gives a warning every time it reads a tree file
        #without a taxa block.  So, shut it up initially and then if it
        #fails let it output whatever error
        $NEXUS_VAL ch.$base*.tre 2> /dev/null
        if [ $? -eq 0 ]
        then
                echo "TREEFILES ARE OK"
        else
                $NEXUS_VAL ch.$base*.tre || exit 1
        fi
done

