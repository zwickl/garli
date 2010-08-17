#!/bin/bash

if [ $# -lt 3 ];
then
	echo Usage: pass three, or \(optionally\) more arguments:
	echo        '$0 <path to tests directory with data subdirectory> <path of GARLI binary> <path of NEXUSvalidator, part of NCL installation> [optional: GARLI command-line arguments]'
	exit 1
fi

TESTS_DIR=$1
GARLI_BIN=$2
NEXUS_VAL=$3

if [ $# -gt 3 ]
then
	shift; shift; shift;
	GARLI_ARGS=$@
fi

#set this to move on to the next test after failing one
#NO_EXIT_ON_ERR=1

rm  -f *.log00.log *.screen.log *.best*.tre *.best*.tre.phy *.boot.tre *.boot.phy *treelog00.tre *treelog00.log *problog00.log *fate00.log .*lock* *swaplog* *.check out.* qout.* mpi_m* *SiteLikes.log *sitelikes.log *best.all.phy *best.phy *current.phy *internalstates.log

echo "Linking to data ...."
if [ -d $TESTS_DIR ];then
	echo "data folder already exists"
else
	ln -sf $TESTS_DIR/data | exit 1
fi

echo "**************************"
echo "Running internal tests ..."
echo "**************************"

if [ -d $TESTS_DIR/internal ];then
	for i in $TESTS_DIR/internal/*.conf
	do
    		base=${i/*\/}
    		base=${base/.conf/}
    		echo "Running internal test $i"
		echo "Running internal test $i" >&2

	 	$GARLI_BIN -t $i $GARLI_ARGS
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi
	done
else
	echo "No internal tests found ..."
fi

echo "**************************"
echo "Running scoring tests ..."
echo "**************************"

if [ -d $TESTS_DIR/scoring ];then
#	line=1
	for i in $TESTS_DIR/scoring/a.conf $TESTS_DIR/scoring/a.G3.conf $TESTS_DIR/scoring/a.G4.conf $TESTS_DIR/scoring/c.conf $TESTS_DIR/scoring/c.M3x2.conf $TESTS_DIR/scoring/n.conf $TESTS_DIR/scoring/n.G4.conf $TESTS_DIR/scoring/n.G5.conf $TESTS_DIR/scoring/p.mk.conf $TESTS_DIR/scoring/p.mkO.conf $TESTS_DIR/scoring/p.mkO.ssr.conf $TESTS_DIR/scoring/p.mk.ssr.conf $TESTS_DIR/scoring/p.mkv.conf $TESTS_DIR/scoring/p.mkvO.conf $TESTS_DIR/scoring/p.mkvO.ssr.conf $TESTS_DIR/scoring/p.mkv.ssr.conf $TESTS_DIR/scoring/p.3diff.conf

	do
		if [ -f $i ];then
			base=${i/*\/}
			base=${base/.conf/}
	
			echo "Running scoring test $i"
			echo "Running scoring test $i" >&2
			$GARLI_BIN $i $GARLI_ARGS
			if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
				exit 1
			fi
	
			#figure out what precision we can expect
			if [ ! `grep "likelihood precision" scr.$base.screen.log | wc -l` -eq 0 ]
			then
				allowed=`grep "likelihood precision" scr.$base.screen.log | awk '{print $6}'`
			else
				#for partitioned models allow a bit more scoring leeway
				if test  "${base:0:2}" = "p." 
				then
					allowed=0.05
				else
					allowed=0.01
				fi
			fi

			score=`tail -1 scr.$base.sitelikes.log | awk '{print $2}'`
			expect=`grep $base.conf data/expected.scr | awk '{print $2}'`
			diff=`echo \($score\) - \($expect\) | bc`
			echo ***********TEST**************
			echo ***Score is $score
			echo ***Expected is $expect
			echo ***SCORE DIFFERENCE IS $diff
			echo ***ALLOWED ERROR IS $allowed
			OK=`echo "$diff < $allowed && $diff > -$allowed" | bc`
		
			if (( $OK ))
			then
				echo "***Scoring OK for $i ***"
			else
				echo "***Scoring test failed for $i ***"
				if [[ ! -n "$NO_EXIT_ON_ERR" ]];then
					exit 1
				fi
			fi
			line=`expr $line + 1`
		fi	
	done
else
	echo "No scoring tests found ..."
fi

echo "**************************"
echo "Running constraint tests ..."
echo "**************************"

if [ -d $TESTS_DIR/const ];then

	for i in $TESTS_DIR/const/*.conf
	do
		base=${i/*\/}
	        base=${base/.conf/}
	        echo "Running constraint test $base"
		echo "Running constraint test $base" >&2

	        $GARLI_BIN $i $GARLI_ARGS
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi

		#NEXUSvalidator gives a warning every time it reads a tree file
		#without a taxa block.  So, shut it up initially and then if it
		#fails let it output whatever error
		$NEXUS_VAL $base*.tre 2> /dev/null
		if [ $? -eq 0 ]
	    	then
	    		echo TREEFILES PASS
	    	else
	    		$NEXUS_VAL $base*.tre
			if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
				exit 1
			fi
	    	fi
	done
else
	echo "No constraint tests found ..."
fi

echo "**************************"
echo "Running output tests ..."
echo "**************************"

if [ -d $TESTS_DIR/output ];then

	for i in $TESTS_DIR/output/*.conf
	do
		base=${i/*\/}
		base=${base/.conf/}
		echo "Running output test $base"
		echo "Running output test $base" >&2

		$GARLI_BIN $i $GARLI_ARGS
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi

	#NEXUSvalidator gives a warning every time it reads a tree file
	#without a taxa block.  So, shut it up initially and then if it
	#fails let it output whatever error
	$NEXUS_VAL $base*.tre 2> /dev/null
	if [ $? -eq 0 ]
		then
		echo "TREEFILES PASS"
		else
			$NEXUS_VAL $base*.tre
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi
    	fi
	done
else
	echo "No output tests found ..."
fi

echo "**************************"
echo "Running checkpoint tests ..."
echo "**************************"

if [ -d $TESTS_DIR/check ];then

	for i in $TESTS_DIR/check/*.conf
	do
		base=${i/*\/}
		base=${base/.conf/}
		echo "Running checkpoint test $i"
		echo "Running checkpoint test $i" >&2
		
		$GARLI_BIN $i $GARLI_ARGS
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi

		$GARLI_BIN $TESTS_DIR/restart/$base.conf $GARLI_ARGS
		if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
			exit 1
		fi

		#NEXUSvalidator gives a warning every time it reads a tree file
		#without a taxa block.  So, shut it up initially and then if it
		#fails let it output whatever error
		$NEXUS_VAL ch.$base*.tre 2> /dev/null
		if [ $? -eq 0 ];then
			echo "TREEFILES PASS"
		else
			$NEXUS_VAL ch.$base*.tre
			if [[ ! $? -eq 0 && ! -n "$NO_EXIT_ON_ERR" ]];then
				exit 1
			fi
		fi
	done
else
	echo "No output tests found ..."
fi
