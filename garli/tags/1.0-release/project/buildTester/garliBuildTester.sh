#!/bin/sh
set -x
build_file="cfgcommand.sh"
num_failures=0
num_attempted=0
failed_variants=""

rm -f garliBuildMasterTesterSuccesses.txt
rm -f garliBuildMasterTesterFailures.txt

parent_dir=`pwd`


for s in `ls | grep ^build`
do
	if test -d "$s"
	then
		if test -f "$s/${build_file}"
		then
			result=0
			num_attempted=`expr ${num_attempted} + 1`
			####################################################################
			# build steps
			####################################################################
			cd "$s"
			echo "Trying build in ${s}"
			make maintainer-clean > /dev/null 2>/dev/null
			rm -f garliBuildTesterLog.txt
			rm -f garliBuildTesterErrorLog.txt
			for cmd in "sh ${build_file}" "make" "make check" "make install" "make installcheck"
			do 
				if ! ${cmd} >> garliBuildTesterLog.txt 2>>garliBuildTesterErrorLog.txt
				then
					result=1
					fail_msg="${cmd} failed"
					break
				fi
			done
			####################################################################
			# go back to the parent dir and write results
			cd -
			if test $result -eq 1
			then
				echo "${parent_dir}/$s" >> garliBuildMasterTesterFailures.txt
				num_failures=`expr ${num_failures} + 1`
			else
				echo "${parent_dir}/$s" >> garliBuildMasterTesterSuccesses.txt
			fi
		else
			echo "$s does not appear to be a Garli build dir ${build_file} not found"
		fi
	fi
done

echo "${num_failures} out of ${num_attempted} Failures" > garliBuildMasterTesterLog.txt

if test $num_failures -gt 0
then
	exit 1
fi
