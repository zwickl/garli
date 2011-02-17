// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
// email: zwickl@nescent.org
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#ifndef STOPWATCH_H
#define STOPWATCH_H


#ifndef UNIX
#include <time.h>
#else
#include <sys/time.h>
#endif

class Stopwatch	{

	public:
	
		Stopwatch()	{
			Restart();
			}
	
		#ifndef UNIX
		void Restart()	{
			time(&start_time);
			}
		void Start()	{
			time(&start_time);
			}		
		int SplitTime()	{
			time(&end_time);
			return (int)(end_time - start_time);		
			}
		//this is for restarting
		void AddPreviousTime(time_t t){
			start_time -= t;
			}
		#else		
		void Restart()	{
			gettimeofday(&start_time, NULL);
			}
		void Start()	{
			gettimeofday(&start_time, NULL);
			}		
		int SplitTime()	{
			gettimeofday(&end_time, NULL);
			return end_time.tv_sec - start_time.tv_sec;
			}
		//this is for restarting
		void AddPreviousTime(int t){
			start_time.tv_sec -= t;
			}
		#endif


	private:
		#ifndef UNIX
		time_t start_time, end_time;
		#else
		timeval start_time, end_time;
		#endif

};

#endif
