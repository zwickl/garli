// GARLI version 2.0 source code
// Copyright 2005-2011 Derrick J. Zwickl
// email: garli.support@gmail.com
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

#ifndef GAML_UTIL_HPP
#define GAML_UTIL_HPP
//some code from Mark Holder for allocating flattened matrices, and other misc stuff

#include <stdlib.h>
#include <cassert>
#include <ostream>
#include <iostream>
#include <iomanip>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <sys/time.h>
#endif

using namespace std;

#include "errorexception.h"

#define DBL_ALIGN 32

template<typename T> T ***New3DArray(unsigned f , unsigned s , unsigned t);
template<typename T> T **New2DArray(unsigned f , unsigned s);
template<typename T> void Delete3DArray	(T ***temp);
template<typename T> void Delete2DArray	(T **temp);

//aligned versions
template<typename T> T ***New3DAlignedArray(unsigned f , unsigned s , unsigned t, unsigned a);
template<typename T> T **New2DAlignedArray(unsigned f , unsigned s, unsigned a);
template<typename T> void Delete3DAlignedArray	(T ***temp);
template<typename T> void Delete2DAlignedArray	(T **temp);

template<typename T> T *NewAlignedArray(unsigned len, unsigned align ){
#ifdef _MSC_VER
	return (T*) _aligned_malloc(sizeof(T)*len, align);
#endif
}

template<typename T> void DeleteAlignedArray(T *a){
#ifdef _MSC_VER
	_aligned_free(a);
#endif
}

/*--------------------------------------------------------------------------------------------------------------------------
| Allocates a three dimensional array of FLOAT_TYPEs as one contiguous block of memory
| the dimensions are f two dimensional arrays that are s by t.  
| the array is set up so that 
| for(i = 0 ; i < f ; i++)
|	for (j = 0 ; j < s ; j++)
|		for (k = 0 ; k < t; k++)
|			array[i][j][k];
| would be the same order of access as: 
| 
|	T *temp = **array;
|	for (i = 0 ; i < f*s*t ; i++)
|		{
|		*temp++;
|		}
*/
template<typename T> T ***New3DArray(unsigned f , unsigned s , unsigned t)
	{
	assert(f > 0 && s > 0 && t> 0);
	T ***temp;
	try{
		temp = new T **[f];
		*temp = new T *[f * s];
		**temp = new T[f * s * t];
		for (unsigned sIt = 1 ; sIt < s ; sIt++)
			temp[0][sIt] = temp[0][sIt-1] + t ;
		for (unsigned fIt = 1 ; fIt < f ; fIt ++)
			{
			temp[fIt] = temp[fIt -1] +  s ;
			temp[fIt][0] = temp[fIt -1][0] + (s*t);
			for (unsigned sIt = 1 ; sIt < s ; sIt++)
				temp[fIt][sIt] = temp[fIt][sIt-1] + t ;
			}
		}
	catch(std::bad_alloc){
		throw ErrorException("Problem allocating 3D array (%d X %d X %d = %.2f MB). Out of mem?", f, s, t, (f * s * t * sizeof(T)) / (1024.0 * 1024.0));
		}
	return temp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Delete a Three Dimensional Array that has been allocated using New3DArray
*/
template<typename T> void Delete3DArray	(T ***temp)
	{
	assert(temp); //these asserts aren't necessary, but right now I can't think of a case in which they'd fail other than following an allocation error
	assert(*temp);
	assert(**temp);
	if (temp)
		{
		if (*temp)
			{
			if (**temp)
				delete [] **temp;
			delete [] * temp;
			}
		delete [] temp;
		}
	
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Allocates a two dimensional array of FLOAT_TYPEs as one contiguous block of memory
| 	the dimensions are f by s.  
| 	the array is set up so that 
| 	
|	for(i = 0 ; i < f ; i++)
|		for (j = 0 ; j < s ; j++)
|			array[i][j];
| 	
|	would be the same order of access as: 
| 
|  	T *temp = **array;
|	for (i = 0 ; i < f*s*t ; i++)
|		*temp++;
*/
template<typename T> T **New2DArray(unsigned f , unsigned s)
	{
	assert(f > 0 && s > 0);
	T **temp;
	try{
		temp = new T *[f];
		*temp = new T [f * s];
		for (unsigned fIt = 1 ; fIt < f ; fIt ++)
			temp[fIt] = temp[fIt -1] +  s ;
		}
	catch(std::bad_alloc){
		throw ErrorException("Problem allocating 2D array (%d X %d = %.2f MB). Out of mem?", f, s, (f * s * sizeof(T)) / (1024.0 * 1024.0));
		}
	return temp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Delete a 2 Dimensional Array New2DArray
*/
template<typename T> inline void Delete2DArray	(T **temp)
	{
	assert(temp); //these asserts aren't necessary, but right now I can't think of a case in which they'd fail other than following an allocation error
	assert(*temp);
	if (temp)
		{
		if (*temp)
			delete [] * temp;
		delete [] temp;
		}
	}

//aligned version
template<typename T> T ***New3DAlignedArray(unsigned f , unsigned s , unsigned t)
	{
	assert(f > 0 && s > 0 && t> 0);
	T ***temp;
	temp = new T **[f];
	*temp = new T *[f * s];
	**temp = new T[f * s * t];
	**temp = NewAlignedArray<T>(f * s * t, DBL_ALIGN);
	for (unsigned sIt = 1 ; sIt < s ; sIt++)
		temp[0][sIt] = temp[0][sIt-1] + t ;
	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
		{
		temp[fIt] = temp[fIt -1] +  s ;
		temp[fIt][0] = temp[fIt -1][0] + (s*t);
		for (unsigned sIt = 1 ; sIt < s ; sIt++)
			temp[fIt][sIt] = temp[fIt][sIt-1] + t ;
		}
	return temp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Delete a Three Dimensional Array that has been allocated using New3DArray
*/
template<typename T> void Delete3DAlignedArray	(T ***temp)
	{
	assert(temp); //these asserts aren't necessary, but right now I can't think of a case in which they'd fail other than following an allocation error
	assert(*temp);
	assert(**temp);
	if (temp)
		{
		if (*temp)
			{
			if (**temp)
				DeleteAlignedArray(**temp);
			delete [] * temp;
			}
		delete [] temp;
		}
	
	}

template<typename T> T **New2DAlignedArray(unsigned f , unsigned s)
	{
	assert(f > 0 && s > 0);
	T **temp;
	temp = new T *[f];
	*temp = NewAlignedArray<T>(f * s, DBL_ALIGN);
	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
		temp[fIt] = temp[fIt -1] +  s ;
	return temp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Delete a 2 Dimensional Array New2DArray
*/
template<typename T> inline void Delete2DAlignedArray	(T **temp)
	{
	assert(temp); //these asserts aren't necessary, but right now I can't think of a case in which they'd fail other than following an allocation error
	assert(*temp);
	if (temp)
		{
		if (*temp)
			DeleteAlignedArray(*temp);
		delete [] temp;
		}
	}



class Profiler{
#ifdef _MSC_VER
	LONGLONG totalTics;
	int numCalls;
	string name;
	bool inuse;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	LARGE_INTEGER ticsPerSec;
#else
	unsigned totalTics;
	int numCalls;
	string name;
	bool inuse;
	timeval start;
	timeval end;
	struct timezone tz;
	int ticsPerSec;
#endif

public:
	Profiler(string n){
		name = n;
		totalTics = 0;
		numCalls = 0;
		inuse=false;
#ifdef _MSC_VER
		QueryPerformanceFrequency(&ticsPerSec);
#else
		ticsPerSec = 1000000;
#endif
		}
	void Start(){
#ifdef ENABLE_CUSTOM_PROFILER
		if(inuse){
			cout << "Error! Don't use this on recursive functions!" << endl;
			exit(1);
			}
		inuse=true;
		numCalls++;
	#ifdef _MSC_VER
		QueryPerformanceCounter(&start);
	#else
		gettimeofday(&start, &tz);
	#endif
#endif
		}
	void Stop(){
#ifdef ENABLE_CUSTOM_PROFILER
		if(!inuse){
			cout << "Error! Profiler was not started!" << endl;
			exit(1);
			}
	#ifdef _MSC_VER
		QueryPerformanceCounter(&end);
		totalTics += end.QuadPart - start.QuadPart;
	#else
		gettimeofday(&end, &tz);
		totalTics += end.tv_usec - start.tv_usec  + (end.tv_sec - start.tv_sec)*1000000;
		
//		cout <<  end.tv_usec - start.tv_usec  + (end.tv_usec - start.tv_usec)*100000 << endl;
	#endif
		inuse=false;
#endif
		}
	void Report(ostream &out, int progTime){
#ifdef ENABLE_CUSTOM_PROFILER
	#ifdef _MSC_VER
		FLOAT_TYPE seconds = totalTics/(FLOAT_TYPE)ticsPerSec.QuadPart;
	#else
		FLOAT_TYPE seconds = totalTics/(FLOAT_TYPE)ticsPerSec;
	#endif
		out << setw( 10 ) << name.c_str() << "\t" << setw( 10 )<< numCalls << "\t";
		out.precision(4);
		out << setw( 10 ) << seconds << "\t" << setw( 10 ) << seconds/(FLOAT_TYPE)numCalls << "\t" << setw( 10 ) << seconds*100/(FLOAT_TYPE)progTime << endl;

#endif
		}
	};

#endif //



