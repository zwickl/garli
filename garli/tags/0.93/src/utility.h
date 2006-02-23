#ifndef GAML_UTIL_HPP
#define GAML_UTIL_HPP
//	code from MTH for allocating flattened matrices

#include "memchk.h"
#include <cassert>

template<typename T> T ***New3DArray(unsigned f , unsigned s , unsigned t);
template<typename T> T **New2DArray(unsigned f , unsigned s);
template<typename T> void Delete3DArray	(T ***temp);
template<typename T> void Delete2DArray	(T **temp);

/*--------------------------------------------------------------------------------------------------------------------------
| Allocates a three dimensional array of doubles as one contiguous block of memory
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
| 	Allocates a two dimensional array of doubles as one contiguous block of memory
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
	temp = new T *[f];
	*temp = new T [f * s];
	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
		temp[fIt] = temp[fIt -1] +  s ;
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
#endif //

