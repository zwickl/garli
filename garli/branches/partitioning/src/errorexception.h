// GARLI version 1.00 source code
// Copyright 2005-2010 Derrick J. Zwickl
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

#ifndef ERROREXCEPTION
#define ERROREXCEPTION

#include <stdarg.h>
#include <cassert>
#include <cstring>
#include <ostream>
#include "outputman.h"

using namespace std;

extern OutputManager outman;

#define BUFFER_LENGTH 500

class ErrorException{
	
	public:
	char *message;
	int messlen;
	//char message[5000];
	//char message[400];
	ErrorException(){
		message = NULL;
		}
/*
	ErrorException(const char *fmt, ...){
		message = new char[500];
		va_list vl;
		va_start(vl, fmt);
		vsprintf(message, fmt, vl);
		assert(strlen(message) < 5000);
		va_end(vl);
		}
*/
	ErrorException(const ErrorException &other){
		messlen = strlen(other.message);
		message = new char[messlen + 1];
		strcpy(message, other.message);
		}

	ErrorException(const char *fmt, ...){
		messlen = BUFFER_LENGTH;
		message = new char[messlen];

		va_list vl;
		va_start(vl, fmt);
		int len = vsnprintf(message, messlen, fmt, vl);
		va_end(vl);

		if((len > -1 && len < messlen) == false){//default buffer is not long enough or there was an error
			delete []message;
			message = NULL;
			//char *longmessage = NULL;
			if(len > -1){//on unix systems vsnprintf returns the required length.  There is some
				//some ambiguity about whether it includes the null termination or not, but
				//the number passed to vsnprintf should definitely include it.
					
				message = new char[len+2];
				va_start(vl, fmt);
				vsnprintf(message, len+1, fmt, vl);
				va_end(vl);
				}
			else{
#if defined(_MSC_VER)
				//on windows a negative value means that the length wasn't engough
				messlen = BUFFER_LENGTH * 2;
				message = new char[messlen+1];
				va_start(vl, fmt);
				while(vsnprintf(message, messlen, fmt, vl) < 0){
					delete []message;
					messlen *= 2;
					message = new char[messlen+1];
					va_end(vl);
					va_start(vl, fmt);
					}
				va_end(vl);
#else
				//otherwise negative means a formatting error
				sprintf(message, "(problem formatting some program output...)");
				return;
#endif
				}
			//message = longmessage;
			}
		}

	~ErrorException(){
		delete []message;
		}

	void Print(ostream &out){
		outman.UserMessage("ERROR!: %s\n\n", message);
		//out << "ERROR!: " << message << endl << endl;
		}

	void Print(FILE *out){
		fprintf(out, "ERROR!: %s\n\n", message);
		}
	};
	
class UnscoreableException{
public:
	UnscoreableException(){};
	};

#endif

