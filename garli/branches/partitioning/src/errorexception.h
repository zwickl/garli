// GARLI version 0.96b8 source code
// Copyright 2005-2008 Derrick J. Zwickl
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
#include <ostream>
#include "outputman.h"

using namespace std;

extern OutputManager outman;

class ErrorException{
	
	public:
	char message[400];
	ErrorException(const char *fmt, ...){
		va_list vl;
		va_start(vl, fmt);
		vsprintf(message, fmt, vl);
		va_end(vl);
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

