#ifndef ERROREXCEPTION
#define ERROREXCEPTION

#include <stdarg.h>
#include <ostream>
#include "outputman.h"

using namespace std;

extern OutputManager outman;

class ErrorException{
	
	public:
	char message[200];
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
	
#endif

