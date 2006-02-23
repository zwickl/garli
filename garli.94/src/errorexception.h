#ifndef ERROREXCEPTION
#define ERROREXCEPTION

#include <stdarg.h>

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
		out << message << endl;
		}
	};
	
#endif

