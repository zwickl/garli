
#ifndef OUTPUTMANAGER
#define OUTPUTMANAGER

#include <stdarg.h>
#include <ostream>
#include <fstream>
#include <iostream>
//#include <xiosbase>

using namespace std;

class fmtflags;

class OutputManager{
	char message[500];
	ostream *defaultOut;
	ofstream logOut;
	bool noOutput;
	bool log;

	public:
		OutputManager(){
			noOutput=false;
			log=false;
			defaultOut=&cout;
			}
			
		~OutputManager(){
			if(log==true) logOut.close();
			}

		void SetOutputStream(ostream &out){
			defaultOut=&out;
			}

		void SetLogFile(char *logname){
			log=true;
			logOut.open(logname);
			}

		void SetNoOutput(bool o){
			noOutput=o;	
			}

		void precision(const int p){
			defaultOut->precision(p);
			if(log==true) logOut.precision(p);
			}

		void setf(const std::ios_base::fmtflags &flags){
			defaultOut->setf(flags);
			if(log==true) logOut.setf(flags);				
			}

		void unsetf(const std::ios_base::fmtflags &flags){
			defaultOut->unsetf(flags);
			if(log==true) logOut.unsetf(flags);				
			}

		void UserMessage(const char *fmt, ...){
			va_list vl;
			va_start(vl, fmt);
			vsprintf(message, fmt, vl);
			va_end(vl);
			Print(*defaultOut);
			}

		void UserMessageNoCR(const char *fmt, ...){
			va_list vl;
			va_start(vl, fmt);
			vsprintf(message, fmt, vl);
			va_end(vl);
			PrintNoCR(*defaultOut);
			}

		void flush(){
			if(noOutput == false) defaultOut->flush();
			if(log==true) logOut.flush();	
			}

		void Print(ostream &out){
			if(noOutput == false) out << message << endl;
			if(log==true) logOut << message << endl;
			}

		void PrintNoCR(ostream &out){
			if(noOutput == false) out << message;
			if(log==true) logOut << message;
			}
	};
	
#endif

