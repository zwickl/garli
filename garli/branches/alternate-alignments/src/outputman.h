// GARLI version 2.1 source code
// Copyright 2005-2014 Derrick J. Zwickl
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

#ifndef OUTPUTMANAGER
#define OUTPUTMANAGER

#include <stdarg.h>
#include <ostream>
#include <fstream>
#include <iostream>
//#include <xiosbase>

using namespace std;

class fmtflags;

#define BUFFER_LENGTH 500

//Default behavior here is to output all UserMessage calls to cout, and if SetLogFile or SetLogFileForAppend are called
//then the exact same messages are output to the filenames specified in those calls.
//flow is 
//UserMessage{
//	construct message from printf style spec
//	call Print or PrintNoCR
//		use << to output to streams	
//
//SetNoOutput sets the noOutput flag, which nixes output to defaultOut, but NOT logOut
//
//debugOut is by default set to the same as logOut, unless SetDebugStream is called

class OutputManager{
	char message[BUFFER_LENGTH+1];
	ostream *defaultOut;
	ofstream logOut;
	ostream *debugOut;
	bool noOutput;
	bool log;

	public:
		OutputManager(){
			noOutput = false;
			log = false;
			defaultOut = &cout;
			debugOut = NULL;
			}
			
		~OutputManager(){
			if(log) 
				logOut.close();
			}

		bool IsLogSet(){
			return log;
			}

		void SetOutputStream(ostream &out){
			defaultOut = &out;
			}

		void SetDebugStream(ostream &out){
			debugOut = &out;
			}

		void SetLogFile(const char *logname){
			log = true;
			if(logOut.is_open()){
				logOut.close();
				logOut.clear();
				}
			logOut.open(logname);
			debugOut = &logOut;
			}
		
		void SetLogFileForAppend(const char *logname){
			log = true;
			if(logOut.is_open()){
				logOut.close();
				logOut.clear();
				}
			logOut.open(logname, ios::app);
			debugOut = &logOut;
			}

		void CloseLogFile(){
			logOut.close();
			}

		ofstream *GetLogStream(){
			if(log) 
				return &logOut;
			else 
				return NULL;
			}

		ostream *GetOutputStream(){
			if(noOutput) 
				return NULL;
			else 
				return defaultOut;
			}

		void SetNoOutput(bool o){
			noOutput = o;	
			}

		void precision(const int p){
			defaultOut->precision(p);
			if(log) 
				logOut.precision(p);
			}

		void setf(const std::ios_base::fmtflags &flags){
			defaultOut->setf(flags);
			if(log) 
				logOut.setf(flags);				
			}

		void unsetf(const std::ios_base::fmtflags &flags){
			defaultOut->unsetf(flags);
			if(log) 
				logOut.unsetf(flags);				
			}

		void UserMessage(const char *fmt, ...){
			va_list vl;
			int len;

			va_start(vl, fmt); len = vsnprintf(message, BUFFER_LENGTH, fmt, vl); va_end(vl);
			
			if(len > -1 && len < BUFFER_LENGTH){
				Print(*defaultOut);
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage = NULL;
				if(len > -1){//on unix systems vsnprintf returns the required length.  There is some
					//some ambiguity about whether it includes the null termination or not, but
					//the number passed to vsnprintf should definitely include it.
					longmessage = new char[len+2];
					va_start(vl, fmt); vsnprintf(longmessage, len+1, fmt, vl); va_end(vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					len = BUFFER_LENGTH * 2;
					longmessage = new char[len+1];
					va_start(vl, fmt);
					while(vsnprintf(longmessage, len, fmt, vl) < 0){
						delete []longmessage;
						len *= 2;
						longmessage = new char[len+1];
						va_end(vl);
						va_start(vl, fmt);
						}
					va_end(vl);
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					if(longmessage) 
						delete []longmessage;
					return;
#endif
					}
				Print(*defaultOut, longmessage);
				if(longmessage) 
					delete []longmessage;
				}
			}

		void UserMessageNoCR(const char *fmt, ...){
			va_list vl;
			int len;

			va_start(vl, fmt); len = vsnprintf(message, BUFFER_LENGTH, fmt, vl); va_end(vl);

			if(len > -1 && len < BUFFER_LENGTH){
				PrintNoCR(*defaultOut);
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage = NULL;
				if(len > -1){//on unix systems vsnprintf returns the required length.  There is some
					//some ambiguity about whether it includes the null termination or not, but
					//the number passed to vsnprintf should definitely include it.
					longmessage = new char[len+2];
					va_start(vl, fmt); vsnprintf(longmessage, len+1, fmt, vl); va_end(vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					len = BUFFER_LENGTH * 2;
					longmessage = new char[len+1];
					va_start(vl, fmt);
					while(vsnprintf(longmessage, len, fmt, vl) < 0){
						delete []longmessage;
						len *= 2;
						longmessage = new char[len+1];
						va_end(vl);
						va_start(vl, fmt);
						}
					va_end(vl);
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					if(longmessage) delete []longmessage;
					return;
#endif
					}
				PrintNoCR(*defaultOut, longmessage);
				if(longmessage) 
					delete []longmessage;
				}
			}

		void DebugMessage(const char *fmt, ...){
#ifdef DEBUG_MESSAGES
			if(!debugOut) 
				return;

			va_list vl;
			int len;
			
			va_start(vl, fmt); len = vsnprintf(message, BUFFER_LENGTH, fmt, vl); va_end(vl);

			if(len > -1 && len < BUFFER_LENGTH){
				//Print(*debugOut);
				*debugOut << message << endl;
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage = NULL;
				if(len > -1){//on unix systems vsnprintf returns the required length.  There is some
					//some ambiguity about whether it includes the null termination or not, but
					//the number passed to vsnprintf should definitely include it.
					longmessage = new char[len+2];
					va_start(vl, fmt); vsnprintf(longmessage, len+1, fmt, vl); va_end(vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					len = BUFFER_LENGTH * 2;
					longmessage = new char[len+1];
					va_start(vl, fmt);
					while(vsnprintf(longmessage, len, fmt, vl) < 0){
						delete []longmessage;
						len *= 2;
						longmessage = new char[len+1];
						va_end(vl);
						va_start(vl, fmt);
						}
					va_end(vl);
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					if(longmessage) 
						delete []longmessage;
					return;
#endif
					}
				//Print(*debugOut, longmessage);
				*debugOut << longmessage << endl;
				if(longmessage) 
					delete []longmessage;
				}
#endif
			}

		void DebugMessageNoCR(const char *fmt, ...){
#ifdef DEBUG_MESSAGES
			if(!debugOut) 
				return;

			va_list vl;
			int len;

			va_start(vl, fmt); len = vsnprintf(message, BUFFER_LENGTH, fmt, vl); va_end(vl);

			if(len > -1 && len < BUFFER_LENGTH){
				//PrintNoCR(*debugOut);
				*debugOut << message;
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage = NULL;
				if(len > -1){//on unix systems vsnprintf returns the required length.  There is some
					//some ambiguity about whether it includes the null termination or not, but
					//the number passed to vsnprintf should definitely include it.
					longmessage = new char[len+2];
					va_start(vl, fmt); vsnprintf(longmessage, len+1, fmt, vl); va_end(vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					len = BUFFER_LENGTH * 2;
					longmessage = new char[len+1];
					va_start(vl, fmt);
					while(vsnprintf(longmessage, len, fmt, vl) < 0){
						delete []longmessage;
						len *= 2;
						longmessage = new char[len+1];
						va_end(vl);
						va_start(vl, fmt);
						}
					va_end(vl);
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					if(longmessage) 
						delete []longmessage;
					return;
#endif
					}
				//PrintNoCR(*debugOut, longmessage);
				*debugOut << longmessage;
				if(longmessage) 
					delete []longmessage;
				}
#endif
			}

		void UserMessage(const string &mess){
			Print(*defaultOut, mess);
			}

		void UserMessageNoCR(const string &mess){
			PrintNoCR(*defaultOut, mess);
			}

		void flush(){
			if(!noOutput) 
				defaultOut->flush();
			if(log) 
				logOut.flush();	
			}

		void Print(ostream &out){
			if(!noOutput) 
				out << message << endl;
			if(log)
				logOut << message << endl;
			}

		void PrintNoCR(ostream &out){
			if(!noOutput) 
				out << message;
			if(log) 
				logOut << message;
			}

		void Print(ostream &out, const string &mess){
			if(!noOutput) 
				out << mess << endl;
			if(log) 
				logOut << mess << endl;
			}

		void PrintNoCR(ostream &out, const string &mess){
			if(!noOutput) 
				out << mess;
			if(log) 
				logOut << mess;
			}
	};
	
#endif

