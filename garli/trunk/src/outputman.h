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

class OutputManager{
	char message[BUFFER_LENGTH+1];
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

		bool IsLogSet(){
			return (log == true);
			}

		void SetOutputStream(ostream &out){
			defaultOut=&out;
			}

		void SetLogFile(const char *logname){
			log=true;
			if(logOut.is_open()){
				logOut.close();
				logOut.clear();
				}
			logOut.open(logname);
			}

		ofstream *GetLogStream(){
			if(log == true) return &logOut;
			else return NULL;
			}

		ostream *GetOutputStream(){
			if(noOutput) return NULL;
			else return defaultOut;
			}

		void SetLogFileForAppend(const char *logname){
			log=true;
			if(logOut.is_open()){
				logOut.close();
				logOut.clear();
				}
			logOut.open(logname, ios::app);
			}

		void CloseLogFile(){
			logOut.close();
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
			
			int len = vsnprintf(message, BUFFER_LENGTH, fmt, vl);

			if(len > -1 && len < BUFFER_LENGTH){
				Print(*defaultOut);
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage;
				if(len > -1){//on unix systems vsnprintf returns the required length
					longmessage = new char[len+1];
					vsnprintf(longmessage, len, fmt, vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					int len2 = BUFFER_LENGTH * 2;
					longmessage = new char[len2+1];
					while(vsnprintf(longmessage, len2, fmt, vl) < 0){
						delete []longmessage;
						len2 *= 2;
						longmessage = new char[len2+1];
						}
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					va_end(vl);
					return;
#endif
					}
				Print(*defaultOut, longmessage);
				delete []longmessage;
				}
			va_end(vl);
			}

		void UserMessageNoCR(const char *fmt, ...){
			va_list vl;
			va_start(vl, fmt);
			
			int len = vsnprintf(message, BUFFER_LENGTH, fmt, vl);

			if(len > -1 && len < BUFFER_LENGTH){
				PrintNoCR(*defaultOut);
				}
			else{//default buffer is not long enough or there was an error
				char *longmessage;
				if(len > -1){//on unix systems vsnprintf returns the required length
					longmessage = new char[len+1];
					vsnprintf(longmessage, len, fmt, vl);
					}
				else{
#if defined(_MSC_VER)
					//on windows a negative value means that the length wasn't engough
					int len2 = BUFFER_LENGTH * 2;
					longmessage = new char[len2+1];
					while(vsnprintf(longmessage, len2, fmt, vl) < 0){
						delete []longmessage;
						len2 *= 2;
						longmessage = new char[len2+1];
						}
#else
					//otherwise negative means a formatting error
					Print(*defaultOut, "(problem formatting some program output...)");
					va_end(vl);
					return;
#endif
					}
				PrintNoCR(*defaultOut, longmessage);
				delete []longmessage;
				}
			va_end(vl);
			}

/*		void UserMessageNoCR(const char *fmt, ...){
			va_list vl;
			va_start(vl, fmt);
			vsprintf(message, fmt, vl);
			va_end(vl);
			PrintNoCR(*defaultOut);
			}
*/
		void UserMessage(const string &mess){
			Print(*defaultOut, mess);
			}

		void UserMessageNoCR(const string &mess){
			PrintNoCR(*defaultOut, mess);
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

		void Print(ostream &out, const string &mess){
			if(noOutput == false) out << mess << endl;
			if(log==true) logOut << mess << endl;
			}

		void PrintNoCR(ostream &out, const string &mess){
			if(noOutput == false) out << mess;
			if(log==true) logOut << mess;
			}
	};
	
#endif

