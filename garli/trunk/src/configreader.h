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

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <stdio.h>
#include <string>
#include <map>

//using std::string;
//using std::map;

typedef map<string, string> Options;
typedef map<string, Options> Sections;

class ConfigReader	{
public:
	ConfigReader();
	ConfigReader(const char*);
	~ConfigReader();

	int Load(const char*);
	int Save(const char*);

	int AddSection(const char*);
	int RemoveSection(const char*);
	int SetSection(const char*);

	int SetOption(const char*, const char*);
	int RemoveOption(const char*);

	int GetStringOption(const char*, string&, bool optional=false);
	int GetBoolOption(const char*, bool&, bool optional=false);
	int GetIntOption(const char*, int&, bool optional=false);
	int GetIntNonZeroOption(const char*, int&, bool optional=false);
	int GetIntRangeOption(const char*, int&, int&);
	int GetUnsignedOption(const char* option, unsigned& val, bool optional=false);
	int GetUnsignedNonZeroOption(const char* option, unsigned& val, bool optional=false);
	int GetFloatOption(const char*, float&);
	int GetFloatRangeOption(const char*, float&, float&);
	int GetDoubleOption(const char*, FLOAT_TYPE&, bool optional=false);
	int GetDoubleRangeOption(const char*, FLOAT_TYPE&, FLOAT_TYPE&);
	int GetPositiveDoubleOption(const char*, FLOAT_TYPE&, bool optional=false);
	int GetPositiveNonZeroDoubleOption(const char* option, FLOAT_TYPE& val, bool optional=false);

	Sections::const_iterator BeginSection() const	{
		return sections.begin();
	}
	Sections::const_iterator EndSection() const	{
		return sections.end();
	}

private:
	static int UNKNOWN;
	static int SECTION;
	static int OPTION;

private:

	int ReadSectionOrOption(FILE* file, string& name, string& val);
	int ReadLine(FILE* file, string& line);

	void TrimWhiteSpace(string& str);

	Sections sections;
	string cur_section;
};


#endif
