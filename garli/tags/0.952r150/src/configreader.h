// GARLI version 0.951 source code
// Copyright  2005-2006 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	National Evolutionary Synthesis Center
//	2024 W. Main Street, Suite A200
//	Durham, NC 27705
//  email: zwickl@nescent.org
//


#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <stdio.h>
#include <string>
#include <map>

using std::string;
using std::map;

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
	int GetIntRangeOption(const char*, int&, int&);
	int GetUnsignedOption(const char* option, unsigned& val, bool optional=false);
	int GetFloatOption(const char*, float&);
	int GetFloatRangeOption(const char*, float&, float&);
	int GetDoubleOption(const char*, double&, bool optional=false);
	int GetDoubleRangeOption(const char*, double&, double&);
	int GetPositiveDoubleOption(const char*, double&, bool optional=false);

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
