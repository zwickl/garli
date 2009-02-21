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

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include <math.h>

using namespace std;

#include "defs.h"
#include "configreader.h"
#include "errorexception.h"
#include "funcs.h"

int ConfigReader::UNKNOWN=0;
int ConfigReader::SECTION=1;
int ConfigReader::OPTION=2; 


bool ParseLineIntoNameValue(const std::string & line, std::string & name, std::string & val) {
	const size_t len = line.length();
	const size_t index = line.find('=');
	if (index < len) {
		val = line.substr(index+1);
		TrimWhiteSpace(val);
		name = line.substr(0, index);
		TrimWhiteSpace(name);
		return true;
	}
	return false; 
}


ConfigReader::ConfigReader()	{

}

ConfigReader::ConfigReader(const char* filename)	{
	Load(filename);
}

ConfigReader::~ConfigReader()	{

}

int ConfigReader::Load(const char* filename)	{
	FILE* file;

	// clear all currently loaded sections
	sections.clear();

#ifndef BOINC
	file = fopen(filename, "r");
#else
	char input_path[512];
    boinc_resolve_filename(filename, input_path, sizeof(input_path));
    file = boinc_fopen(input_path, "r");
#endif

	if (file == NULL) throw ErrorException("could not open file \"%s\".", filename);

	int type;
	string sectionName, name, val;
	Options options;
	bool first = true;

	while (!feof(file))	{

		type = ReadSectionOrOption(file, name, val);

		if (type == SECTION)	{
			if (!first)	{
				sections.insert(pair<string, Options>(sectionName, options));
				options.clear();
			}
			else
				first = false;
			sectionName = name;
		}
		else if (type == OPTION)	{
			options.insert(pair<string, string>(name, val));
		}

	}

	// insert the last section
	sections.insert(pair<string, Options>(sectionName, options));

	fclose(file);

	return 0;
}

/*****************************************************************************************
** Save() **
*****************************************************************************************/
int ConfigReader::Save(const char* filename)	{
	FILE* file = fopen(filename, "w");

	if (file == NULL)	{
		printf("Error opening file \"%s\" for writing.\n", filename);
		return -1;
	}

	map<string, Options>::iterator sit = sections.begin();
	map<string, string>::iterator oit;

	while (sit != sections.end())	{

		// write the section
		fprintf(file, "[%s]\n", sit->first.c_str());

		// write the options
		oit = sit->second.begin();
		while (oit != sit->second.end())	{
			fprintf(file, "%s = %s\n", oit->first.c_str(), oit->second.c_str());
			++oit;
		}

		fprintf(file, "\n");
		++sit;

	}

	fclose(file);

	return 0;
}

/****************************************************************************************/
/*** AddSection() ***/
/****************************************************************************************/
int ConfigReader::AddSection(const char* _name)	{
	int rv;
	string name;
	Sections::iterator it;

	name = _name;
	TrimWhiteSpace(name);
	it = sections.find(name);
	if (it == sections.end())
		rv = 0;
	else
		rv = 1;
	sections.insert(pair<string, Options>(name, Options()));
	return rv;
}

/****************************************************************************************/
/*** RemoveSection() ***/
/****************************************************************************************/
int ConfigReader::RemoveSection(const char* _name)	{
	int rv;
	string name;
	Sections::iterator it;

	name = _name;
	TrimWhiteSpace(name);
	it = sections.find(name);
	if (it == sections.end())	// section doesn't exist, bomb out
		rv = -1;
	else	{
		sections.erase(it);
		rv = 0;
	}
	return rv;
}

/****************************************************************************************/
/*** SetSection() ***/
/****************************************************************************************/
int ConfigReader::SetSection(const char* name)	{
	cur_section = name;
	TrimWhiteSpace(cur_section);
	Sections::iterator sit = sections.find(cur_section);
	if (sit == sections.end())
		return -1;
	return 0;
}

/****************************************************************************************/
/*** SetOption() ***/
/****************************************************************************************/
int ConfigReader::SetOption(const char* _option, const char* _val)	{
	int rv;
	string option, val;
	Sections::iterator sit;
	Options::iterator oit;

	sit = sections.find(cur_section);
	if (sit == sections.end())	// section doesn't exist...bomb out
		rv = -1;
	else	{
		option = _option;
		val = _val;
		TrimWhiteSpace(option);
		TrimWhiteSpace(val);
		oit = sit->second.find(option);
		if (oit == sit->second.end())	{	// option doesn't exist, create it
			rv = 0;
		}
		else	{	// option exists, overwrite
			sit->second.erase(oit);
			rv = 1;
		}
		sit->second.insert(pair<string, string>(option, val));
	}

	return rv;
}

/****************************************************************************************/
/*** RemoveOption() ***/
/****************************************************************************************/
int ConfigReader::RemoveOption(const char* _option)	{
	int rv;
	string option;
	Sections::iterator sit;
	Options::iterator oit;

	sit = sections.find(cur_section);
	if (sit == sections.end())	// section doesn't exist, bomb out
		rv = -2;
	else	{
		option = _option;
		TrimWhiteSpace(option);
		oit = sit->second.find(option);
		if (oit == sit->second.end())	{	// option doesn't exist, bomb out
			rv = -1;
		}
		else	{	// option exists, remove it
			sit->second.erase(oit);
			rv = 0;
		}
	}

	return rv;
}

/****************************************************************************************/
/*** GetStringOption() ***/
/****************************************************************************************/
int ConfigReader::GetStringOption(const char* _option, string& val, bool optional /*=false*/)	{
	int rv;
	string option;
	Sections::iterator sit;
	Options::iterator oit;

	sit = sections.find(cur_section);
	if (sit == sections.end())	// section doesn't exist, bomb out
		rv = -2;
	else	{
		option = _option;
		TrimWhiteSpace(option);
		oit = sit->second.find(option);
		if (oit == sit->second.end())	{	// option doesn't exist, bomb out
			rv = -1;
			if(!optional) throw ErrorException("could not find string configuration entry \"%s\"", option.c_str());
		}
		else	{	// option exists, get the value
			val = oit->second;
			rv = 0;
		}
	}

	return rv;
}


bool ParseStringAsBool(const std::string & str, const char * optName) {
	if (str == "true")
		return true;
	if (str == "false")
		return false;
	if(isdigit(str[0]) != 0)
		return (atoi(str.c_str()) != 0);
	throw ErrorException("expecting boolean (0 or 1) for entry \"%s\", found %s", optName, str.c_str());
}

int  ParseStringAsInt(const std::string & str, const char * optName) {
	long dummy;
	if (!NxsString::to_long(str.c_str(), &dummy))
		throw ErrorException("entry for option \"%s\" (%s) is not an integer" , optName, str.c_str());
	if(dummy > long(INT_MAX-1))
		throw ErrorException("entry for option \"%s\" (%s) is greater than its max (%u)" , optName, str.c_str(), (INT_MAX-1));
	int val = (int) dummy;
	return val;
}

FLOAT_TYPE ParseStringAsDouble(const std::string & str, const char * optName) {
	double dummy;
	if (!NxsString::to_double(str.c_str(), &dummy))
		throw ErrorException("entry for option \"%s\" (%s) is not an number" , optName, str.c_str());
	return FLOAT_TYPE(dummy);
}

unsigned ParseStringAsUnsigned(const std::string & str, const char * optName) {
	int i = ParseStringAsInt(str, optName);
	if (i < 0)
		ErrorException("entry for option \"%s\" must be >=0", optName);
	return unsigned (i);
}


/****************************************************************************************/
/*** GetBoolOption() ***/
/****************************************************************************************/
int ConfigReader::GetBoolOption(const char* option, bool& val, bool optional /*=false*/)	{
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
 		val = ParseStringAsBool(str, option);
		return 0;
	}
	if (!optional)
		throw ErrorException("could not find boolean configuration entry \"%s\"", option);
	return -1;
}

/****************************************************************************************/
/*** GetIntOption() ***/
/****************************************************************************************/
int ConfigReader::GetIntOption(const char* option, int& val, bool optional /*=false*/)	{
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = ParseStringAsInt(str, option);
		return 0;
	}
	if(!optional)
		throw ErrorException("could not find integer configuration entry \"%s\"", option);
	return -1;
}

/****************************************************************************************/
/*** GetIntNonZeroOption() ***/
/****************************************************************************************/
int ConfigReader::GetIntNonZeroOption(const char* option, int& val, bool optional /*=false*/)	{
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		int dummy = ParseStringAsInt(str, option);
		if (dummy == 0)
			throw ErrorException("entry for option \"%s\" cannot be zero", option);
		val = dummy;
		return 0;
	}
	if(!optional)
		throw ErrorException("could not find integer configuration entry \"%s\"", option);
	return -1;
}

/****************************************************************************************/
/*** GetUnsignedOption() ***/
/****************************************************************************************/
int ConfigReader::GetUnsignedOption(const char* option, unsigned& val, bool optional /*=false*/)	{
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = ParseStringAsUnsigned(str, option);
		return 0;
	}
	if (!optional)
		throw ErrorException("could not find unsigned integer configuration entry \"%s\"", option);
	return -1;
}

/****************************************************************************************/
/*** GetUnsignedOption() ***/
/****************************************************************************************/
int ConfigReader::GetUnsignedNonZeroOption(const char* option, unsigned& val, bool optional /*=false*/)	{
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		unsigned u = ParseStringAsUnsigned(str, option);
		if (u == 0)
			throw ErrorException("entry for option \"%s\" must be >0", option);
		val = u;
		return 0;
	}
	if (!optional)
		throw ErrorException("could not find unsigned integer configuration entry \"%s\"", option);
	return -1;
}


/****************************************************************************************/
/*** GetIntRangeOption() ***/
/****************************************************************************************/
int ConfigReader::GetIntRangeOption(const char* option, int& val1, int& val2)	{
	int rv;
	string str;
	if (GetStringOption(option, str) == 0)	{  // option exists

		// split up the string
		int len = (int)str.length();
		int i = (int)str.find(' ', 0);
		if (i < 0)
			rv = -1;
		else	{
			val1 = atoi(str.substr(0, i).c_str());
			val2 = atoi(str.substr(i+1, len-i).c_str());
			rv = 0;
		}
	}
	else{
		rv = -1;
		throw ErrorException("could not find integer range configuration entry \"%s\"", option);
		}

	return rv;
}

/****************************************************************************************/
/*** GetFloatOption() ***/
/****************************************************************************************/
int ConfigReader::GetFloatOption(const char* option, float& val)	{
	int rv;
	string str;
	if (GetStringOption(option, str) == 0)	{  // option exists
		val = (float)atof(str.c_str());
		rv = 0;
	}
	else
		rv = -1;

	return rv;
}

/****************************************************************************************/
/*** GetFloatRangeOption() ***/
/****************************************************************************************/
int ConfigReader::GetFloatRangeOption(const char* option, float& val1, float& val2)	{
	int rv;
	string str;
	if (GetStringOption(option, str) == 0)	{  // option exists

		// split up the string
		int len = (int)str.length();
		int i = (int)str.find(' ', 0);
		if (i < 0)
			rv = -1;
		else	{
			val1 = (float)atof(str.substr(0, i).c_str());
			val2 = (float)atof(str.substr(i+1, len-i).c_str());
			rv = 0;
		}
	}
	else
		rv = -1;

	return rv;
}

/****************************************************************************************/
/*** GetDoubleOption() ***/
/****************************************************************************************/
int ConfigReader::GetDoubleOption(const char* option, FLOAT_TYPE& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = (FLOAT_TYPE) atof(str.c_str());
		rv = 0;
	}
	else{
		rv = -1;
		if(!optional) throw ErrorException("error: could not find float configuration entry \"%s\"", option);
		}

	return rv;
}

/****************************************************************************************/
/*** GetPositiveDoubleOption() ***/
/****************************************************************************************/
//this is just a version of GetDoubleOption that checks that the value is non-negative
int ConfigReader::GetPositiveDoubleOption(const char* option, FLOAT_TYPE& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = (FLOAT_TYPE) atof(str.c_str());
		if(val < 0.0) throw ErrorException("configuration entry \"%s\" cannot be negative", option);
		rv = 0;
	}
	else{
		rv = -1;
		if(!optional) throw ErrorException("could not find float configuration entry \"%s\"", option);
		}

	return rv;
}

/****************************************************************************************/
/*** GetPositiveNonZeroDoubleOption() ***/
/****************************************************************************************/
//this is just a version of GetDoubleOption that checks that the value is non-negative, and that it is not
//zero.  atof returns zero when it encounters an error, which is very annoying behavior.  When the entry must
//be nonzero, at least we can check for that
int ConfigReader::GetPositiveNonZeroDoubleOption(const char* option, FLOAT_TYPE& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = (FLOAT_TYPE) atof(str.c_str());
		if(val == ZERO_POINT_ZERO) throw ErrorException("configuration entry \"%s\" cannot be zero (possible problems reading this entry)", option);
		if(val < 0.0) throw ErrorException("configuration entry \"%s\" cannot be negative", option);
		rv = 0;
	}
	else{
		rv = -1;
		if(!optional) throw ErrorException("could not find float configuration entry \"%s\"", option);
		}

	return rv;
}

/****************************************************************************************/
/*** GetDoubleRangeOption() ***/
/****************************************************************************************/
int ConfigReader::GetDoubleRangeOption(const char* option, FLOAT_TYPE& val1, FLOAT_TYPE& val2)	{
	int rv;
	string str;
	if (GetStringOption(option, str) == 0)	{  // option exists

		// split up the string
		int len = (int)str.length();
		int i = (int)str.find(' ', 0);
		if (i < 0)
			rv = -1;
		else	{
			val1 = (FLOAT_TYPE) atof(str.substr(0, i).c_str());
			val2 = (FLOAT_TYPE) atof(str.substr(i+1, len-i).c_str());
			rv = 0;
		}
	}
	else{
		rv = -1;
		throw ErrorException("could not find float range configuration entry \"%s\"", option);
		}

	return rv;
}

/*****************************************************************************************
** PRIVATE METHODS ***********************************************************************
*****************************************************************************************/

int ConfigReader::ReadSectionOrOption(FILE* file, string& name, string& val)	{

	string line;
	size_t index;
	size_t len;
	int type = UNKNOWN;

	do	{
		len = ReadLine(file, line);
		if (line.find('=') < len)
			type = OPTION;
		else if (line.find('[') < len && line.find(']') < len)
			type = SECTION;
	}
	while (type == UNKNOWN && !feof(file));

	if (type == SECTION)	{
		line.erase(line.find('['), 1);
		line.erase(line.find(']'), 1);
		name = line;
	}
	else if (type == OPTION)	{
		index = line.find('=');
		val = line.substr(index+1);
		name = line.substr(0, index);
	}

	TrimWhiteSpace(name);
	if (type == OPTION)
		TrimWhiteSpace(val);

	return type;
}

int ConfigReader::ReadLine(FILE* file, string& line)	{
	char ch;

	line = "";

	fread(&ch, sizeof(char), 1, file);
	while (ch != '\n' && ch != '\r' && !feof(file))	{
		line += ch;
		fread(&ch, sizeof(char), 1, file);
	}

	return (int)line.length();
}

void TrimWhiteSpace(string& str)	{
	int index;

	if (str.length() == 0)
		return;

	index = (int)str.find(' ', 0);
	while (index != -1 && index < (int)str.length())	{
		while (index < (int)str.length()-1 && str[index+1] == ' ')
			str.erase(index+1, 1);
		index = (int)str.find(' ', index+1);
	}

	if (str.find(' ', 0) == 0)
		str.erase(0, 1);
 	if ( (str.length() > 0) && (str.find(' ', str.length()-1) == str.length()-1) )
 		str.erase(str.length()-1, 1);

}
