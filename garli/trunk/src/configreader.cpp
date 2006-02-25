// GARLI version 0.94 source code
// Copyright  2005 by Derrick J. Zwickl
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//  Derrick Zwickl
//	Integrative Biology, UT
//	1 University Station, C0930
//	Austin, TX  78712
//  email: zwickl@mail.utexas.edu
//
//	Note: In 2006  moving to NESCENT (The National
//	Evolutionary Synthesis Center) for a postdoc

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

#include "configreader.h"

using std::pair;

int ConfigReader::UNKNOWN=0;
int ConfigReader::SECTION=1;
int ConfigReader::OPTION=2; 

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

	file = fopen(filename, "r");

	if (file == NULL)	{
		printf("error opening file \"%s\".", filename);
		return -1;
	}

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
			if(!optional)
				cout << "error: could not find string configuration entry \"" << option << "\"" << endl;
		}
		else	{	// option exists, get the value
			val = oit->second;
			rv = 0;
		}
	}

	return rv;
}

/****************************************************************************************/
/*** GetBoolOption() ***/
/****************************************************************************************/
int ConfigReader::GetBoolOption(const char* option, bool& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		// lower case it
		for (int i = 0; i < (int)str.length(); ++i)
			str[i] = tolower(str[i]);

		if (str == "true")
			val = true;
		else if (str == "false")
			val = false;
		else if (atoi(str.c_str()) != 0)
			val = true;
		else
			val = false;

		rv = 0;
	}
	else{
		rv = -1;
		if(!optional)
			cout << "error: could not find boolean configuration entry \"" << option << "\"" << endl;
		}
		
	return rv;
}

/****************************************************************************************/
/*** GetIntOption() ***/
/****************************************************************************************/
int ConfigReader::GetIntOption(const char* option, int& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = atoi(str.c_str());
		rv = 0;
	}
	else{
		rv = -1;
		if(!optional)
			cout << "error: could not find integer configuration entry \"" << option << "\"" << endl;
		}

	return rv;
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
		cout << "error: could not find integer range configuration entry \"" << option << "\"" << endl;
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
int ConfigReader::GetDoubleOption(const char* option, double& val, bool optional /*=false*/)	{
	int rv;
	string str;
	if (GetStringOption(option, str, optional) == 0)	{  // option exists
		val = atof(str.c_str());
		rv = 0;
	}
	else{
		rv = -1;
		if(!optional) cout << "error: could not float configuration entry \"" << option << "\"" << endl;
		}

	return rv;
}

/****************************************************************************************/
/*** GetDoubleRangeOption() ***/
/****************************************************************************************/
int ConfigReader::GetDoubleRangeOption(const char* option, double& val1, double& val2)	{
	int rv;
	string str;
	if (GetStringOption(option, str) == 0)	{  // option exists

		// split up the string
		int len = (int)str.length();
		int i = (int)str.find(' ', 0);
		if (i < 0)
			rv = -1;
		else	{
			val1 = atof(str.substr(0, i).c_str());
			val2 = atof(str.substr(i+1, len-i).c_str());
			rv = 0;
		}
	}
	else{
		rv = -1;
		cout << "error: could not float range configuration entry \"" << option << "\"" << endl;
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

void ConfigReader::TrimWhiteSpace(string& str)	{
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
