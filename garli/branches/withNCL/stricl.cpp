// set.h
// Copyright © 1998 by Paul O. Lewis
// All rights reserved.
//
// This code may be used and modified for non-commercial purposes
// but redistribution in any form requires written permission.
// Please contact:
//
//   Paul O. Lewis, Assistant Professor
//   167 Castetter Hall
//   Department of Biology
//   The University of New Mexico
//   Albuquerque, NM 87131-1091
//   Phone: (505) 277-6681
//   Fax: (505) 277-0304
//   email: lewisp@unm.edu
//   http://biology.unm.edu/~lewisp/pol.html
//
// Note: moving January 1, 1999, to the Department of Ecology and
// Evolutionary Biology, University of Connecticut
//
// Associated header file: "stricl.h"
//

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <cassert>

using namespace std;

#include "stricl.h"
#include "defs.h"

unsigned MyStri::flags = 0;
unsigned NxsMyString::tokenpos = 0;

MyStri::MyStri()
{
	 len = 0;
	 maxSize = DEF_STRI_SIZE;

	 //str = new char[maxSize];
	 MEM_NEW_ARRAY(str,char,maxSize);

	 *str = '\0';
	 strErr = NO_ERR;
	 badindex = -1;
	 sizeIncr = DEF_STRI_INCR;
}

MyStri::MyStri(const unsigned size, const unsigned sizeIncrement)
{
	 len = 0;
	 maxSize = (size < 2) ? 2 : size;

	 //str = new char[maxSize];
	 MEM_NEW_ARRAY(str,char,maxSize);

	 *str = '\0';
	 strErr = NO_ERR;
	 badindex = -1;
	 sizeIncr = (sizeIncrement == 0) ? 1 : sizeIncrement;
}

MyStri::MyStri(const char* s, const unsigned sizeIncrement)
{
	 if (s != NULL && s[0] != '\0') {
		 len = (int)strlen(s);
		 maxSize = len+1;

		 //str = new char[maxSize];
		 MEM_NEW_ARRAY(str,char,maxSize);

		 strcpy(str, s);
		 strErr = NO_ERR;
	 }
	 else {
		 len = 0;
		 maxSize = 2;

		 //str = new char[maxSize]; //POL
		 MEM_NEW_ARRAY(str,char,maxSize);

		 *str = '\0';    //POL
		 strErr = NULL_STRING;
	 }
	 badindex = -1;
	 sizeIncr = (sizeIncrement == 0) ? 1 : sizeIncrement;
}

MyStri::MyStri(const MyStri &s, const unsigned sizeIncrement)
{
	if( s.len > 0 ) {
		len = s.len;
		maxSize = len+1;

		//str = new char[maxSize];
		MEM_NEW_ARRAY(str,char,maxSize);

		strcpy(str, s.str);
	 strErr = NO_ERR;
	 badindex = -1;
  }
  else {
	 len = 0;
	 maxSize = 2;

	 //str = new char[maxSize]; //POL
	 MEM_NEW_ARRAY(str,char,maxSize);

	 *str = '\0'; //POL
	 strErr = NULL_STRING;
	 badindex = -1;
  }
  sizeIncr = (sizeIncrement == 0) ? 1 : sizeIncrement;
}

MyStri::~MyStri(void)
{
	len = 0;
	assert( str );

	//delete [] str;
	MEM_DELETE_ARRAY(str);	//POL-1/3/98 (no. elem. = maxSize)

	maxSize = 0;
	sizeIncr = 0;
}

void MyStri::ToUpper()
{
	char* s = str;
	while(*s) {
		*s = (char)toupper(*s); //build-5> added (char) cast
		++s;
	}
}

unsigned MyStri::hash1(unsigned modulo) const
// return hashing function for a given table size
{
   unsigned y = 0;
   char* p = str;

   while (*p)
      y = (y + len + 13 * *p++) % modulo;
   return y;
}

void MyStri::reallocate()
{
	 len = 0;
	 maxSize = DEF_STRI_SIZE;

	 //str = new char[maxSize];
	 MEM_NEW_ARRAY(str,char,maxSize);

	 *str = '\0';
	 strErr = NO_ERR;
	 badindex = -1;
	 sizeIncr = DEF_STRI_INCR;
}

MyStri& MyStri::resize(unsigned new_size)
// resize the string object to new_size+1 bytes
{
   char* strcopy;
	unsigned current_len = len;

   // if new_size is not a perfect multiple of sizeIncr,
   // increase new_size
   if (((new_size / sizeIncr) * sizeIncr) < new_size)
       new_size = ((new_size / sizeIncr) + 1) * sizeIncr;

	 if (len >= new_size || maxSize >= new_size) return *this;

	 //strcopy = new char[len+1];
	 MEM_NEW_ARRAY(strcopy,char,len+1);

	 strcpy(strcopy, str);

	 //delete [] str;
	 MEM_DELETE_ARRAY(str);	//POL-1/3/98 (no. elem. = maxSize)

	 maxSize = new_size;

	 //str = new char[maxSize];
	 MEM_NEW_ARRAY(str,char,maxSize);

	 len = current_len;
	 strcpy(str, strcopy);

	 //delete [] strcopy;
	 MEM_DELETE_ARRAY(strcopy);	//POL-1/3/98 (no. elem. = len+1)

	 return *this;
}

MyStri& MyStri::stringupr()
{
	ToUpper();
	return *this;
}

MyStri& MyStri::operator =(const MyStri &s)
{
	 len = s.len;
	 if ((len+1) > maxSize)
			 resize(len+1);
	 *str = '\0';
	 if( s.str ) strcpy(str, s.str);
	 return *this;
}

MyStri& MyStri::operator =(const char* s)
{
	len = (s ? (int)strlen(s) : 0);
	if ((len+1) > maxSize)
		resize(len+1);
	*str = '\0';
	if( s ) strcpy(str, s);
	return *this;
}

MyStri& MyStri::operator =(const char c)
{
  len = 1;
  *str = c;
  *(str+1) = '\0';
  return *this;
}

char& MyStri::operator [](const unsigned i)
{
	if (checkBounds(i)) {
		return *(str + i);
	}
	else {
		strErr = BAD_INDEX;
		return (char&)(badindex=0);
	}
}

char MyStri::operator [](const unsigned i) const
{
	if (checkBounds(i)) {
		return *(str + i);
	}

	return '\0';
}

int MyStri::operator >(const MyStri& s) const
{
	if( flags & IGNORECASE ) {
	  flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s.str) > 0;
#else
	  return strcmp(str,s.str) > 0;
#endif
   }
   else
		 return strcmp(str,s.str) > 0;
}

int MyStri::operator >=(const MyStri& s) const
{
	 if( flags & IGNORECASE ) {
		 flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
		 return strcasecmp(str,s.str) >= 0;
#else
		 return strcmp(str,s.str) >= 0;
#endif
	 }
	 else
		 return strcmp(str,s.str) >= 0;
}

int MyStri::operator ==(const MyStri& s) const
{
	 if( flags & IGNORECASE ) {
		 flags ^= IGNORECASE;
#if defined( USING_UNISTD_H ) 
		 return strcasecmp(str,s.str) == 0;
#else
		 return strcmp(str,s.str) == 0;
#endif
	 }
	 else
		 return strcmp(str,s.str) == 0;
}

int MyStri::operator <(const MyStri& s) const
{
	 if( flags & IGNORECASE ) {
		 flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
		 return strcasecmp(str,s.str) < 0;
#else
		 return strcmp(str,s.str) < 0;
#endif
	 }
	 else
		 return strcmp(str,s.str) < 0;
}

int MyStri::operator <=(const MyStri& s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s.str) <= 0;
#else
	  return strcmp(str,s.str) <= 0;
#endif
   }
   else
     return strcmp(str,s.str) <= 0;
}

int MyStri::operator !=(const MyStri& s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s.str) != 0;
#else
	  return strcmp(str,s.str) != 0;
#endif
   }
   else
     return strcmp(str,s.str) != 0;
}

int MyStri::operator >(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s) > 0;
#else
	  return strcmp(str,s) > 0;
#endif
   }
   else
     return strcmp(str,s) > 0;
}

int MyStri::operator >=(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp( str, (char*)s ) >= 0;
#else
	  return strcmp(str,s) >= 0;
#endif
   }
   else
     return strcmp(str,s) >= 0;
}

int MyStri::operator ==(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s) == 0;
#else
	  return strcmp(str,s) == 0;
#endif
   }
   else
     return strcmp(str,s) == 0;
}

int MyStri::operator <(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H ) 
	  return strcasecmp(str,s) < 0;
#else
	  return strcmp(str,s) < 0;
#endif
   }
   else
     return strcmp(str,s) < 0;
}

int MyStri::operator <=(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s) <= 0;
#else
	  return strcmp(str,s) <= 0;
#endif
   }
   else
     return strcmp(str,s) <= 0;
}

int MyStri::operator !=(const char* s) const
{
   if( flags & IGNORECASE ) {
     flags ^= IGNORECASE;
#if defined( USING_UNISTD_H )
	  return strcasecmp(str,s) != 0;
#else
	  return strcmp(str,s) != 0;
#endif
   }
   else
     return strcmp(str,s) != 0;
}

MyStri operator +(MyStri &s1, MyStri &s2)
{
	 MyStri result(s1.len + s2.len + 1);

	 if( s1.str ) strcpy(result.str, s1.str);
	 if( s2.str ) strcat(result.str, s2.str);
	 result.len = s1.len + s2.len;

	 return result;
}

MyStri operator +(MyStri &s1, char* s2)
{
	 MyStri result(s1.len + strlen(s2) + 1);
	 unsigned length = (unsigned) strlen(s2);

	 if( s1.str ) strcpy(result.str, s1.str);
   if( s2 ) strcat(result.str, s2);
   result.len = s1.len + length;
   return result;
}

MyStri operator +(char* s1, MyStri &s2)
{
   MyStri result((unsigned) strlen(s1) + s2.len + 1);
   unsigned length = (unsigned) strlen(s1);

   if( s1 ) strcpy(result.str, s1) ;
   if( s2.str ) strcat(result.str, s2.str);
   result.len = s2.len + length;
   return result;
}

MyStri operator +(MyStri &s, char c)
{
   MyStri result(s.len + 2);

   if( s.str ) strcpy(result.str,s.str);
   result.len = s.len;
   *(result.str+result.len+1) = '\0';
   *(result.str+result.len) = c;
   result.len++;
   return result;
}

MyStri operator +(char c, MyStri &s)
{
   MyStri result(2 + s.len) ;

   *(result.str+1) = '\0';
   *(result.str) = c;
   if( s.str ) strcat(result.str, s.str);
   result.len = s.len + 1;
   return result;
}

MyStri& MyStri::operator +=(const MyStri &s)
{
	 unsigned newlen = len + s.len;

	 if ((newlen+1) > maxSize)
			 resize(newlen+1);
	 if( s.str ) strcat(str, s.str);
	 len += s.len;

	 return *this;
}

MyStri& MyStri::operator +=(const char* s)
{
	 unsigned newlen = len + strlen(s);

	 if ((newlen+1) > maxSize)
			 resize(newlen+1);
	 if( s ) strcat(str, s);
	 len += (unsigned) strlen(s);
	 return *this;
}

MyStri& MyStri::operator +=(const char c)
{
	 unsigned newlen = len + 1;

	 if ((newlen+1) > maxSize)
			 resize(newlen+1);
	 str[len] = c;
	 str[len++ +1] = '\0';

	 return *this;
}

ostream& operator <<(ostream& os, MyStri& s)
{
		os << s.str;
		return os;
}

istream& operator >>(istream& is, MyStri& s)
{
		if (s.str) {
			//delete [] s.str; // delete current string
			MEM_DELETE_ARRAY(s.str);	//POL-1/3/98 (no. elem. = s.maxSize)
		}
		//BUGBUG: this needs work!!
		// read the maximum and current string sizes
		is >> s.maxSize >> s.len;

		// allocate the space for the new string
		//s.str = new char[s.maxSize];
		MEM_NEW_ARRAY(s.str,char,s.maxSize);

		// read the new string
		is.getline(s.str, s.len+1);
		return is;
}

NxsMyString::NxsMyString(const NxsMyString &s, const unsigned sizeIncrement)
		: MyStri (s.str, sizeIncrement)
{
	//len = s.len;
	//maxSize = len+1;
	//str = new char[maxSize];
	//strcpy(str, s.str);
	//sizeIncr = (sizeIncrement == 0) ? 1 : sizeIncrement;
	//strErr = NO_ERR;
	//badindex = -1;
}

NxsMyString& NxsMyString::operator -=(const char c)
{
	 if(str[len-1] == c)
		 str[--len] = '\0';

	 return *this;
}

NxsMyString& NxsMyString::operator =(const int i)
{	
	char s[6];
	//itoa(i, s, 10);
//	sprintf(s, "%d", i);
	*s=i;
	// below modified from MyStri::operator =(const char*)
	len = 1;
//	if ( ( len + 1 ) > maxSize )
//		resize( len + 1 );
	*(s+1) = '\0';
	strcpy( str, s );

	return *this;
}

NxsMyString& NxsMyString::operator +=(const int i)
{
	char s[6];
	sprintf(s, "%d", i);  

	// below modified from MyStri::operator +=(const MyStri&)
	int slen = (int)strlen(s);
	unsigned newlen = len + slen;

	if( ( newlen + 1 ) > maxSize )
		resize( newlen + 1 );
	strcat( str, s );

	len += slen;

	return *this;
}

NxsMyString& NxsMyString::operator +=(const double i)
{
	char s[40];
	sprintf(s, "%lf", i);  

	// below modified from MyStri::operator +=(const MyStri&)
	int slen = (int)strlen(s);
	unsigned newlen = len + slen;

	if( ( newlen + 1 ) > maxSize )
		resize( newlen + 1 );
	strcat( str, s );

	len += slen;

	return *this;
}

NxsMyString& NxsMyString::operator +=( const NxsMyString& s )
{
	// borrowed from MyStri::operator +=(const MyStri&)
	unsigned newlen = len + s.len;

	if ((newlen+1) > maxSize)
		resize(newlen+1);
	if( s.str ) strcat(str, s.str);
	len += s.len;

	return *this;
}

NxsMyString& NxsMyString::operator +=( const char* s )
{
	if( !s || s[0] == '\0' ) return *this;

	// borrowed from MyStri::operator +=(const char*)
	int slen = (int)strlen(s);
	unsigned newlen = len + slen;

	if ( ( newlen + 1 ) > maxSize )
		resize( newlen + 1 );
	strcat( str, s );
	len += slen;

	return *this;
}

NxsMyString& NxsMyString::operator +=( const char ch )
{
	// borrowed from MyStri::operator +=(const char)
	unsigned newlen = len + 1;

	if( ( newlen + 1 ) > maxSize )
		resize( newlen + 1 );
	str[len] = ch;
	str[len++ +1] = '\0';

	return *this;
}

char* NxsMyString::cut_trash(char* trashSet)
{
	if( !str || !trashSet)
		return NULL;
	char* tough = str;
	char* r     = str;
	while( *r ) {
		if( strchr(trashSet, *r) ) { r++; len--; continue; }
		*tough++ = *r++;
	}
	*tough = '\0';
	return str;
}

void NxsMyString::trim()
{
	for(;;) {
		char ch = str[len-1];
		if( ch != ' ' ) break;
		str[--len] = '\0';
	}

	char* p = str;
	while( p && *p == ' ' ) p++;
	if( p != str ) {
		char* q = str;
		len -= (unsigned)(p - q);
		while( *p ) { *q++ = *p++; }
		*q = *p;	// copy the '\0'
	}
}

void NxsMyString::underscores_to_blanks()
{
	char* p;
	for( p = str; *p; p++) {
		if( *p == '_' ) *p = ' ';
	}
}

void NxsMyString::blanks_to_underscores()
{
	char* p;
	for( p = str; *p; p++) {
		if( *p == ' ' ) *p = '_';
	}
}

//
// cf. strtok(str, delim)
//
char* FirstToken(NxsMyString& Str, const char* delims, int ignoreLeadingWS /* = 1 */)
{
	Str.tokenpos = 0;
	char* start = Str.str;
	char* s     = Str.str;

	if( ignoreLeadingWS ) {
		for(; *s && isspace(*s); s++) {
			start++;
			Str.tokenpos++;
		}
	}

	for(; *s && !strchr(delims, *s); s++ )
		Str.tokenpos++;

	if( *s ) {
		*s = '\0';
		Str.tokenpos++;
	}

	// Str.tokenpos now indicates beginning of next token
	// _or_ the terminal null character in Str.str
	return start;
}

//
// cf. strtok(0, delim)
//
char* NextToken(NxsMyString& Str, const char* delims, int ignoreLeadingWS /* = 1 */)
{
	if( Str.tokenpos > Str.len ) {
		Str.tokenpos = Str.len;
		return NULL;
	}

	char* start = Str.str + Str.tokenpos;
	char* s     = start;

	if( ignoreLeadingWS ) {
		for(; *s && isspace(*s); s++) {
			start++;
			Str.tokenpos++;
		}
	}

	for(; *s && !strchr(delims, *s); s++ )
		Str.tokenpos++;

	if( *s ) {
		*s = '\0';
		Str.tokenpos++;
	}

	return ( *start ? start : 0 );
}

//
// emulation of strstr - looks for substr within str
// if found, returns pointer to substr in str
// if not found, returns NULL
// if ignoreCase = 1, comparison is case insensitive
//
char* Strstr( const NxsMyString& Str, const char* substr, int ignoreCase /* = 0 */)
{
	if( !substr ) return NULL;

	if( !ignoreCase )
		return ::strstr(Str.str, substr);

	// make copies of Str.str and substr
	NxsMyString Str_tmp = Str;
	NxsMyString substr_tmp = substr;

	// convert both to all upper case
	Str_tmp.to_upper();
	substr_tmp.to_upper();

	return ::strstr(Str_tmp.str, substr_tmp.str);
}

char* Strstr( const NxsMyString& Str, const NxsMyString& substr, int ignoreCase /* = 0 */)
{
	return Strstr(Str, substr.str, ignoreCase);
}

unsigned NxsMyString::HashValue() const
{
	 unsigned y = 0;
	 char* p = str;

	 while (*p) {
			y  += (y + len + 13 * *p++);
	 }

	 return y;
}


