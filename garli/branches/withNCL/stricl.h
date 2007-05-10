// stricl.h
// Copyright © 1995 by Dmitri Zaykin (class MyStri)
// and Paul O. Lewis (class NxsMyString)
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
// Associated source code file: "stricl.cpp"
//

#ifndef __STRICL_H
#define __STRICL_H

#include <iosfwd>
#include <fstream>
#include <iostream>
#include <string.h>		// strcmp, stricmp

using namespace std;

const unsigned DEF_STRI_SIZE = 128;
const unsigned DEF_STRI_INCR = 10;

class MyStri
{
	public:
		enum stringErrors { NO_ERR, NULL_STRING, BAD_INDEX };

	protected:
		enum {
			IGNORECASE = 0x01
		};
		MyStri& replaceOneChar(char finchar,
		char replacechar,
		unsigned start = 0);

		int badindex;     // POL was private, but I needed it for NxsMyString
		char *str;         // pointer to characters
		unsigned maxSize;  // total number of bytes allocated
		unsigned len;      // current string size
		unsigned sizeIncr; // increment
		stringErrors strErr;

		// currently only flag is IGNORECASE, which is automatically reset to zero
		// after the next comparison; it is static so that we don't have to worry
		// about whether the correct string's flag was set
		static unsigned flags;

	public:
		static void ignore_case()  { flags |= IGNORECASE; }
		static void respect_case() { if( flags & IGNORECASE ) flags ^= IGNORECASE; }

	public:
		MyStri();
		MyStri(const unsigned size,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		MyStri(const char* s,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		MyStri(const MyStri &,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		virtual ~MyStri(void);
		int				checkBounds(unsigned i) const { return (i < len) ? 1 : 0; }
		void			reallocate();	// cjb - NEVER call this function, its a cheap hack for Parameter serialization
		MyStri&			resize(unsigned new_size);
		void			clear(void) { badindex = -1; }
		char*			getstr(void) const { return str; }
		unsigned	hash1(unsigned modulo) const;
		unsigned	stringlen(void) const { return len; }
		MyStri&			stringupr();
		void			ToUpper();
		stringErrors getError(void)
			{ stringErrors temp = strErr; strErr = NO_ERR; return temp; }
		friend unsigned hash2(MyStri& s, unsigned modulo)
			{ return s.hash1(modulo); }
		MyStri& operator =(const MyStri &);
		MyStri& operator =(const char*);
		MyStri& operator =(const char);
		char& operator[](const unsigned);
		char operator[](const unsigned) const;
		int operator >(const MyStri&) const;
		int operator >=(const MyStri&) const;
		int operator ==(const MyStri&) const;
		int operator <(const MyStri&) const;
		int operator <=(const MyStri&) const;
		int operator !=(const MyStri&) const;
		int operator >(const char*) const;
		int operator >=(const char*) const;
		int operator ==(const char*) const;
		int operator <(const char*) const;
		int operator <=(const char*) const;
		int operator !=(const char*) const;
		friend MyStri operator +(MyStri&, MyStri&);
		friend MyStri operator +(MyStri&, char*);
		friend MyStri operator +(char*, MyStri&);
		friend MyStri operator +(MyStri&, char);
		friend MyStri operator +(char, MyStri&);
		MyStri& operator +=(const MyStri&);
		MyStri& operator +=(const char*);
		MyStri& operator +=(const char);
		friend ostream& operator <<(ostream& os, MyStri& s);
		friend istream& operator >>(istream& is, MyStri& s);
};

class NxsMyString : public MyStri
{
	protected:
		static unsigned tokenpos;
	public:
		NxsMyString() : MyStri() {}
		NxsMyString(const unsigned size
			, const unsigned sizeIncrement = DEF_STRI_INCR)
			: MyStri(size, sizeIncrement) {}
		NxsMyString(const char* s
			, const unsigned sizeIncrement = DEF_STRI_INCR)
			: MyStri(s, sizeIncrement) {}
		NxsMyString(const NxsMyString&
			, const unsigned = DEF_STRI_INCR);

		// removes last char from str iff that char same as cf
		NxsMyString& operator -=(const char);

		// makes a string out of an integer (i.e.,
		//		NxsMyString s = 56;
		// results in s.str being "56"
		NxsMyString& operator =(const int);

		// adds a string made from an integer to *this (i.e.,
		//		NxsMyString s = "Locus-";
		//		s += 4;
		// results in s.str being "Locus-4"
		NxsMyString& operator +=(const int);
		NxsMyString& operator +=(const double i);

		// must override the other += operators since we overrode the one above
		NxsMyString& operator +=( const char );
		NxsMyString& operator +=( const char* );
		NxsMyString& operator +=( const NxsMyString& );

		// removes chars from str that are in trash_set
		char*    cut_trash(char* /*trash_set*/);

		// trims leading and trailing blanks
		void		 trim();

		// conversions to and from underscores
		void		 underscores_to_blanks();
		void		 blanks_to_underscores();

		// these added for compatability with Borland class string
		void     to_upper() { ToUpper(); }
		unsigned length(void) const { return len; }
		char*    c_str(void) const { return str; }

		// added for compatability with Borland template TDDAssociation
		unsigned HashValue() const;

		// for emulation of functions in "string.h" and "iostream.h"
		friend char* Strstr( const NxsMyString&, const char*, int );
		friend char* Strstr( const NxsMyString&, const NxsMyString&, int );
		friend char* MyStristr( const NxsMyString&, const char* );
		friend char* MyStristr( const NxsMyString&, const NxsMyString& );
		friend int Strcmp( const NxsMyString&, const char* );
		friend int Strcmp( const NxsMyString&, const NxsMyString& );
		friend int MyStricmp( const NxsMyString&, const char* );
		friend int MyStricmp( const NxsMyString&, const NxsMyString& );
		friend char* FirstToken( NxsMyString&, const char*, int);
		friend char* NextToken( NxsMyString&, const char*, int);
};

// functions emulating strstr (use MyStristr for case insensitive comparison)
char* Strstr( const NxsMyString& /*str*/, const char* /*substr*/, int /*ignoreCase*/ = 0 );
char* Strstr( const NxsMyString& /*str*/, const NxsMyString& /*substr*/, int /*ignoreCase*/ = 0 );
inline char* MyStristr( const NxsMyString& s, const char* ss ) { return Strstr(s, ss, 1); }
inline char* MyStristr( const NxsMyString& s, const NxsMyString& ss ) { return Strstr(s, ss, 1); }

// functions emulating strcmp and stricmp
inline int Strcmp( const NxsMyString& s, const char* ss) { return strcmp(s.str, ss); }
inline int Strcmp( const NxsMyString& s, const NxsMyString& ss) { return strcmp(s.str, ss.str); }
inline int MyStricmp( const NxsMyString& s, const char* ss) { 
#if defined( USING_UNISTD_H )
	return strcasecmp(s.str, ss); 
#else
	return strcmp(s.str, ss); 
#endif
}
inline int MyStricmp( const NxsMyString& s, const NxsMyString& ss) {
#if defined( USING_UNISTD_H )
	return strcasecmp(s.str, ss.str); 
#else
	return strcmp(s.str, ss.str); 
#endif
}

// functions emulating strtok but with the option to ignore leading
// white space (using isspace) before looking for the delimiters
char* FirstToken(NxsMyString&, const char* /*delims*/, int = 1); // cf. strtok(str, delim)
char* NextToken(NxsMyString&, const char* /*delims*/, int = 1);  // cf. strtok(0, delim)

#endif

