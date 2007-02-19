// stricl.h
// Copyright © 1995 by Dmitri Zaykin (class Stri)
// and Paul O. Lewis (class NxsString)
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

class Stri
{
	public:
		enum stringErrors { NO_ERR, NULL_STRING, BAD_INDEX };

	protected:
		enum {
			IGNORECASE = 0x01
		};
		Stri& replaceOneChar(char finchar,
		char replacechar,
		unsigned start = 0);

		int badindex;     // POL was private, but I needed it for NxsString
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
		Stri();
		Stri(const unsigned size,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		Stri(const char* s,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		Stri(const Stri &,
			const unsigned sizeIncrement = DEF_STRI_INCR);
		virtual ~Stri(void);
		int				checkBounds(unsigned i) const { return (i < len) ? 1 : 0; }
		void			reallocate();	// cjb - NEVER call this function, its a cheap hack for Parameter serialization
		Stri&			resize(unsigned new_size);
		void			clear(void) { badindex = -1; }
		char*			getstr(void) const { return str; }
		unsigned	hash1(unsigned modulo) const;
		unsigned	stringlen(void) const { return len; }
		Stri&			stringupr();
		void			ToUpper();
		stringErrors getError(void)
			{ stringErrors temp = strErr; strErr = NO_ERR; return temp; }
		friend unsigned hash2(Stri& s, unsigned modulo)
			{ return s.hash1(modulo); }
		Stri& operator =(const Stri &);
		Stri& operator =(const char*);
		Stri& operator =(const char);
		char& operator[](const unsigned);
		char operator[](const unsigned) const;
		int operator >(const Stri&) const;
		int operator >=(const Stri&) const;
		int operator ==(const Stri&) const;
		int operator <(const Stri&) const;
		int operator <=(const Stri&) const;
		int operator !=(const Stri&) const;
		int operator >(const char*) const;
		int operator >=(const char*) const;
		int operator ==(const char*) const;
		int operator <(const char*) const;
		int operator <=(const char*) const;
		int operator !=(const char*) const;
		friend Stri operator +(Stri&, Stri&);
		friend Stri operator +(Stri&, char*);
		friend Stri operator +(char*, Stri&);
		friend Stri operator +(Stri&, char);
		friend Stri operator +(char, Stri&);
		Stri& operator +=(const Stri&);
		Stri& operator +=(const char*);
		Stri& operator +=(const char);
		friend ostream& operator <<(ostream& os, Stri& s);
		friend istream& operator >>(istream& is, Stri& s);
};

class NxsString : public Stri
{
	protected:
		static unsigned tokenpos;
	public:
		NxsString() : Stri() {}
		NxsString(const unsigned size
			, const unsigned sizeIncrement = DEF_STRI_INCR)
			: Stri(size, sizeIncrement) {}
		NxsString(const char* s
			, const unsigned sizeIncrement = DEF_STRI_INCR)
			: Stri(s, sizeIncrement) {}
		NxsString(const NxsString&
			, const unsigned = DEF_STRI_INCR);

		// removes last char from str iff that char same as cf
		NxsString& operator -=(const char);

		// makes a string out of an integer (i.e.,
		//		NxsString s = 56;
		// results in s.str being "56"
		NxsString& operator =(const int);

		// adds a string made from an integer to *this (i.e.,
		//		NxsString s = "Locus-";
		//		s += 4;
		// results in s.str being "Locus-4"
		NxsString& operator +=(const int);
		NxsString& operator +=(const FLOAT_TYPE i);

		// must override the other += operators since we overrode the one above
		NxsString& operator +=( const char );
		NxsString& operator +=( const char* );
		NxsString& operator +=( const NxsString& );

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
		friend char* Strstr( const NxsString&, const char*, int );
		friend char* Strstr( const NxsString&, const NxsString&, int );
		friend char* Stristr( const NxsString&, const char* );
		friend char* Stristr( const NxsString&, const NxsString& );
		friend int Strcmp( const NxsString&, const char* );
		friend int Strcmp( const NxsString&, const NxsString& );
		friend int Stricmp( const NxsString&, const char* );
		friend int Stricmp( const NxsString&, const NxsString& );
		friend char* FirstToken( NxsString&, const char*, int);
		friend char* NextToken( NxsString&, const char*, int);
};

// functions emulating strstr (use Stristr for case insensitive comparison)
char* Strstr( const NxsString& /*str*/, const char* /*substr*/, int /*ignoreCase*/ = 0 );
char* Strstr( const NxsString& /*str*/, const NxsString& /*substr*/, int /*ignoreCase*/ = 0 );
inline char* Stristr( const NxsString& s, const char* ss ) { return Strstr(s, ss, 1); }
inline char* Stristr( const NxsString& s, const NxsString& ss ) { return Strstr(s, ss, 1); }

// functions emulating strcmp and stricmp
inline int Strcmp( const NxsString& s, const char* ss) { return strcmp(s.str, ss); }
inline int Strcmp( const NxsString& s, const NxsString& ss) { return strcmp(s.str, ss.str); }
inline int Stricmp( const NxsString& s, const char* ss) { 
#if defined( USING_UNISTD_H )
	return strcasecmp(s.str, ss); 
#else
	return strcmp(s.str, ss); 
#endif
}
inline int Stricmp( const NxsString& s, const NxsString& ss) {
#if defined( USING_UNISTD_H )
	return strcasecmp(s.str, ss.str); 
#else
	return strcmp(s.str, ss.str); 
#endif
}

// functions emulating strtok but with the option to ignore leading
// white space (using isspace) before looking for the delimiters
char* FirstToken(NxsString&, const char* /*delims*/, int = 1); // cf. strtok(str, delim)
char* NextToken(NxsString&, const char* /*delims*/, int = 1);  // cf. strtok(0, delim)

#endif

