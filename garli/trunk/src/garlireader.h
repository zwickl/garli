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

//  This file was adapted from from the BasicCmdLine example provided as
//  part of the NCL

//	Copyright (C) 1999-2002 Paul O. Lewis
//
//	This file is part of NCL (Nexus Class Library).
//
//	NCL is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	NCL is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with NCL; if not, write to the Free Software Foundation, Inc., 
//	59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef NCL_GarliReader_H
#define NCL_GarliReader_H

#define COMMAND_MAXLEN  255

#include "ncl.h"
#include "nxsmultiformat.h"

class ModelSpecification;

//the reader is no longer derived from NexusBlock itself which was done such that it was it's own
//custom block (a bit weird).  Garli block is separate entity now.
class GarliBlock: public NxsBlock{
    public:   
        GarliBlock():NxsBlock(){
			id ="GARLI";
			}
		NxsString modelString;
		char *next_command;
        void Read(NxsToken &token);
        void HandleEndblock(NxsToken &token);
		const NxsString GetModelString(){return modelString;}
		bool ModelStringWasRead(){return modelString.empty() == false;}
		void Clear() { modelString.clear();}
		void Report(ostream &out) const;
	//	void HandleNextCommand();
//		void NexusError(NxsString msg, file_pos pos, long line, long col);
    }; 

/*----------------------------------------------------------------------------------------------------------------------
|	GarliReader provides a template for creating a program that reads NEXUS data files and provides a basic command 
|	line. After compiling GarliReader, you will already have a program that understands the following commands, either 
|	typed in at the console or provided in a GarliReader block in a NEXUS data file (exception is the execute command,
|	which can only be entered at the console). Keywords in the descriptions below are given in uppercase, however the
|	commands themselves are case-insensitive. Lower-case indicates a parameter supplied by the user (e.g., "filename" 
|	would be replaced by the actual name of the file). Square brackets indicate optional keywords or subcommands.
|>
|	EXECUTE filename;
|	
|	LOG [options];
|	
|	  Option         Action
|	  ------------------------------------------------------
|	  FILE=filename  specifies name of log file to start
|	  START          indicates logging is to be started
|	  STOP           indicates logging is to be stopped
|	  APPEND         append to log file if it already exists
|	  REPLACE        replace log file without asking
|	
|	QUIT;
|>
|	See the Read function for details and to add other commands.
|	
|	To change the name of the program (which is also the prompt name and the name of the program's private NEXUS 
|	block), replace all occurrences of GarliReader with the name of your program (also search for the string 
|	"GarliReader" and replace with an appropriate string at each occurrence).
|	
|	This class handles reading and storage for the NxsReader block GarliReader. It also serves as the main class for 
|	the program GarliReader, acting as both a NxsReader object (in order to be capable of parsing data files) as well 
|	as a NxsBlock object (in order to be able to process commands in a GarliReader block). 
|	
|	Adding a new data member? Don't forget to:
|~
|	o Describe it in the class header comment at the top of "GarliReader.h"
|	o Initialize it (unless it is self-initializing) in the constructor and reinitialize it in the Reset function
|	o Describe the initial state in the constructor documentation
|	o Delete memory allocated to it in both the destructor and Reset function
|	o Report it in some way in the Report function
|~
*/

class GarliReader
  : public MultiFormatReader
	{
	friend class NxsBlock;
	public:
		static GarliReader & GetInstance();
		enum UserQueryEnum		/* enumeration used with UserQuery member function to specify which choices to provide the user */
			{
			uq_cancel = 0x01,	/* provide opportunity to cancel */
			uq_ok	  = 0x02,	/* provide opportunity to answer ok */
			uq_yes	  = 0x04,	/* provide opportunity to answer yes */
			uq_no 	  = 0x08	/* provide opportunity to answer no */
			};

							GarliReader();
		virtual				~GarliReader();

		bool				EnteringBlock(NxsString blockName);
		void				ExitingBlock(NxsString blockName);
		void				ExecuteStarting();
		void				ExecuteStopping();
		void				OutputComment(const NxsString &msg);
		void				HandleNextCommand();
		void				NexusError(NxsString msg, file_pos pos, long line, long col){
							NexusError(msg, pos, line, col, true);
							}
		void				NexusError(NxsString msg, file_pos pos, long line, long col, bool throwExcept);
		void				PreprocessNextCommand();
		void				PrintMessage(bool linefeed = true);
	//	virtual void		Report(ostream &out);
		void				Run(char *infile_name);
		void				SkippingBlock(NxsString blockName);
		void				SkippingCommand(NxsString commandName);
		void				SkippingDisabledBlock(NxsString blockName);
		virtual bool		UserQuery(NxsString mb_message, NxsString mb_title, GarliReader::UserQueryEnum mb_choices = GarliReader::uq_ok);

		//a bunch of hacky stuff here got removed when going to the new NLC Factory API and deriving the reader
		//from Multiformat Reader- I don't need to worry about multiple charblocks and such myself.  Many functions
		//that were here are now further up in the inheritance chain and not part of my code
		//char blocks
	
	protected:
		bool				inf_open;			/* true if `inf' is currently open */
		bool				logf_open;			/* true if `logf' is currently open */
		bool				quit_now;			/* set to false at beginning of Run and turns true only when QUIT command processed */
		ofstream			logf;				/* the log file output stream */
		NxsString			message;			/* workspace for composing output strings */
		//none of these should be getting used with the new Factory/MultiformatReader system
		NxsTreesBlock		*trees;				/* pointer to NxsTreesBlock object */
		NxsTaxaBlock		*taxa;				/* pointer to NxsTaxaBlock object */
		NxsAssumptionsBlock	*assumptions;		/* pointer to NxsAssumptionsBlock object */
		NxsDistancesBlock	*distances;			/* pointer to NxsDistancesBlock object */
		NxsCharactersBlock	*characters;		/* pointer to NxsCharactersBlock object */
		NxsDataBlock		*data;				/* pointer to NxsDataBlock object */
		//this still is being used
		GarliBlock			*garliBlock;

		NxsString			errormsg;
		char				*next_command;		/* workspace for processing next command entered interactively by user */

		unsigned			CharLabelToNumber(NxsString s) const;
		bool				FileExists(const char* fn) const;
		bool				FileIsNexus(const char *name) const;
		bool				FileIsFasta(const char *name) const;
		int					GetToken( FILE *in, char* tokenbuf, int maxlen) const;
		NxsString			GetFileName(NxsToken& token);
		void				HandleEndblock(NxsToken& token);
		void				HandleShow(NxsToken& token);
		void				HandleHelp(NxsToken& token);
		void				HandleLog(NxsToken& token);
		void				HandleExecute(NxsToken& token);
		void				HandleGarliReader(NxsToken &token);

	public:
		int				HandleExecute(const char *filename, bool purge);
		string			GetModelString(){
							return garliBlock->GetModelString();
							}
		bool FoundModelString() {return garliBlock->ModelStringWasRead();}
		void ClearModelString() {garliBlock->Clear();}

		//this removes and deallocates everything in the reader and gets it ready
		//for further reading
		void ClearContent();
		bool ReadData(const char* filename, const ModelSpecification &modspec);
		const NxsCharactersBlock *CheckBlocksAndGetCorrectCharblock(const ModelSpecification &modspec) const;
		static void GetDefaultIntWeightSet(const NxsCharactersBlock *charblock, vector<int> &wset);
		};

/*----------------------------------------------------------------------------------------------------------------------
|	The MyNexusToken class provides a NxsToken-derived object that can display output comments as it encounters them.
|	The virtual function NxsToken::OutputComment is overridden in this class for this purpose.
*/
class MyNexusToken
  : public NxsToken
	{
	public:
				MyNexusToken(istream &i);

		void	OutputComment(const NxsString &msg);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Will be called by NxsReader::Execute after the initial "#NEXUS" keyword is found in a NEXUS file but before other
|	tokens are read. Add code here if you need to do any initializations prior to encountering any NEXUS blocks in a
|	NEXUS data file.
*/
inline void GarliReader::ExecuteStarting()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Will be called by NxsReader::Execute just before it exits after reading to the end of a NEXUS data file (or until
|	encountering a LEAVE command between NEXUS blocks. Add code here if you need to clean up any memory allocated in
|	ExecuteStarting.
*/
inline void GarliReader::ExecuteStopping()
	{
	}


#endif

