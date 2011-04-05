// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//	do not re-disribute
//	by Paul O. Lewis, Mark T. Holder and David L. Swofford
//	
//
#ifndef __MEMORYACCOUNTANT
#define __MEMORYACCOUNTANT


#if defined (MONITORING_ALLOCATION) && !defined(NDEBUG)
	#include <vector>
	#include <string>
	#include <cassert>
	#include <iostream>
	#include <algorithm>
	#define CREATE_MEMCHK					\
		MemoryAccountant memAccountant;	\
		memTracker = &memAccountant;			\
		memTracker->StartRecording();
	#define MEMCHK_REPORT(o)					\
		memTracker->StopRecording();			\
		memAccountant.Summarize(o);

	using namespace std;
	
	class MemoryAccountant;

	#if defined (INSTANTIATE_MEMCHK)
		//struct dmalloc_t dmalloc;
		MemoryAccountant* memTracker = NULL;
	#else
		//extern struct dmalloc_t dmalloc;
		extern MemoryAccountant* memTracker;
	#endif

	/*----------------------------------------------------------------------------------------------------------------------
	|	Structure used to store accounting information for one individual memory allocation.
	*/
	class MemoryInfo
		{
		public:
		void		*ptr;
		unsigned 	filename_index;
		unsigned 	line_number;
		unsigned 	num_bytes;
		enum	{	Mem_Array 				= 0x01,	/* was allocated using the new [] operator */
					Mem_Op_Err_Free_Array 	= 0x02, /* was allocated using the new, but freed with delete [] */
					Mem_Op_Err_New_Array 	= 0x04, /* was allocated using the new[], but freed with delete */
					Mem_Mid_Array_Delete	= 0x08  /* delete called with an address in the middle of array */
				};
		unsigned	flag;
		bool		array_allocation;
		MemoryInfo(bool is_arr = false);
		};
		
	inline MemoryInfo::MemoryInfo(
	  bool is_arr)	/*true if the memory was allocated using new [] */
		{
		flag = is_arr ? (unsigned) Mem_Array : 0U;
		ptr = (void *)NULL;
		filename_index = 0L;
		line_number = 0L;
		num_bytes = 0L;
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Keeps track of memory allocations and deletions if AddMemoryInfo is called after every call to the new operator and 
	|	MarkDeleted is called for every call to the delete or delete [] operator, which can be done using macros so that 
	|	memory accounting is not done in the release version.
	*/
	class MemoryAccountant
		{
		typedef vector<MemoryInfo>	MemInfoVector;
		typedef vector<string>				FileNameVector;

		FileNameVector		filenames;			/* vector of source file names used to provide filename_index element in the MemoryInfo struct */
		MemoryInfo	tmp;				/* workspace used when adding a new mem_info element */
		MemInfoVector		mem_info;			/* vector of allocation records */
		unsigned			nBad;				/* number of remaining allocs that have not yet been deleted */
		unsigned			nUnknown;			/* number of deletes on non-NULL pointers that are not in our database */
		unsigned			nAllocs;			/* total number of allocs that we have caught with our overload of new */
		unsigned long		numBytesAllocated;	/* total number of bytes from allocs that we have caught with our overload of new */
		unsigned long		currentlyAllocated;	/* the number of bytes that are allocated but haven't been freed (only catchs allocations through our new operator) */
		unsigned long		peakAllocation;		/* the most bytes that are ever allocated at one time (only catchs allocations through our new operator) */
		public:
		static bool			recording;			/* if true, records allocations and deletions; otherwise, ignores them */
			MemoryAccountant();
			~MemoryAccountant();

			void StartRecording();
			void StopRecording();

			void AddMemoryInfo(void *p, string fn, unsigned ln, unsigned long b = 0L, bool is_arr = false);
			void MarkDeleted(void *p, bool is_arr = false);

			void Summarize(ostream &out);
		};

	inline void MemoryAccountant::StartRecording()
		{
		recording = true;
		}

	inline void MemoryAccountant::StopRecording()
		{
		recording = false;
		}

	inline void *operator new (size_t size, const char *file, int line)
		{
		void *p = malloc(size);
		if (p == NULL)
			throw std::bad_alloc();
		if (memTracker != NULL && MemoryAccountant::recording)
			memTracker->AddMemoryInfo(p, file, (unsigned)line, size, false);
		return p;
		}

	inline void operator delete (void *p)
		{
		if (p != NULL)
			{
			if (memTracker != NULL && MemoryAccountant::recording)
				memTracker->MarkDeleted(p, false);
			free(p);
			}
		}

	#if !defined (_MSC_VER)
		inline void * operator new [] (size_t size, const char *file, int line)
			{
		 	void *p = malloc (size);
			if (p == NULL)
				throw std::bad_alloc();
			if (memTracker != NULL && MemoryAccountant::recording)
					memTracker->AddMemoryInfo(p, file, (unsigned) line, size, true);
			return p;
			}

		inline void operator delete [](void *p)
			{
			if (p != NULL)
				{
				if (memTracker != NULL && MemoryAccountant::recording)
					memTracker->MarkDeleted(p, true);
				free(p);
				}
			}
	#endif	//!defined (_MSC_VER)

	#define NEW new (__FILE__, __LINE__)
	#define new NEW
	#if defined (INSTANTIATE_MEMCHK)
		bool MemoryAccountant::recording = false;

		/*----------------------------------------------------------------------------------------------------------------------
		|	Default constructor does nothing currently.
		*/
		MemoryAccountant::MemoryAccountant() 
			: nBad(0L), 
			nUnknown(0L), 
			nAllocs(0L), 
			numBytesAllocated(0L),
			peakAllocation(0L),
			currentlyAllocated(0L)
			{
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	Destructor does nothing currently.
		*/
		MemoryAccountant::~MemoryAccountant()
			{
			filenames.erase(filenames.begin(), filenames.end());
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	Fills in fields of temporary MemoryInfo structure tmp and pushes it onto the mem_info vector.
		*/
		void MemoryAccountant::AddMemoryInfo(
		  void			*p,		/* pointer to the allocated memory */
		  string		fn,		/* name of file where allocation occurred */
		  unsigned		ln,		/* line number where allocation occurred */
		  unsigned long b,		/* number of bytes allocated (defaults to 0L, which means number of bytes is not being tracked) */
		  bool 			is_arr)	/* true if the allocated using the new [] operator */
			{
			if (!MemoryAccountant::recording) 
				return;
			nAllocs++;
			numBytesAllocated += b;
			currentlyAllocated += b;
			if (currentlyAllocated > peakAllocation)
				peakAllocation = currentlyAllocated;
				
			// Attempt to find fn in the stored vector of file names
			//
			FileNameVector::iterator i = find(filenames.begin(), filenames.end(), fn);
			if (i == filenames.end())
				{
				// fn has not previously been encountered
				//
				tmp.filename_index = filenames.size();
				filenames.push_back(fn);
				}
			else
				{
				// fn has not previously been encountered
				//
				tmp.filename_index = (unsigned) (i - filenames.begin());
				}
			
			tmp.ptr = p;
			tmp.flag = is_arr ? (unsigned) MemoryInfo::Mem_Array : 0U;
			tmp.num_bytes = b;
			tmp.line_number = ln;
			mem_info.push_back(tmp);
			nBad++;
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	Finds pointer p in mem_info vector and sets it to NULL, marking it as having been deleted.
		*/
		void MemoryAccountant::MarkDeleted(
		  void *p,		/* the pointer to be marked */
		  bool is_arr)	/*true if the memory is being freed with the delete [] operator */
			{
			if (!recording || p == NULL)
				return;

			bool found = false;
			bool middleArrayError = false;
			MemInfoVector::iterator i;
			for (i = mem_info.begin(); i != mem_info.end(); i++)
				{
				if (i->ptr == p)
					{
					if (is_arr)
						{
						if (!(i->flag & MemoryInfo::Mem_Array))
							i->flag |= MemoryInfo::Mem_Op_Err_Free_Array;
						}
					else
						{
						if (i->flag & MemoryInfo::Mem_Array)
							i->flag |= MemoryInfo::Mem_Op_Err_New_Array;
						}
					i->ptr = NULL;
					nBad--;
					found = true;
					currentlyAllocated -= i->num_bytes;
					break;
					}
				else if (i->ptr < p && ((unsigned long) i->ptr + (unsigned long) i->num_bytes > (unsigned long) p))
					{
					middleArrayError = true;
					i->flag |= MemoryInfo::Mem_Mid_Array_Delete;
					assert(middleArrayError == false);
					}
				}
			assert(! (found && middleArrayError));	// shouldn't be possible to trip the middle of an array error and later find the mem object
			if (!found)
				nUnknown++;
			}

		/*----------------------------------------------------------------------------------------------------------------------
		|	Summarizes memory allocations recorded, displaying total number of allocations, number of allocations currently
		|	not deleted, and file name and line number for remaining undeleted elements.
		*/
		void MemoryAccountant::Summarize(
		  ostream &out)	/* output stream for showing summary */
			{
			out << "\n\nMemory Report" << endl;

			out << "\nTotal allocations: " << nAllocs << " ("<< numBytesAllocated <<" bytes)"<<endl;
			out << "\nLargest Memory Requirement: " << peakAllocation << "  bytes"<<endl;
			if (nBad == 0L)
				out << "0 undeleted memory elements!" << endl;
			else if (nBad == 1L)
				out << "There was 1 undeleted memory element" << endl;
			else
				out << "There were " << nBad << " undeleted memory elements" << endl;

			// Compute approximate amount of memory required to just track allocs
			//
			unsigned z = mem_info.size() * sizeof(MemoryInfo);
			z += sizeof(MemoryAccountant);
			FileNameVector::iterator j;
			for (j = filenames.begin(); j != filenames.end(); j++)
				z += (*j).length();
			out << "The memory tracker itself required at least " << z << " bytes" << endl;

			MemInfoVector::iterator i;
			for (i = mem_info.begin(); i != mem_info.end(); i++)
				{
				MemoryInfo &mi = (*i);
				if (mi.ptr != NULL || mi.flag > 1)
					{
					out << "  ";
					if (mi.num_bytes > 0L)
						out << mi.num_bytes << " bytes ";
					else
						out << "unknown number of bytes ";
					if (mi.ptr != NULL)
						out << "remaining from allocation at ";
					if (mi.flag & MemoryInfo::Mem_Op_Err_Free_Array)
						out << "allocated with new but freed with delete [].  Allocation at ";
					else if (mi.flag & MemoryInfo::Mem_Op_Err_New_Array)
						out << "allocated with new [] but freed with delete.  Allocation at ";
					else if (mi.flag & MemoryInfo::Mem_Mid_Array_Delete)
						out << "allocated but deletion occurred from the middle. Allocation at ";
					string s = filenames[mi.filename_index];
					out << s;
					out << " (" << mi.line_number << ")" << endl;
					}
				}
			}
#		endif //defined (INSTANTIATE_MEMCHK)
#	else	// #if defined(MONITORING_ALLOCATION) && !defined(NDEBUG
#		define CREATE_MEMCHK
#		define MEMCHK_REPORT(o)
#	endif	// #if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)

#endif // #ifndef __MEMORYACCOUNTANT
