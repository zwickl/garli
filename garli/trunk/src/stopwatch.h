#ifndef STOPWATCH_H
#define STOPWATCH_H


#ifndef UNIX
#include <time.h>
#else
#include <sys/time.h>
#endif

class Stopwatch	{

	public:
	
		Stopwatch()	{
			Restart();
			}
	
		#ifndef UNIX
		void Restart()	{
			time(&start_time);
			}
		void Start()	{
			time(&start_time);
			}		
		int SplitTime()	{
			time(&end_time);
			return (int)(end_time - start_time);		
			}
		//this is for restarting
		void AddPreviousTime(time_t t){
			start_time -= t;
			}
		#else		
		void Restart()	{
			gettimeofday(&start_time, NULL);
			}
		void Start()	{
			gettimeofday(&start_time, NULL);
			}		
		int SplitTime()	{
			gettimeofday(&end_time, NULL);
			return end_time.tv_sec - start_time.tv_sec;
			}
		//this is for restarting
		void AddPreviousTime(int t){
			start_time.tv_sec -= t;
			}
		#endif


	private:
		#ifndef UNIX
		time_t start_time, end_time;
		#else
		timeval start_time, end_time;
		#endif

};

#endif
