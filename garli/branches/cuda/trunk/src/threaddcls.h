#ifndef THREADDCLS_H
#define THREADDCLS_H

#include <pthread.h>

class MasterGamlConfig;
class Population;

struct transferred_data_t	{
	char *tree_strings;
	unsigned int ts_size;
	double *kappas;
	unsigned int k_size;
	double score;
	int tag;
};

struct thread_arg_t	{
	int nprocs;
	MasterGamlConfig *conf;
	Population *pop;
};

extern transferred_data_t *node_results;
extern pthread_mutex_t lock_pm;
extern pthread_mutex_t lock_pop;
extern pthread_cond_t cond_pm;
extern bool g_quit_time;
extern bool g_processing_message;

void *master_poller(void *varg);
void *thread_func2(void *varg);
void purge_results(transferred_data_t *r);
void copy_results(transferred_data_t *lhs, transferred_data_t rhs);
bool valid_results(transferred_data_t r);
void send_quit_messages(int);
int process_message(char *buf, int size, int who, int tag, thread_arg_t *targ);
void DoMasterSM(char *buf, int size, int who, int tag, thread_arg_t *targ);
void DoMasterAMR(char *buf, int size, int who, int tag, thread_arg_t *targ);

int DoMasterSW(char *buf, int size, int who, int tag, thread_arg_t *targ);

#endif
