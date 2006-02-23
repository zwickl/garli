#undef NO_ERROR
#undef ERROR
#define NO_ERROR	0
#define ERROR		1

#undef FALSE
#undef TRUE
#define FALSE		0
#define TRUE		1

#define NO			0
#define YES			1

#define UP			0
#define DOWN		1

#define	pos1(i,j,n)			((i)*(n)+(j))
#define pos2(i,c,s,k,C,S,K)	((i)*(C)*(S)*(K)+(c)*(S)*(K)+(s)*(K)+(k))
#define pos3(i,j,k,I,J)		((k)*(I)*(J)+(i)*(J)+(j))
