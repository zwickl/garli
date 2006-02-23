#ifndef _SUBSET_
#define _SUBSET_

class subset{
 public:
  int total;
  int seednumber ;
  int element[1024];
  int front[1024];
  double pathlength[1024];
  subset();
  void setseed(int, double dist=0.0);
  void addelement(int,int, double);
  void elementremove(int);
  void removennis();
  void setfront(int,int);
  int getfront(int);
  void compact();
  void clear();
  void print();
};


#endif

