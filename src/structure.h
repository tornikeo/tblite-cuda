#ifndef STRUCTURE_H
#define STRUCTURE_H
#include "limits.h"

typedef struct {
  int nat;
  int nid;
  int nbd;
  bool periodic = false; 
  int id[MAX_NAT];
  int num[MAX_NELEM];
  float xyz[MAX_NAT][3];
  float charge = 0;
  int uhf = 0;
} structure_type;

#endif