#ifndef LIMITS_H
#define LIMITS_H

/* same as bas.maxl, maximum angular number l */
#define MAXL 2 
/* Twice maximum angular number l * 2 */
#define MAXL2 MAXL*2

/* Maximum contraction length of basis functions */
#define MAXG 12
/* Max number of atoms in molecule */
#define MAX_NAT 9
/* Max number of shells in molecule */
#define MAX_NSH 14
/* Max shell an each atom can have */
#define MSHELL 3
/* Max number of atomic orbitals */
#define MAX_NAO 26

#define MAX_NSPIN 1

/* Max number of different atoms in molecule */
#define MAX_NELEM 4


// Steepness of the first counting function
#define KA 10.0

// Steepness of the second counting function
#define KB 20.0

// Offset of the second counting function
#define R_SHIFT 2.0

/* Atoms separated more than this value do not interact */
#define DEFAULT_CUTOFF 25.0

#endif