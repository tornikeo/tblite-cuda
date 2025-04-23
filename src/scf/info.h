#ifndef SCF_INFO_H
#define SCF_INFO_H

enum resolution_type
{
  not_used, atom_resolved, shell_resolved, orbital_resolved
};

typedef struct
{
  /* data */
  resolution_type charge = not_used;
  resolution_type dipole = not_used;
  resolution_type quadrupole = not_used;
} scf_info;


#endif