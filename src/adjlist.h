#ifndef ADJLIST_H
#define ADJLIST_H
#include "limits.h"
#include "structure.h"
#include <cassert>

/* type :: adjacency_list
 !> Offset index in the neighbour map
 integer, allocatable :: inl(:)
 !> Number of neighbours for each atom
 integer, allocatable :: nnl(:)
 !> Index of the neighbouring atom
 integer, allocatable :: nlat(:)
 !> Cell index of the neighbouring atom
 integer, allocatable :: nltr(:)
 end type adjacency_list */

 typedef struct
 {
     int inl[MAX_NAT] = {0};  // Offset index in the neighbour map
     int nnl[MAX_NAT] = {0};  // Number of neighbours for each atom
     int nlat[MAX_NAT * 5] = {0}; // Index of the neighbouring atom
     int nltr[MAX_NAT * 5] = {0}; // Cell index of the neighbouring atom
 } adjacency_list;
 
 /*
 subroutine generate(mol, trans, cutoff, inl, nnl, nlat, nltr, complete)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Translation vectors for all images
      real(wp), intent(in) :: trans(:, :)
      !> Realspace cutoff for neighbourlist generation
      real(wp), intent(in) :: cutoff
      !> Offset index in the neighbour map
      integer, intent(inout) :: inl(:)
      !> Number of neighbours for each atom
      integer, intent(inout) :: nnl(:)
      !> Index of the neighbouring atom
      integer, allocatable, intent(out) :: nlat(:)
      !> Cell index of the neighbouring atom
      integer, allocatable, intent(out) :: nltr(:)
      !> Whether a complete or a symmetrical reduced map should be generated
      logical, intent(in) :: complete

 */

 __device__
 void generate(
  const structure_type &mol,
  /* const float (&trans)[3][1], */
  const float &cutoff,
  int (&inl)[MAX_NAT],
  int (&nnl)[MAX_NAT],
  int (&nlat)[5 * MAX_NAT],
  int (&nltr)[5 * MAX_NAT]
  // const bool &complete
 );
 /* 
   !> Create new neighbourlist for a given geometry and cutoff
   subroutine new_adjacency_list(self, mol, trans, cutoff, complete)
      !> Instance of the neighbourlist
      type(adjacency_list), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Translation vectors for all images
      real(wp), intent(in) :: trans(:, :)
      !> Realspace cutoff for neighbourlist generation
      real(wp), intent(in) :: cutoff
      !> Whether a complete or a symmetrical reduced map should be generated
      logical, intent(in), optional :: complete

      logical :: cmplt

      cmplt = .false.
      if (present(complete)) cmplt = complete

      allocate(self%inl(mol%nat), source=0)
      allocate(self%nnl(mol%nat), source=0)
      call generate(mol, trans, cutoff, self%inl, self%nnl, self%nlat, self%nltr, cmplt)
   end subroutine new_adjacency_list
      */
__device__
void new_adjacency_list(
  adjacency_list &self,
  const structure_type &mol,
  /* const float (&trans)[3][1], */
  const float &cutoff,
  const bool complete
);

#endif