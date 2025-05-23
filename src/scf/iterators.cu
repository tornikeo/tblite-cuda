#include <stdio.h>
#include <cassert>
#include "iterators.h"
#include "../utils/array.h"
#include "../lapack/sygvd.h"
#include "../wavefunction/fermi.h"
#include "../wavefunction/mulliken.h"
#include "potential.h"

__device__
// constexpr 
int get_mixer_dimension(
  const structure_type &mol, 
  const basis_type &bas,
  const scf_info &info
)
{
  int ndim = 0;

  switch (info.charge) {
    case atom_resolved:
      ndim += mol.nat;
      break;
    case shell_resolved:
      ndim += bas.nsh;
      break;
  }

  switch (info.dipole) {
    case atom_resolved:
      ndim += 3 * mol.nat;
      break;
  }

  switch (info.quadrupole) {
    case atom_resolved:
      ndim += 6 * mol.nat;
      break;
  }

  return ndim;
}

/*
subroutine set_mixer(mixer, wfn, info)
   use tblite_scf_info, only : atom_resolved, shell_resolved
   type(broyden_mixer), intent(inout) :: mixer
   type(wavefunction_type), intent(in) :: wfn
   type(scf_info), intent(in) :: info

   select case(info%charge)
   case(atom_resolved)
      call mixer%set(wfn%qat)
   case(shell_resolved)
      call mixer%set(wfn%qsh)
   end select

   select case(info%dipole)
   case(atom_resolved)
      call mixer%set(wfn%dpat)
   end select

   select case(info%quadrupole)
   case(atom_resolved)
      call mixer%set(wfn%qpat)
   end select
end subroutine set_mixer
*/

__device__ 
void set_mixer(
  broyden_mixer &mixer,
  const wavefunction_type &wfn,
  const scf_info &info
)
{
  // printf("info.charge: %d\n", info.charge);
  // printf("info.dipole: %d\n", info.dipole);
  // printf("info.quadrupole: %d\n", info.quadrupole);
  switch (info.charge)
  {
  case atom_resolved:
    mixer.set(&wfn.qat[0][0], MAX_NAT);
    break;
  case shell_resolved:
    mixer.set(&wfn.qsh[0][0], MAX_NSH);
  default:
    break;
  }
  switch (info.dipole)
  {
  case atom_resolved:
    mixer.set(&wfn.dpat[0][0][0], 3 * MAX_NAT);
    break;
  default:
    break;
  }

  switch (info.quadrupole)
  {
  case atom_resolved:
    mixer.set(&wfn.qpat[0][0][0], 6 * MAX_NAT);
    break;
  default:
    break;
  }

}

/*
subroutine get_qat_from_qsh(bas, qsh, qat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: qsh(:, :)
   real(wp), intent(out) :: qat(:, :)

   integer :: ish, ispin

   qat(:, :) = 0.0_wp
   ! $omp parallel do schedule(runtime) collapse(2) default(none) &
   ! $omp reduction(+:qat) shared(bas, qsh) private(ish)
   do ispin = 1, size(qsh, 2)
      do ish = 1, size(qsh, 1)
         qat(bas%sh2at(ish), ispin) = qat(bas%sh2at(ish), ispin) + qsh(ish, ispin)
      end do
   end do
end subroutine get_qat_from_qsh

*/

__device__
void get_qat_from_qsh
(
 const basis_type &bas,
 const float (&qsh)[MAX_NSPIN][MAX_NSH],
 float (&qat)[MAX_NSPIN][MAX_NAT]  
)
{
  for (int ispin = 0; ispin < MAX_NSPIN; ++ispin) {
    for (int ish = 0; ish < MAX_NSH; ++ish) {
      int atom_index = bas.sh2at[ish];
      qat[ispin][atom_index] += qsh[ispin][ish];
    }
  }
}

__device__
void diff_mixer(
  broyden_mixer &mixer,
  const wavefunction_type &wfn,
  const scf_info &info
)
{
  switch (info.charge)
  {
  case atom_resolved:
    assert(false&&"Unimplemented");
    break;
  case shell_resolved:
    mixer.diff(&wfn.qsh[0][0], SIZEOF_ARRAY(wfn.qsh));
  default:
    break;
  }
}

/*

subroutine get_density(wfn, solver, ints, ts, error)
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(inout) :: wfn
   !> Solver for the general eigenvalue problem
   class(solver_type), intent(inout) :: solver
   !> Integral container
   type(integral_type), intent(in) :: ints
   !> Electronic entropy
   real(wp), intent(out) :: ts
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
...
*/
__device__
void get_density(
  wavefunction_type &wfn,
  sygvd_solver &solver,
  const integral_type &ints,
  float &ts
  /* error &err*/
)
{
  /*    
  real(wp) :: e_fermi, stmp(2)
  real(wp), allocatable :: focc(:)
  integer :: spin
   */
  float e_fermi = 0;
  float stmp[2] = {0};
  float focc[MAX_NAO] = {0};
  int spin = 0;

  switch (wfn.nspin)
  {
  case 1:
    solver.solve(
      wfn.coeff[0], 
      ints.overlap,
      wfn.emo[0]
      /*error,*/
    );
    /* wfn%focc(:, :) = 0.0_wp */
    for (int i = 0; i < MAX_NSPIN; ++i) {
      for (int j = 0; j < MAX_NAO; j++)
      {
        wfn.focc[i][j] = 0.0f;
      }
    }


    /* do spin = 1, 2
         call get_fermi_filling(wfn%nel(spin), wfn%kt, wfn%emo(:, 1), &
            & wfn%homo(spin), focc, e_fermi, stmp(spin))
         wfn%focc(:, 1) = wfn%focc(:, 1) + focc
      end do
    */
   for (int spin = 0; spin < 2; spin++)
   {
     get_fermi_filling(
       wfn.nel[spin], 
       wfn.kt,
       wfn.emo[0],
       wfn.homo[spin],
       focc,
       e_fermi,
       stmp[spin]
     );
     for (int i = 0; i < MAX_NAO; i++)
     {
      wfn.focc[0][i] += focc[i];
     }
   }
   ts += (stmp[0] + stmp[1]);
   /*       call get_density_matrix(wfn%focc(:, 1), wfn%coeff(:, :, 1), wfn%density(:, :, 1)) */
    get_density_matrix(
      wfn.focc[0],
      wfn.coeff[0],
      wfn.density[0]
    );
    break;
  case 2:
    assert(false&&"Unimplemented");
  default:
    break;
  }
}

/*
subroutine get_electronic_energy(h0, density, energies)
   real(wp), intent(in) :: h0(:, :)
   real(wp), intent(in) :: density(:, :, :)
   real(wp), intent(inout) :: energies(:)

   integer :: iao, jao, spin

   ! $omp parallel do collapse(3) schedule(runtime) default(none) &
   ! $omp reduction(+:energies) shared(h0, density) private(spin, iao, jao)
   do spin = 1, size(density, 3)
      do iao = 1, size(density, 2)
         do jao = 1, size(density, 1)
            energies(iao) = energies(iao) + h0(jao, iao) * density(jao, iao, spin)
         end do
      end do
   end do
end subroutine get_electronic_energy
*/
__device__
void get_electronic_energy(
  const float (&h0)[MAX_NAO][MAX_NAO], 
  const float (&density)[MAX_NSPIN][MAX_NAO][MAX_NAO], 
  float (&energies)[MAX_NAO]
)
{
  for (int spin = 0; spin < MAX_NSPIN; ++spin) {
    for (int iao = 0; iao < MAX_NAO; ++iao) {
      for (int jao = 0; jao < MAX_NAO; ++jao) {
        energies[iao] += h0[iao][jao] * density[spin][iao][jao];
      }
    }
  }
}

/*
subroutine reduce(reduced, full, map)
   real(wp), intent(inout) :: reduced(:)
   real(wp), intent(in) :: full(:)
   integer, intent(in) :: map(:)

   integer :: ix

   do ix = 1, size(map)
      reduced(map(ix)) = reduced(map(ix)) + full(ix)
   end do
end subroutine reduce
*/
__device__
void reduce(
  float (&reduced)[MAX_NAT],
  const float (&full)[MAX_NAO],
  const int (&map)[MAX_NAO]
)
{
  for (int ix = 0; ix < MAX_NAO; ++ix) {
    int mapped_index = map[ix];
    if (mapped_index >= 0 && mapped_index < MAX_NAT) {
      reduced[mapped_index] += full[ix];
    }
  }
}

__device__
void next_scf(
  int &iscf,
  const structure_type &mol,
  const basis_type &bas,
  wavefunction_type wfn,
  sygvd_solver &solver,
  broyden_mixer &mixer,
  const scf_info &info,
  const tb_coulomb &coulomb,
  /*const dispersion_type &dispersion,*/
  /*const container_list &interactions,*/
  const integral_type &ints,
  potential_type &pot,
  coulomb_cache &cache, /* NOTE: Not an abstract container_cache since we don't mess with the polymorphism for now */
  /*container_cache &dcache,*/
  /*container_cache &icache,*/
  float (&energies)[MAX_NAT]
  /*error_type &error*/
) {
  // printf("NEXT SCF %i\n", iscf);
  float eao[MAX_NAO] = {0};
  float ts = 0;
  if (iscf > 0)
  {
    // next(mixer);
    // get_mixer(mixer, bas, wfn, info); // Update wfn from mixer, pretty much
    // assert(false && "Unimplemented");
  }
  pot.reset(); 
  iscf++;

  /* if (present(coulomb)) then */
    // call coulomb%get_potential(mol, cache, wfn, pot)
    coulomb.get_potential(mol, cache, wfn, pot);
  /*end if*/
  /*if (present(dispersion)) then
    call dispersion%get_potential(mol, dcache, wfn, pot) // Not yet done
  end if*/
  /*if (present(interactions)) then
    call interactions%get_potential(mol, icache, wfn, pot) // Not relevant for EIMS
  end if*/

  add_pot_to_h1(bas, ints, pot, wfn.coeff);

  set_mixer(mixer, wfn, info);

  get_density(wfn, solver, ints, ts /*,error*/);

  /* call get_mulliken_shell_charges(bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)*/
  get_mulliken_shell_charges(
    bas,
    ints.overlap,
    wfn.density,
    wfn.n0sh,
    wfn.qsh
  );

  /*    call get_qat_from_qsh(bas, wfn%qsh, wfn%qat) */
  get_qat_from_qsh(bas, wfn.qsh, wfn.qat);

  /*call get_mulliken_atomic_multipoles(bas, ints%dipole, wfn%density, &
      & wfn%dpat)*/
  get_mulliken_atomic_multipoles(
    bas, ints.dipole, wfn.density, wfn.dpat);

  // call get_mulliken_atomic_multipoles(bas, ints%quadrupole, wfn%density, &
  //     & wfn%qpat
  get_mulliken_atomic_multipoles(bas, 
    ints.quadrupole, 
    wfn.density, 
    wfn.qpat);

  diff_mixer(mixer, wfn, info);

  /* allocate(eao(bas%nao), source=0.0_wp) */
  for (int i = 0; i < bas.nao; i++)
  {
    eao[i] = 0;
  }
  
  /*call get_electronic_energy(ints%hamiltonian, wfn%density, eao)*/
  get_electronic_energy(ints.hamiltonian, wfn.density, eao);

  /* energies(:) = ts / size(energies) */
  for (int i = 0; i < mol.nat; i++)
  {
    energies[i] = ts / mol.nat;
  }
  
  /* call reduce(energies, eao, bas%ao2at) */
  reduce(energies, eao, bas.ao2at);

  /* if (present(coulomb)) then */
    /*call coulomb%get_energy(mol, cache, wfn, energies)*/
    coulomb.get_energy(mol, cache, wfn, energies);
  /* end if */
  /* if (present(dispersion)) then
      call dispersion%get_energy(mol, dcache, wfn, energies)
   end if */ // Dispersion is ignored for now
  /*
  if (present(interactions)) then
      call interactions%get_energy(mol, icache, wfn, energies)
   end if
  */ // Interactions are ignored for now
  
}