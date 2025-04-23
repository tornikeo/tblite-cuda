#include <cuda_runtime.h>
#include <stdio.h>
#include <cassert>
#include <cstdlib>
#include "xtb.h"
#include "basis/type.h"
#include "potential.h"
#include "adjlist.h"
#include "wavefunction/type.h"
#include "integral/type.h"
#include "scf/broyden.h"
#include "scf/iterators.h"
#include "lapack/sygvd.h"

template<int N>
__device__ inline float sum_square(const float (&arr)[N]) {
  float val = 0;
  for (int i = 0; i < N; i++)
  { 
    val += arr[i] * arr[i];
  }
  return val;
}

template<int N>
__device__ inline float sum(const float (&arr)[N]) {
  float val = 0;
  for (int i = 0; i < N; i++)
  { 
    val += arr[i];
  }
  return val;
}

/*   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Exponents for all element pairs
   real(wp), intent(in) :: alpha(:, :)
   !> Effective nuclear charges for all element pairs
   real(wp), intent(in) :: zeff(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: kexp(:, :)
   !> Pairwise parameters for all element pairs
   real(wp), intent(in) :: rexp(:, :)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), intent(inout) :: sigma(:, :)
*/

__device__ 
inline void get_repulsion_derivs(
  const structure_type &mol,
  const float cutoff,
  const float (&alpha)[MAX_NELEM][MAX_NELEM],
  const float (&zeff)[MAX_NELEM][MAX_NELEM],
  const float (&kexp)[MAX_NELEM][MAX_NELEM],
  const float (&rexp)[MAX_NELEM][MAX_NELEM],
  float (&energies)[MAX_NAT],
  float (&gradient)[MAX_NAT][3],
  float (&sigma)[3][3]
) {
  int iat, jat, izp, jzp, itr;
  float r1, r2, rij[3], r1k, r1r, exa, cutoff3, dE, dG[3], dS[3][3];
  float cutoff2 = cutoff*cutoff;
  /*
     do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            r1k = r1**kexp(jzp, izp)
            exa = exp(-alpha(jzp, izp)*r1k)
            r1r = r1**rexp(jzp, izp)
            dE = zeff(jzp, izp) * exa/r1r
            dG = -(alpha(jzp, izp)*r1k*kexp(jzp, izp) + rexp(jzp, izp)) * dE * rij/r2
            dS = spread(dG, 1, 3) * spread(rij, 2, 3)
            energies(iat) = energies(iat) + 0.5_wp * dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + 0.5_wp * dE
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma + dS
            else
               sigma(:, :) = sigma + 0.5_wp * dS
            end if
         end do
      end do
   end do
  */
  for (int iat = 0; iat < MAX_NAT; iat++)
  {
    izp = mol.id[iat];
    for (int jat = 0; jat < iat; jat++)
    {
      jzp = mol.id[jat];
      for (int i = 0; i < 3; i++)
      {
        rij[i] = mol.xyz[iat][i] - mol.xyz[jat][i]; // trans is assumed to be zero
      }

      r2 = sum_square(rij);
      if(r2 > cutoff2 || r2 < 1.0e-12) continue;
      r1 = sqrt(r2);

      // r1k = r1**kexp(jzp, izp)
      // exa = exp(-alpha(jzp, izp)*r1k)
      // r1r = r1**rexp(jzp, izp)
      // dE = zeff(jzp, izp) * exa/r1r

      r1k = powf(r1, kexp[izp][jzp]);
      exa = expf(-alpha[izp][jzp] * r1k);
      r1r = powf(r1, rexp[izp][jzp]);
      dE = zeff[izp][jzp] * exa / r1r;

      for (int i = 0; i < 3; i++) {
        dG[i] = -(alpha[izp][jzp] * r1k * kexp[izp][jzp] + rexp[izp][jzp]) * dE * rij[i] / r2;
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          dS[i][j] = dG[i] * rij[j];
        }
      }

      energies[iat] += 0.5f * dE;
      
      // dG = -(alpha(jzp, izp)*r1k*kexp(jzp, izp) + rexp(jzp, izp)) * dE * rij/r2
      // dS = spread(dG, 1, 3) * spread(rij, 2, 3)
      // if (iat /= jat) then
      //     energies(jat) = energies(jat) + 0.5_wp * dE
      //     gradient(:, iat) = gradient(:, iat) + dG
      //     gradient(:, jat) = gradient(:, jat) - dG
      //     sigma(:, :) = sigma + dS
      // else
      //     sigma(:, :) = sigma + 0.5_wp * dS
      // end if
      if (iat != jat) {
        energies[jat] += 0.5f * dE;
        for (int i = 0; i < 3; i++) {
          gradient[iat][i] += dG[i];
          gradient[jat][i] -= dG[i];
        }
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            sigma[i][j] += dS[i][j];
          }
        }
      } else {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            sigma[i][j] += 0.5f * dS[i][j];
          }
        }
      }
    }
  }
}

__device__ void get_engrad(
    const tb_repulsion &self,
    const structure_type &mol,
    float (&energies)[MAX_NAT],
    float (&graident)[MAX_NAT][3],
    float (&sigma)[3][3])
{
  // call get_repulsion_derivs(mol, trans, self%cutoff, self%alpha, self%zeff, &
  //   & self%kexp, self%rexp, energies, gradient, sigma)
  get_repulsion_derivs(
    mol,
    self.cutoff,
    self.alpha,
    self.zeff,
    self.kexp,
    self.rexp,
    energies,
    graident,
    sigma
  );
}


// __device__
// void get_amat_0d()

__device__ void xtb_singlepoint(
    const structure_type &mol,
    const xtb_calculator &calc,
    wavefunction_type &wfn,
    const float accuracy,

    float energy,
    float (&gradient)[MAX_NAT][3],
    float (&sigma)[3][3],
    const int verbosity)
{
  bool grad, converged, econverged, pconverged;
  float econv, pconv, cutoff, elast, dpmom[3], qpmom[6], nel;
  int iscf;
  /* 
   if (present(verbosity)) then
      prlevel = verbosity
   else
      prlevel = ctx%verbosity
   end if
  */

  int prlevel = verbosity;
  // allocate(energies(mol%nat), source=0.0_wp)
  // allocate(erep(mol%nat), source=0.0_wp)
  // allocate(edisp(mol%nat), source=0.0_wp)
  // allocate(eint(mol%nat), source=0.0_wp)
  // allocate(exbond(mol%nat), source=0.0_wp)
  // allocate(cn(mol%nat))
  // allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
  // allocate(selfenergy(calc%bas%nsh), dsedcn(calc%bas%nsh))
  // allocate(eelec(mol%nat), source=0.0_wp)
  // allocate(dEdcn(mol%nat))
  // allocate(wdensity(calc%bas%nao, calc%bas%nao, wfn%nspin))

  float energies[MAX_NAT] = {0};
  float erep[MAX_NAT] = {0};
  float edisp[MAX_NAT] = {0};
  float eint[MAX_NAT] = {0};
  float exbond[MAX_NAT] = {0};
  float cn[MAX_NAT] = {0};
  float dcndr[MAX_NAT][MAX_NAT][3] = {0};
  float dcndL[MAX_NAT][3][3] = {0};
  float selfenergy[MAX_NSH] = {0};
  float dsedcn[MAX_NSH] = {0};
  float eelec[MAX_NAT] = {0};
  float dEdcn[MAX_NAT] = {0};
  float wdensity[MAX_NSPIN][MAX_NAO][MAX_NAO] = {0};

  adjacency_list list;
  potential_type pot;
  coulomb_cache ccache;
  integral_type ints;
  broyden_mixer mixer;
  sygvd_solver sygvd;

  // gradient(:, :) = 0.0_wp
  // sigma(:, :) = 0.0_wp
  for (size_t i = 0; i < MAX_NAT; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      gradient[i][j] = 0;
    }
  }

  for (size_t i = 0; i < 3; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      sigma[i][j] = 0;
    }
  }

  /* if (allocated(calc%repulsion)) then */
    get_engrad(calc.repulsion, mol, erep, gradient, sigma);
    for (int i = 0; i < MAX_NAT; i++)
    {
      energies[i] += erep[i];
    }
    if(prlevel > 1)
    {
      printf("repulsion energy %d Eh\n", sum(erep));
    }  
  /* end if */

  /* if (allocated(calc%dispersion)) then  */
    /* TODO: priority medium 

  dispersion_get_engrad(calc.dispersion, mol, edisp, gradient, sigma); */ 

  /* if (allocated(calc%interactions)) then TODO */
    /* TODO: priority low */

  /* call new_potential(pot, mol, calc%bas, wfn%nspin) */
  // new_potential(pot, mol, calc.bas, wfn.nspin);
  new_potential(pot, mol, calc.bas, 1);
  /* if (allocated(calc%coulomb)) then */
  /* call calc%coulomb%update(mol, ccache) */
  update(calc.coulomb, mol, ccache);

  /* call get_occupation(mol, calc%bas, calc%h0, wfn%nocc, wfn%n0at, wfn%n0sh) */
  get_occupation(mol, calc.bas, calc.h0, wfn.nocc, wfn.n0at, wfn.n0sh);

  /* nel = sum(wfn%n0at) - mol%charge
   if (mod(mol%uhf, 2) == mod(nint(nel), 2)) then
      wfn%nuhf = mol%uhf
   else
      wfn%nuhf = mod(nint(nel), 2)
   end if
   */

  nel = 0.0f;
  for (int i = 0; i < MAX_NAT; i++) {
    nel += wfn.n0at[i];
  }
  nel -= mol.charge;

  if (static_cast<int>(mol.uhf) % 2 == static_cast<int>(roundf(nel)) % 2) {
    wfn.nuhf = mol.uhf;
  } else {
    wfn.nuhf = static_cast<int>(roundf(nel)) % 2;
  }

  /* call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2)) */
  get_alpha_beta_occupation(wfn.nocc, wfn.nuhf, wfn.nel[0], wfn.nel[1]);
  /*    if (prlevel > 1) print *, property("number of electrons", wfn%nocc, "e") */

  if(prlevel > 1)
  {
    printf("Number of electrons: %f e\n", wfn.nocc);
  }

  /* call timer%push("hamiltonian")
   if (allocated(calc%ncoord)) then
      allocate(cn(mol%nat))
      if (grad) then
         allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      end if
      call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
   end if
  */
  ncoord_get_cn(calc.ncoord, mol, cn, dcndr, dcndL);

  /* allocate(selfenergy(calc%bas%nsh), dsedcn(calc%bas%nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=selfenergy, dsedcn=dsedcn)
  */  
  get_selfenergy(calc.h0, mol.id, calc.bas.ish_at, calc.bas.nsh_id, 
  /*cn=*/cn, /*selfenergy=*/selfenergy, /*dsedcn=*/dsedcn);

  /* cutoff = get_cutoff(calc%bas, accuracy) */
  cutoff = get_cutoff(calc.bas, accuracy);

  /* TODO: Priority low. lattices are unsupported */
  /* call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr) */
  // float lattr[3][1] = {0};
  /* call new_adjacency_list(list, mol, lattr, cutoff) */
  new_adjacency_list(list, mol, /*lattr,*/ cutoff, false);

  /*    if (prlevel > 1) then
      print *, property("integral cutoff", cutoff, "bohr")
      print *
   end if 
   */
  if(prlevel > 1)
  {
    printf("Integral cutoff: %f bohr\n", cutoff);
  }

  /* call new_integral(ints, calc%bas%nao) */  
  new_integral(ints, calc.bas.nao);

  /* call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian) */
  get_hamiltonian(
    mol, /* lattr, */list, calc.bas, calc.h0, selfenergy, 
    ints.overlap, ints.dipole, ints.quadrupole, ints.hamiltonian
  );

  /*allocate(eelec(mol%nat), source=0.0_wp)*/
  iscf = 0;
  converged = false;
  /*info = calc%variable_info()*/
  /*call new_broyden(mixer, calc%max_iter, wfn%nspin*get_mixer_dimension(mol, calc%bas, info), &
      & calc%mixer_damping)*/
  constexpr scf_info info {
    /*charge = */shell_resolved,
    /*dipole = */atom_resolved,
    /*quadrupole = */atom_resolved
  };
  int broyden_ndim = wfn.nspin * get_mixer_dimension(mol, calc.bas, info);
  printf("Broyden dimension %i\n", broyden_ndim);
  new_broyden(
    mixer, 
    calc.max_iter, 
    broyden_ndim,
    calc.mixer_damping
  );

  // for (int i = 0; i < BROYDEN_NDIM; i++)
  // {
  //   for (int j = 0; j < MAX_ITER_DEFAULT; j++)
  //   {
  //     printf("%i,", mixer.df[i][j]);
  //   }
  //   printf("\n");
  // }

  /*   if (prlevel > 0) then
      call ctx%message(repeat("-", 60))
      call ctx%message("  cycle        total energy    energy error   density error")
      call ctx%message(repeat("-", 60))
   end if
  */
 if(prlevel > 1)
 {
  printf("----------------------------------------------------------------\n");
  printf("  cycle        total energy    energy error   density error\n");
  printf("----------------------------------------------------------------\n");
 }
 while(!converged && iscf < calc.max_iter)
 {
  /* call next_scf(iscf, mol, calc%bas, wfn, sygvd, mixer, info, &
         & calc%coulomb, calc%dispersion, calc%interactions, ints, pot, &
         & ccache, dcache, icache, eelec, error)
  */
  // next_scf(
  //   iscf, mol, calc.bas, wfn, sygvd, mixer, info,
  //   calc.coulomb, /*calc.dispersion,*/ /*calc.interactions,*/ ints,
  //   pot, ccache, /*icache,*/ eelec /*error*/
  // );

  if (prlevel > 0)
  {
    printf("%7d %24.13f %16.7e \n", /*%16.7e*/
         iscf, 
         sum(eelec) + sum(energies), 
         sum(eelec) - elast
        /*mixer.get_error()*/);
  }
 
  iscf++;
 }
}

__global__ void test_xtb_singlepoint()
{
  structure_type mol = {0};
  xtb_calculator calc = {0};
  wavefunction_type wfn = {0};
  /*
  -exec print mol%xyz
  ((0.98819706000000118, -0.79388594999999362, 1.8932948000000023) (-1.1710541999999968, 0.074677960000005192, 3.3566230999999975) (1.6169206000000027, -3.8529364000000079, 3.2129892000000067) (-1.2120191999999921, 2.6964877999999981, 3.6464740000000031) (-1.6371701000000018, 3.4899488999999972, 2.1451410999999907) (-2.9580292000000052, -0.57318774000000117, 2.5348827999999992) (-0.97500824000000375, -0.79775967000000325, 5.2966186000000048) (0.40421389000000885, -1.4110043000000081, -0.10239630999999662) (2.7573567000000088, 0.39135772999999358, 1.9942773999999945))
  */
  int id[MAX_NAT] = {1, 1, 2, 3, 4, 4, 4, 4, 4};
  float xyz[MAX_NAT][3] = {
    {0.98819706000000118, -0.79388594999999362, 1.8932948000000023},
    {-1.1710541999999968, 0.074677960000005192, 3.3566230999999975},
    {1.6169206000000027, -3.8529364000000079, 3.2129892000000067},
    {-1.2120191999999921, 2.6964877999999981, 3.6464740000000031},
    {-1.6371701000000018, 3.4899488999999972, 2.1451410999999907},
    {-2.9580292000000052, -0.57318774000000117, 2.5348827999999992},
    {-0.97500824000000375, -0.79775967000000325, 5.2966186000000048},
    {0.40421389000000885, -1.4110043000000081, -0.10239630999999662},
    {2.7573567000000088, 0.39135772999999358, 1.9942773999999945}
  };

  for (int i = 0; i < MAX_NAT; i++) {
    for (int j = 0; j < 3; j++) {
      mol.xyz[i][j] = xyz[i][j];
    }
  }

  for (int i = 0; i < MAX_NAT; i++) {
    mol.id[i] = id[i];
  }

  mol.nat = MAX_NAT;
  calc.bas.nsh = MAX_NSH;
  calc.max_iter = MAX_ITER_DEFAULT;

  float accuracy = 0.01f;
  float energy = 0.0f;

  float gradient[MAX_NAT][3] = {0};
  float sigma[3][3] = {0.0};
  int verbosity = 2;

  xtb_singlepoint(mol, calc, wfn, accuracy, energy, gradient, sigma, verbosity);
  // printf("Energy: %f\n", energy);
}

void xtb_test()
{
  cudaDeviceSetLimit(cudaLimitStackSize, 1024 * sizeof(float));
  cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128 * 1024 * 1024);
  test_xtb_singlepoint<<<1, 1>>>();
  cudaDeviceSynchronize();
}