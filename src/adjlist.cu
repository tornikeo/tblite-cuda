#include "adjlist.h"

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
)
{
 /* 
     integer :: iat, jat, itr, img
     real(wp) :: r2, vec(3), cutoff2

     img = 0
     cutoff2 = cutoff**2

     call resize(nlat, 10*mol%nat)
     call resize(nltr, 10*mol%nat)

     do iat = 1, mol%nat
        inl(iat) = img
        do jat = 1, merge(mol%nat, iat, complete)
           do itr = 1, size(trans, 2)
              vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
              r2 = sum(vec**2)
              if (r2 < epsilon(cutoff2) .or. r2 > cutoff2) cycle
              img = img + 1
              if (size(nlat) < img) call resize(nlat)
              if (size(nltr) < img) call resize(nltr)
              nlat(img) = jat
              nltr(img) = itr
           end do
        end do
        nnl(iat) = img - inl(iat)
     end do

     call resize(nlat, img)
     call resize(nltr, img)
*/

 int iat = 0, jat = 0, itr = 0, img = 0;
 float r2 = 0, vec[3] = {0}, cutoff2 = 0;

 for (int iat = 0; iat < mol.nat; iat++)
 {
   inl[iat] = img;
   for (int jat = 0; jat < mol.nat; jat++)
   {
     // for (int itr = 0; itr < 1; itr++) // Assuming trans has a fixed size of 1 in the second dimension
     // {
     for (int k = 0; k < 3; k++)
     {
       vec[k] = mol.xyz[k][iat] - mol.xyz[k][jat]; // - trans[k][itr]
     }
     r2 = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
     cutoff2 = cutoff * cutoff;

     if (r2 < 1e-6f || r2 > cutoff2)
     {
       continue;
     }

     nlat[img] = jat;
     nltr[img] = itr;
     img++;
     // }
   }
   nnl[iat] = img - inl[iat];
 }
}

__device__
void new_adjacency_list(
  adjacency_list &self,
  const structure_type &mol,
  /* const float (&trans)[3][1], */
  const float &cutoff,
  const bool complete
)
{
  bool cmplt = complete;

  // Initialize inl array
  for (int i = 0; i < MAX_NAT; i++)
  {
    self.inl[i] = 0;
  }
  for (int i = 0; i < MAX_NAT; i++)
  {
    self.nnl[i] = 0;
  }

  // Generate the neighbor list
  assert(!cmplt);
  generate(mol, /* trans, */ cutoff, self.inl, self.nnl, self.nlat, self.nltr /*, cmplt*/);  
}
