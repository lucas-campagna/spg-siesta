!
! Program to generate output for BoltzTrap
! 


subroutine siesta2bt(cell, xa, na)

 use spglib_f08

 use m_find_kgrid,   only : find_kgrid

 use fdf
 use precision,    only : dp, sp
 use siesta_geom,  only : isa
 use sys,          only : die

 use chemical,     only : number_of_species
 use atomlist,     only : qtot

 ! Needed for band call
 use atomlist,            only : no_s, no_u, no_l, indxuo
 use m_spin,              only : h_spin_dim, spinor_dim
 use sparse_matrices,     only : maxnh, listh, listhptr, numh
 use sparse_matrices,     only : H, S, xijo
 use siesta_options,      only : occtol
 use m_energies,          only : ef
 use band,                only : bands
 use m_get_kpoints_scale, only: get_kpoints_scale


 use files, only : slabel

 implicit none

! *******************************************************************
! This subroutine prepare input files for BoltzTraP simulations generating
! the following files:
! <slabel>.intras
! <slabel>.struct
! <slabel>.enegy
!
! where <slabel> will be replaced by the SystemLabel.
! The .intrans file may require some aditional adjusts according to your study purpose.
!
! Written by Lucas Prett Campagna under the supervision of
! Rodrigo Garcia Amorim and Marcos Verissimo Alves, January 2018.
! ********* INPUT ***************************************************
! integer na                : number of atoms
! double precision cell(3,3): Lattice (supercell) vectors
! double precision xa(3,na) : atomic coordinates in Bohr cartesian
! integer isa(na)           : atomic species of different atoms
! *******************************************************************

 integer                    :: na
 real(dp)                   :: cell(3,3)
 real(dp)                   :: xa(3,na)
 real(dp)                   :: tcell(3,3)
 real(dp)                   :: fxa(3,na)

 ! reciprocal cell, irreducible k-points, number of (i)k-points, weight of k-points, index of bands to be saved

 real(dp)                   :: rcell(3,3)
 real(dp),      allocatable :: kpoint_frac(:,:)
 real(dp),      allocatable :: kpoint_cart(:,:)
 integer                    :: nk, ink
 integer                    :: first_band, last_band
 real(dp),      allocatable :: kweight(:)

 ! spglib variables
 
 logical                   :: time_reversal
 real(dp)                  :: symprec
 type(SpglibDataset)       :: ds

 ! fdf variables

 type(block_fdf)            :: bfdf
 type(parsed_line), pointer :: pline

 ! internal variables

 integer                    :: i, ii, iii
 logical                    :: bt_calc

 !mesh and shift of reciprocal space

 integer                    :: mesh(3)
 integer                    :: shift(3)

 ! Band energies
 real(dp),    allocatable :: ebk(:,:,:)

 ! transposing cell for c-fortran wrapper used in spglib
 tcell = transpose(cell)

 ! check if block BT.kgrid_Monkhorst_Pack
 bt_calc = fdf_block('BT.kgrid_Monkhorst_Pack',bfdf)

 if ( .not. bt_calc ) then
  return
 endif
 
 ! ===================================================================
 ! (1) Get input parameters:
 ! 1) mesh and shift of recipricral space
 ! 2) symprec
 ! 3) First band index to save
 ! 4) Last band index to save
 ! 5) Time reversal
 ! ===================================================================

 ! read block BT.kgrid_Monkhorst_Pack
 do i=1, 3
   if (.not. fdf_bline(bfdf,pline)) then
     call die('siesta2bt: ERROR in BT.kgrid_Monkhorst_Pack block')
   end if
     mesh(i) = fdf_bintegers(pline,i)
   if ( fdf_bnvalues(pline) > 3 ) then
     shift(i) = nint(fdf_bvalues(pline,4))
   else
     shift(i) = 0
   end if
 enddo

 symprec = fdf_single('BT.Symprec',1.e-2)

 first_band = fdf_integer('BT.First_Band',1)
 last_band  = fdf_integer('BT.Last_Band',no_u)

 if (first_band > last_band) call die('BT.Last_Band argument must be bigger than BT.First_Band')
 if (last_band > no_u) last_band = no_u

 time_reversal  = fdf_integer('BT.Time_Reversal',1)

 ! ===================================================================
 ! (2) Find rotation matrices
 ! ===================================================================

 nk = product(mesh)

 call cart2frac(na, xa(1,:), xa(2,:), xa(3,:), cell, fxa)

 ! ds usage:
 ! https://atztogo.github.io/spglib/api.html?highlight=kpoint#spg-get-dataset-and-spg-get-dataset-with-hall-number
 ds = spg_get_dataset(tcell, fxa, isa(1:na), na, symprec)

 ! ===================================================================
 ! (3) Find irreducible k-points
 ! ===================================================================

 ink = nk
 
 allocate(kpoint_frac(3,nk), kweight(nk))

 call find_irreducible_kpoints(mesh,             &
                               shift,            &
                               time_reversal,    &
                               ds%rotations,     &
                               ds%n_operations,  &
                               kpoint_frac,      &
                               kweight,          &
                               ink)

 allocate( kpoint_cart(3,ink) )

 if(nk /= ink) then
   kpoint_cart = kpoint_frac
   deallocate(kpoint_frac)
   allocate(kpoint_frac(3,ink))
   kpoint_frac = kpoint_cart
 end if

 call reclat( cell, rcell, 1 )
 kpoint_cart = matmul(rcell,kpoint_frac)

 ! ===================================================================
 ! (4) Calculate bands
 ! ===================================================================

 allocate( ebk(last_band, spinor_dim, nk) )

 ! calculate eigenenergies
 call bands( no_s, h_spin_dim, spinor_dim, no_u, no_l, maxnh, nk, &
             numh, listhptr, listh, H, S, ef, xijo, indxuo, &
             .false., ink, kpoint_cart, ebk, occtol, .false. )

 ! ===================================================================
 ! (5) Save files in BoltzTrap format:
 ! ===================================================================
 ! 1) intrans file 
 
 open(101,file=trim(slabel)//'.intrans')
 
 write(101,'(A)') 'GENE                      # Format of DOS' 
 write(101,'(A)') '0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap'
 write(101,'(F7.5,A,F6.0,A)') Ef, ' 0.0005, 0.4 ',qtot,&
            ' # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons'
 write(101,'(A)') 'CALC                      # CALC (calculate expansion coeff), NOCALC read from file'
 write(101,'(A)') '5                         # lpfac, number of latt-points per k-point'
 write(101,'(A)') 'BOLTZ                     # run mode (only BOLTZ is supported)'
 write(101,'(A)') '0.15                      # (efcut) energy range of chemical potential'
 write(101,'(A)') '300. 10.                  # Tmax, temperature grid'
 write(101,'(A)') '-1                        # energyrange of bands given individual DOS output'
 write(101,'(A)') 'HISTO                     # sig_xxx and dos_xxx (xxx is band number)'

 close(101) 

 ! ===================================================================
 ! 2) struct file
 open(102,file=trim(slabel)//'.struct')

 ! title 
 write(102,'(A)') trim(slabel)
 
 ! lattice vectors (Bohr)
 write(102,'(3F16.9)') cell(:,1)
 write(102,'(3F16.9)') cell(:,2)
 write(102,'(3F16.9)') cell(:,3)

 ! number of operations
 write(102,'(I3)') ds%n_operations

 ! rotation matrices
 do i=1, ds%n_operations
   write(102,'(9I3)') ds%rotations(:,:,i)
 enddo

 close(102)

 ! ===================================================================
 ! 3) energy file

 open(103,file=trim(slabel)//'.energy')
 
 ! title 
 write(103,*) trim(slabel)

 ! number of irreducible kpoints
 write(103,'(I8)') ink
 
 ! write ik-points
 do i = 1, ink
   write(103,'(3F14.9,I5)') kpoint_frac(:,i), spinor_dim*(last_band - first_band) + 1
   
   ! write energies
   do ii = 1, spinor_dim
     do iii = first_band, last_band
       write(103,'(F14.10)') ebk(iii,ii,i)
     enddo
   enddo
 enddo

 close(103)

end subroutine siesta2bt

subroutine find_irreducible_kpoints(mesh, shift, is_time_reversal,rot, nrot, kpoint, kweight, ink)
 use precision,    only : dp
 implicit none
 ! ***************** INPUT **********************************************
 ! integer  mesh(3)    : Division of reciprocal space
 ! integer  shift(3)   : Supercell reciprocal of k-grid unit cell
 ! integer  nrot       : Numer of rotation matrixes
 ! integer  rot        : Rotations matrixes
 ! ***************** IN/OUTPUT *********************************************
 ! integer  ink        : In:  Number of total k-points
 !                     : Out: Number of irreducible k-points 
 ! ***************** OUTPUT *********************************************
 ! real(dp) kpoints(3,ink) : K-points in terms of reciprocal vectors
 ! real(dp) kweight(ink)   : weights
 !
 ! Input
 integer,   intent(in) :: mesh(3)
 integer,   intent(in) :: shift(3)
 integer,   intent(in) :: nrot      
 integer,   intent(in) :: rot(3,3,nrot) 
 integer,   intent(in) :: is_time_reversal
 !
 ! Output
 integer,  intent(inout) :: ink
 real(dp), intent(out) :: kpoint(3,ink)
 real(dp), intent(out) :: kweight(ink)
 ! 
 ! Local
 real(dp)              :: kpoint_rot(3)
 real(dp)              :: ds2
 integer               :: equiv(ink)
 integer               :: kindex_rot(3)
 integer               :: i, j, k, l, m, n
 integer               :: nk
 real(dp)              :: rr(3)
 logical               :: in_the_list
 !
 ds2 = 1.0E-8
 !
 nk = product(mesh)
 !
 ! TODO: check invertion matrix
 !
 do i=1, mesh(1)
   do j=1, mesh(2)
     do k=1, mesh(3)
       n = (k-1) + (j-1)*mesh(3) + (i-1)*mesh(2)*mesh(3) + 1
       kpoint(:,n) = ( 2*dble((/i,j,k/)-1) + dble(shift) )/(2*mesh)
     enddo 
   enddo 
 enddo 
 do n=1,nk
   equiv(n) = n
 enddo
 kweight = 1.d0
 do l=1, nk
   if (equiv(l) == l) then
     do m = 1, nrot
       kpoint_rot = matmul(rot(:,:,m),kpoint(:,l))
       kpoint_rot = kpoint_rot - nint(kpoint_rot)
       rr = kpoint_rot*mesh - 0.5d0*shift
       in_the_list = maxval(dabs(rr - nint(rr))/2.d0) <= ds2
       if (in_the_list) then
         kindex_rot = mod ( nint ( kpoint_rot*mesh - 0.5d0*shift + 2*mesh ) , mesh ) + 1
         n = (kindex_rot(3)-1) + (kindex_rot(2)-1)*mesh(3) + (kindex_rot(1)-1)*mesh(2)*mesh(3) + 1
         if (n > l .and. equiv(n) == n) then
           equiv(n) = l
           kweight(l) = kweight(l) + 1.0d0
         else
           IF (equiv(n)/=l .or. n<l ) write(*,*) 'siesta2bt: Erro kpoint_grid!'
         endif
         if (is_time_reversal == 1) then
           rr = -kpoint_rot*mesh - 0.5d0*shift
           in_the_list = maxval(dabs(rr - nint(rr))/2.d0) <= ds2
           if (in_the_list) then
             kindex_rot = mod ( nint ( - kpoint_rot*mesh - 0.5d0*shift + 2*mesh ) , mesh ) + 1
             n = (kindex_rot(3)-1) + (kindex_rot(2)-1)*mesh(3) + (kindex_rot(1)-1)*mesh(2)*mesh(3) + 1
             if (n > l .and. equiv(n) == n) then
               equiv(n) = l
               kweight(l) = kweight(l) + 1.0d0
             else
               IF (equiv(n)/=l .or. n<l ) write(*,*) 'siesta2bt: Erro kpoint_grid in is time reversal!'
             endif
           endif
         endif
       endif
     enddo
   endif
 enddo

 ink = 0
 do n = 1, nk
   if (equiv(n) == n) then
     ink = ink + 1
     kweight(ink) = kweight(n)
     kpoint(:,ink) = kpoint(:,n) - nint(kpoint(:,n))
   endif
 enddo
 kweight(1:ink) = 2*kweight(1:ink)/sum(kweight(1:ink))

end subroutine find_irreducible_kpoints
