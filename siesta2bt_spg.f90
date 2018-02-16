!
! Program to generate output for BoltzTrap
! 


subroutine siesta2bt(cell, xa, na)

 use spglib_f08

 use fdf
 use precision,    only : dp, sp
 use siesta_geom,  only : isa
 use sys,          only : die

 use chemical,     only : number_of_species
 use atomlist,     only : qtot

 ! Needed for band call
 use atomlist,        only : no_s, no_u, no_l, indxuo
 use m_spin,          only : h_spin_dim, spinor_dim
 use sparse_matrices, only : maxnh, listh, listhptr, numh
 use sparse_matrices, only : H, S, xijo
 use siesta_options,  only : occtol
 use m_energies,      only : ef
 use band,            only : bands

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
! inversion of lattice vectors for fractional coordinates is done
! through subroutine reclat.
! *******************************************************************

 integer                    :: na, na4
 real(dp)                   :: cell(3,3), cell2(3,3)
 real(dp)                   :: xa(3,na)
 real(dp)                   :: tcell(3,3)
 real(dp)                   :: fxa(3,na)
 real(dp), allocatable      :: fxa4(:,:)
 integer,  allocatable      :: isa4(:)
! integer                    :: na_conv

 ! irreducible k-points (ikp), number of k-points, num. of ikp

 real(dp),      allocatable :: ikp(:,:)
 integer                    :: nk
 integer                    :: ink
 integer                    :: maxnk

 ! spglib variables
 
 integer,      allocatable :: grid_address(:,:)
 integer,      allocatable :: map(:)
 integer                   :: mesh(3)
 integer                   :: is_shift(3)
 integer                   :: is_time_reversal
 real(dp)                  :: symprec
 type(SpglibDataset)       :: ds

 ! fdf variables

 type(block_fdf)            :: bfdf
 type(parsed_line), pointer :: pline

 ! internal variables

 integer                    :: i, ii, iii, counter_kp_ibz, nb2write
 logical                    :: bt_calc

 ! reciprocral cell and weights

 integer                    :: bt_rcell(3,3)
 real(dp)                   :: bt_wrcell(3)

 ! module and angule of lattice vectors

 real(dp)                   :: cellm(3), celang(3)
 real(dp),        parameter ::  pi = 3.1415926d0

 ! Band energies
 real(dp),    allocatable :: ebk(:,:,:)

 ! transposing cell for c-fortran wrapper used in spglib
 tcell = transpose(cell)

 ! check if block BT.kgrid_Monkhorst_Pack
 bt_calc = fdf_block('BT.kgrid_Monkhorst_Pack',bfdf)

 if ( .not. bt_calc ) then
  return
 endif
 
 ! read block BT.kgrid_Monkhorst_Pack
 do i=1, 3
   if (.not. fdf_bline(bfdf,pline)) then
     call die('siesta2bt: ERROR in BT.kgrid_Monkhorst_Pack block')
   end if
     bt_rcell(1,i) = fdf_bintegers(pline,1)
     bt_rcell(2,i) = fdf_bintegers(pline,2)
     bt_rcell(3,i) = fdf_bintegers(pline,3)
   if ( fdf_bnvalues(pline) > 3 ) then
     bt_wrcell(i) = fdf_bvalues(pline,4)
   else
     bt_wrcell(i) = 0._dp
   end if
 enddo
 maxnk = abs( bt_rcell(1,1) * bt_rcell(2,2) * bt_rcell(3,3) + &
              bt_rcell(2,1) * bt_rcell(3,2) * bt_rcell(1,3) + &
              bt_rcell(3,1) * bt_rcell(1,2) * bt_rcell(2,3) - &
              bt_rcell(1,1) * bt_rcell(3,2) * bt_rcell(2,3) - &
              bt_rcell(2,1) * bt_rcell(1,2) * bt_rcell(3,3) - &
              bt_rcell(3,1) * bt_rcell(2,2) * bt_rcell(1,3) )   

 allocate( ebk(no_u, spinor_dim, maxnk) )

 is_shift = bt_wrcell

 ! read symprec from fdf file
 symprec = fdf_single('SPG.Symprec',1.e-2)
 
 call cart2frac(na, xa(1,:), xa(2,:), xa(3,:), cell, fxa)

 ! TODO: symmetrize structure
 ! spg_find_primitive(cell, fxa, isa, na, symprec) => ?
 cell2 = tcell

 allocate(fxa4(3,na*4), isa4(na*4), grid_address(3,maxnk), map(maxnk))

 fxa4 = 0.0
 isa4 = 0
 fxa4(:,1:na) = fxa
 isa4(1:na) = isa

 is_time_reversal = 0

 ! TODO: calculate rotation matrices, translation vectors
 ! and lattice type.
 ! ds usage:
 ! https://atztogo.github.io/spglib/api.html?highlight=kpoint#spg-get-dataset-and-spg-get-dataset-with-hall-number
 ds = spg_get_dataset(tcell, fxa, isa, na, symprec)

 print*, 'siesta2bt:  Space Group:', ds%spacegroup_number

 print*, 'siesta2bt:  Number of symops:', ds%n_operations
 ! get irreducible k-points (ik-points)
 ! What do in the case of non-cartesian basis?
 mesh(1) = bt_rcell(1,1) 
 mesh(2) = bt_rcell(2,2) 
 mesh(3) = bt_rcell(3,3) 

 nk = maxnk

 ink = spg_get_ir_reciprocal_mesh(grid_address, &
                                  map,          &
                                  mesh,         &
                                  is_shift,     &
                                  is_time_reversal,  &
                                  tcell,        &
                                  fxa,          &
                                  isa,          &
                                  na,           &
                                  symprec)

 allocate(ikp(3,ink))
 ikp=0.d0
!
! To account for the difference of indices in C and Fortran, add 1 to values of all
! elements of map
!
 map = map + 1

!
! Extracts k-points in the ibz. map contains the index of a point in
! grid_address; if iteration variable i and map(i) are equal, then this is one
! point in the IBZ.
!
 
 counter_kp_ibz = 0
 do i = 1, maxnk
    if (i == map(i)) then
       counter_kp_ibz = counter_kp_ibz + 1
       ikp(:,counter_kp_ibz) = dble(is_shift + 2*grid_address(:,map(i)))/dble(2*mesh)
!       print*, ikp(:,counter_kp_ibz)
    endif
 end do

 print*, ' siesta2bt: Number of k-points in IBZ, from spglib and extracted from grid_address: ',ink, counter_kp_ibz

 ! calculate eigenenergies

 call bands( no_s, h_spin_dim, spinor_dim, no_u, no_l, maxnh, maxnk, &
             numh, listhptr, listh, H, S, ef, xijo, indxuo, &
             .true., ink, ikp, ebk, occtol, .false. )

 ! Write BoltzTrap input files:
 ! ===============================================
 ! (1) intrans file
 
 open(101,file=trim(slabel)//'.intrans')
 
 write(101,'(A)') 'GENE                      # Format of DOS' 
 write(101,'(A)') '0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap'
 write(101,'(F7.5,A,F6.0,A)') Ef, ' 0.0005, 0.4 ',qtot,&
            ' # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons'
 write(101,'(A)') 'CALC                      # CALC (calculate expansion coeff), NOCALC read from file'
 write(101,'(A)') '5                         # lpfac, number of latt-points per k-point'
 write(101,'(A)') 'BOLTZ                     # run mode (only BOLTZ is supported)'
 write(101,'(A)') '0.15                      # (efcut) energy range of chemical potential'
 write(101,'(A)') '800. 50.                  # Tmax, temperature grid'
 write(101,'(A)') '-1                        # energyrange of bands given individual DOS output'
 write(101,'(A)') '                          # sig_xxx and dos_xxx (xxx is band number)'
 write(101,'(A)') 'TETRA'

 close(101) 

 ! ===============================================
 ! (2) struct file
 open(102,file=trim(slabel)//'.struct')

 ! title 
 write(102,'(A)') trim(slabel)
 
 ! lattice vectors (Bohr)
 write(102,'(3F16.9)') cell(1,:)
 write(102,'(3F16.9)') cell(2,:)
 write(102,'(3F16.9)') cell(3,:)

 ! number of operations
 write(102,'(I2)') ds%n_operations

 ! rotation matrices are transposed because C ordering and Fortran ordering are
 ! different
 do i=1, ds%n_operations
   ds%rotations(:,:,i) = transpose(ds%rotations(:,:,i))
   write(102,'(3I4)') ds%rotations(1,:,i)
   write(102,'(3I4)') ds%rotations(2,:,i)
   write(102,'(3I4)') ds%rotations(3,:,i)
   write(102,'(I4)') i
 enddo

 close(102)

 ! ===============================================
 ! (3) energy file. k-points are written **in fracional coordinates**

 open(103,file=trim(slabel)//'.energy')
 
 ! title 
 write(103,*) trim(slabel)

 ! number of irreducible kpoints
 write(103,'(I8)') ink
 
 nb2write = 20  ! Trocar por NumberOfEigenstates para economizar poder computacional
 ! write k-points
 do i = 1, ink
   write(103,'(3F14.9,I5)') ikp(:,i), nb2write
   ! write energies
   do ii = 1,spinor_dim
     do iii = 1, nb2write
       write(103,'(F14.10)') ebk(iii,ii,i)
     enddo
   enddo
 enddo

 close(103)

 ! ===============================================
 ! (4) debug file

 open(104,file=trim(slabel)//'.debug')
 write(104,*) 'cell:'
 write(104,*) cell(:,1)*0.529177
 write(104,*) cell(:,2)*0.529177
 write(104,*) cell(:,3)*0.529177
 write(104,*) 'tcell:'
 write(104,*) tcell(:,1)*0.529177
 write(104,*) tcell(:,2)*0.529177
 write(104,*) tcell(:,3)*0.529177
 write(104,*) na
 do i=1, na
   write(104,*) fxa(:,i)
 enddo
 write(104,*) 'space group num:',ds%spacegroup_number
 write(104,*) 'cell2:'
 write(104,*) cell2(1,:)
 write(104,*) cell2(2,:)
 write(104,*) cell2(3,:)
 close(104)

 
end subroutine siesta2bt

