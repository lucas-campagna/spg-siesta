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
 use atomlist,        only : no_s, no_u ,no_l, indxuo
 use m_spin,          only : h_spin_dim, spinor_dim
 use sparse_matrices, only : maxnh, listh, listhptr, numh
 use sparse_matrices, only : H, S, xijo
 use siesta_options,  only : occtol
 use m_energies,      only : ef
 use band,            only : bands

 use files, only : slabel
 use cellsubs,      only : reclat

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

 integer                    :: na
 real(dp)                   :: cell(3,3)
 real(dp)                   :: xa(3,na)
 real(dp)                   :: tcell(3,3)
 real(dp)                   :: fxa(3,na)

 ! reciprocal cell, irreducible k-points, number of (i)k-points, index of bands to be saved

 real(dp)                   :: rcell(3,3)
 real(dp),      allocatable :: kpoint_frac(:,:)
 real(dp),      allocatable :: kpoint_cart(:,:)
 integer                    :: nk, ink
 integer                    :: first_band, last_band

 ! spglib variables
 
 integer                   :: time_reversal
 real(dp)                  :: symprec
 type(SpglibDataset)       :: ds
 integer,      allocatable :: grid_address(:,:)
 integer,      allocatable :: map(:)

 ! fdf variables

 type(block_fdf)            :: bfdf
 type(parsed_line), pointer :: pline

 ! internal variables

 integer                    :: i, ii, iii, counter_kp_ibz
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
 ds = spg_get_dataset(tcell, fxa, isa, na, symprec)

 print*, 'siesta2bt:  Space Group:', ds%spacegroup_number

 print*, 'siesta2bt:  Number of symops:', ds%n_operations

 ! ===================================================================
 ! (3) Find irreducible k-points
 ! ===================================================================

 allocate(grid_address(3,nk), map(nk))

 ! get irreducible k-points (ik-points)
 ink = spg_get_ir_reciprocal_mesh(grid_address, &
                                  map,          &
                                  mesh,         &
                                  shift,     &
                                  time_reversal,  &
                                  tcell,        &
                                  fxa,          &
                                  isa,          &
                                  na,           &
                                  symprec)

 allocate(kpoint_frac(3,ink),kpoint_cart(3,ink))
 kpoint_frac=0.d0

 ! To account for the difference of indices in C and Fortran, add 1 to values of all
 ! elements of map

 map = map + 1

 ! Extracts k-points in the ibz. map contains the index of a point in
 ! grid_address; if iteration variable i and map(i) are equal, then this is one
 ! point in the IBZ.

 counter_kp_ibz = 0
 do i = 1, nk
    if (i == map(i)) then
      counter_kp_ibz = counter_kp_ibz + 1
      kpoint_frac(:,counter_kp_ibz) = dble(shift + 2*grid_address(:,map(i)))/dble(2*mesh)
    end if
 end do

 call reclat( cell, rcell, 1 )
 kpoint_cart = matmul(rcell,kpoint_frac)

 ! ===================================================================
 ! (4) Calculate bands
 ! ===================================================================

 allocate( ebk(last_band, spinor_dim, nk) )

 call bands( no_s, h_spin_dim, spinor_dim, no_u, no_l, maxnh, ink, &
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
   write(103,'(3F14.9,I5)') kpoint_frac(:,i), spinor_dim*(last_band - first_band)+1

   ! write energies
   do ii = 1, spinor_dim
     do iii = first_band, last_band
       write(103,'(F14.10)') ebk(iii,ii,i)
     enddo
   enddo
 enddo

 close(103)

end subroutine siesta2bt

