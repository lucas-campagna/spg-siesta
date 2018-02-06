!
! Program to generate output for BoltzTrap
! 


subroutine siesta2bt(cell, xa, na)

 use spglib_f08

 use Kpoint_grid,  only : kscell, nkpnt, kpoint


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

! integer                    :: cna        ! number of atoms of conventional cell
! real(dp)                   :: ccell(3,3) ! conventional cell
 real(dp)                   :: tcell(3,3) ! fractional cell
 real(dp)                   :: fxa(3,na)  ! fractional atomic position
! real(dp), allocatable      :: cfxa(:,:)  ! fractional atomic position of conventional cell
! integer,  allocatable      :: cisa(:)    ! atomic type of conventional cell

 ! irreducible k-points, number of k-points, num. of *_ikp

 real(dp),      allocatable :: cart_ikp(:,:), frac_ikp(:,:)
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

 integer                    :: i, ii, iii
 logical                    :: bt_calc
 integer                    :: ierr
 integer                    :: cont
 real(dp)                   :: deltae
 real(dp)                   :: maxbv, minbc

 ! cell, mesh and shift of reciprocal space

 real(dp)                   :: rcell(3,3)
 integer                    :: bt_rcell(3,3)
 real(dp)                   :: bt_srcell(3)

 ! module and angule of lattice vectors

 real(dp)                   :: cellm(3), celang(3)
 real(dp),        parameter ::  pi = 3.1415926d0

 ! Band energies
 real(dp),    allocatable :: ebk(:,:,:)

 ! transposing cell for c-fortran wrapper used in spglib
 tcell = transpose(cell)

 ! read bloc BT.kfrid_Monkhosrt_Pack
 bt_calc = fdf_block('BT.kgrid_Monkhorst_Pack',bfdf)

 if ( .not. bt_calc ) then
  return
 endif
 
 ! read block BT.kfrid_Monkhosrt_Pack
 do i=1, 3
   if (.not. fdf_bline(bfdf,pline)) then
     call die('siesta2bt: ERROR in BT.kgrid_Monkhorst_Pack block')
   end if
     bt_rcell(1,i) = fdf_bintegers(pline,1)
     bt_rcell(2,i) = fdf_bintegers(pline,2)
     bt_rcell(3,i) = fdf_bintegers(pline,3)
   if ( fdf_bnvalues(pline) > 3 ) then
     bt_srcell(i) = fdf_bvalues(pline,4)
   else
     bt_srcell(i) = 0._dp
   end if
 enddo
 ! TODO: Conferir contribuição do bt_srcell nessa conta
 maxnk = abs( bt_rcell(1,1) * bt_rcell(2,2) * bt_rcell(3,3) + &
              bt_rcell(2,1) * bt_rcell(3,2) * bt_rcell(1,3) + &
              bt_rcell(3,1) * bt_rcell(1,2) * bt_rcell(2,3) - &
              bt_rcell(1,1) * bt_rcell(3,2) * bt_rcell(2,3) - &
              bt_rcell(2,1) * bt_rcell(1,2) * bt_rcell(3,3) - &
              bt_rcell(3,1) * bt_rcell(2,2) * bt_rcell(1,3) )   

 allocate( ebk(no_u, spinor_dim, maxnk) )

 ! read symprec from fdf file
 symprec = fdf_single('SPG.Symprec',1.e-2)
 
! if(coord_format/='Fractional') then
   call cart2frac(na, xa(1,:), xa(2,:), xa(3,:), cell, fxa)
! else
!   fxa = xa
! endif

 ! TODO: simetrize structure
 ! spg_find_primitive(cell, fxa, isa, na, symprec)
 ! ccell = tcell

 ! normalizing to spg_refine_cell (find conventional cell)
! do i = 1, 3
!  cellm(i) = 0.d0
!  do ii = 1, 3
!    cellm(i) = cellm(i) + ccell(ii,i)*ccell(ii,i)
!  enddo
!  cellm(i) = sqrt(cellm(i))
!enddo
!ccell = ccell/maxval(cellm) * 0.9
 
! allocate(cfxa(3,na*4), cisa(na*4))

! cfxa = 0.0
! cisa = 0
! cfxa(:,1:na) = fxa
! cisa(1:na) = isa

! write(*,*) 'siesta2bt: About to run spg_refine_cell'
! cna = spg_refine_cell(ccell, cfxa, cisa, na, symprec)
! write(*,*) 'siesta2bt: cell refined'

! ccell = ccell*maxval(cellm)/0.9

! TODO: calculate rotation matrices, translation vectors
 ! and lattice type.
 ! ds usage:
 ! https://atztogo.github.io/spglib/api.html?highlight=kpoint#spg-get-dataset-and-spg-get-dataset-with-hall-number
! write(*,*) 'siesta2bt: About to run spg_get_dataset'
 ds = spg_get_dataset(tcell, fxa, isa(1:na), na, symprec)
! write(*,*) 'siesta2bt: Dataset ok'

 ! get irredutible k-points (ik-points)
 mesh(1) = bt_rcell(1,1) 
 mesh(2) = bt_rcell(2,2) 
 mesh(3) = bt_rcell(3,3) 

 nk = maxnk
 allocate(grid_address(3,nk),map(nk))

 is_shift = bt_srcell*2
 is_time_reversal = 0

 write(*,*) 'siesta2bt: input of spg_get_ir_reciprocal_mesh'
 write(*,'(A18,3I3)') 'mesh: ', mesh 
 write(*,'(A18,3I3)') 'is_shift: ', is_shift 
 write(*,'(A18,3I3)') 'is_time_reversal: ', is_time_reversal
 write(*,'(A18/,(3F16.9))') 'lattice: ', cell 
 write(*,'(A18/,(3F16.9))') 'position: ', fxa 
 write(*,'(A18,(100I2))') 'type: ', isa(1:na)
 write(*,'(A18,I2)') 'na: ', na
 write(*,'(A18,E16.9)') 'symprec: ', symprec
 write(*,*)

 ink = spg_get_ir_reciprocal_mesh(grid_address, &
                                  map,          &
                                  mesh,         &
                                  is_shift,     &
                                  is_time_reversal,  &
                                  cell,         &
                                  fxa,          &
                                  isa(1:na),    &
                                  na,           &
                                  symprec)
 write(*,*) 'siesta2bt: output of spg_get_ir_reciprocal_mesh'
 write(*,'(A18,3I3)') 'mesh: ', mesh 
 write(*,'(A18,3I3)') 'is_shift: ', is_shift 
 write(*,'(A18,3I3)') 'is_time_reversal: ', is_time_reversal
 write(*,'(A18/,(3F16.9))') 'lattice: ', cell 
 write(*,'(A18/,(3F16.9))') 'position: ', fxa 
 write(*,'(A18,(100I2))') 'type: ', isa(1:na)
 write(*,'(A18,I2)') 'na: ', na
 write(*,'(A18,E16.9)') 'symprec: ', symprec
 write(*,'(A18,I4)') 'ink: ', ink
 write(*,*)

 maxnk = ink
 allocate(cart_ikp(3,ink), frac_ikp(3,ink))

 ! get reciprical cell vectors
 call get_kpoints_scale('BandLinesScale',rcell,ierr)

 if (ierr /= 0) return

 ! calculate ik-prints from spglib's output as in the spglib's manual
 do i=1, ink
   cart_ikp(:,i) = matmul(rcell,(grid_address(:,i) + is_shift/2.)/mesh)
   frac_ikp(:,i) = (grid_address(:,i) + is_shift/2.)/mesh
 enddo
 
 ! calculate eigenenergies
 call bands( no_s, h_spin_dim, spinor_dim, no_u, no_l, maxnh, maxnk, &
             numh, listhptr, listh, H, S, ef, xijo, indxuo, &
             .false., ink, cart_ikp, ebk, occtol, .false. )

 ! calculate modules and angules between the cell-vectors
 ! Modules
 do i = 1, 3
   cellm(i) = 0.d0
   do ii = 1, 3
     cellm(i) = cellm(i) + cell(ii,i)*cell(ii,i)
   enddo
   cellm(i) = sqrt(cellm(i))
 enddo

 ! Angles
 celang(1) = 0.d0
 do i = 1, 3
   celang(1) = celang(1) + cell(i,1)*cell(i,2)
 enddo
 celang(1) = acos(celang(1)/(cellm(1)*cellm(2)))*180.d0/pi
 celang(2) = 0.d0
 do i = 1, 3
   celang(2) = celang(2) + cell(i,1)*cell(i,3)
 enddo
 celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180.d0/pi
 celang(3) = 0.d0
 do i = 1, 3
   celang(3) = celang(3) + cell(i,2)*cell(i,3)
 enddo
 celang(3) = acos(celang(3)/(cellm(2)*cellm(3)))*180.d0/pi


 ! save BoltzTrap format:
 ! ===============================================
 ! (1) intrans file 
 open(101,file=trim(slabel)//'.intrans')
 
 write(101,'(A)') 'GENE                      # Format of DOS' 
 write(101,'(A)') '0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap'
 write(101,'(F7.5,A,F6.0,A)') Ef, ' 0.0005 0.4 ',qtot,&
            ' # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons'
 write(101,'(A)') 'CALC                      # CALC (calculate expansion coeff), NOCALC read from file'
 write(101,'(A)') '3                         # lpfac, number of latt-points per k-point'
 write(101,'(A)') 'BOLTZ                     # run mode (only BOLTZ is supported)'
 write(101,'(A)') '0.15                      # (efcut) energy range of chemical potential'
 write(101,'(A)') '800. 50.                  # Tmax, temperature grid'
 write(101,'(A)') '-1                        # energyrange of bands given individual DOS output'
 write(101,'(A)') 'HISTO                     # sig_xxx and dos_xxx (xxx is band number)'

 close(101) 

 ! ===============================================
 ! (2) struct file
 open(102,file=trim(slabel)//'.struct')

 ! title 
 write(102,'(A)') trim(slabel)
 
 ! lattice vectors (Bohr)
 !write(102,'(6F20.14)') cellm, celang
 !TODO!!!
 write(102,'(3F16.9)') cell(:,1)
 write(102,'(3F16.9)') cell(:,2)
 write(102,'(3F16.9)') cell(:,3)

 ! number of operations
 write(102,'(I3)') ds%n_operations

 do i=1, ds%n_operations
   write(102,'(9I3)') ds%rotations(:,:,i)
!   write(102,'(3I4)') ds%rotations(:,1,i)
!   write(102,'(3I4)') ds%rotations(:,2,i)
!   write(102,'(3I4)') ds%rotations(:,3,i)
!   write(102,'(I4)') i
 enddo

 close(102)
 ! ===============================================
 ! (3) energy file
 open(103,file=trim(slabel)//'.energy')
 open(104,file=trim(slabel)//'.denergy')
 
 ! title 
 write(103,*) trim(slabel)

 ! number of irredutible kpoints
 write(103,'(I8)') ink
 
 ! write ik-points
 !TODO: implementar leitura do no_u pela flag BT.NumberOfEigenStates
 !TODO: implementar um range de energia ou range de bands pra serem escritos neste arquivo
 no_u=8
 minbc = 200
 maxbv = -200
 do i=1, ink
   write(103,'(3F14.9,I5)') frac_ikp(:,i), spinor_dim*no_u
   
   ! write energies
   do ii=1,spinor_dim
     do iii=1, no_u
       write(103,'(F14.10)') ebk(iii,ii,i)
     enddo
     deltae=ebk(1,ii,i) - Ef
     cont=1
     do while(deltae <= 0.000)
       deltae = ebk(cont,ii,i) - Ef
       cont = cont + 1
     enddo

     write(104,'(3F14.9)') frac_ikp(:,i)
     write(104,'(A,F14.9,A,F14.9,A,F14.9,A,I5)')'maxbv: ', ebk(cont-1,ii,i),'   minbc: ', ebk(cont,ii,i), &
     '   egap: ', ebk(cont,ii,i) - ebk(cont-1,ii,i), '    cont:  ', cont
     write(104,*) 
     minbc = min(minbc,ebk(cont,ii,i))
     maxbv = max(maxbv,ebk(cont-1,ii,i))
   enddo
 enddo
 write(104,*) 
 write(104,*) 'Global'
 write(104,'(A,F14.9,A,F14.9,A,F14.9)')'maxbv: ', maxbv,'   minbc: ', minbc, '   egap: ', minbc-maxbv

 close(104)
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
 write(104,*) 'mesh:'
 write(104,*) mesh
 write(104,*) ''
 write(104,*) 'rcell'
 write(104,'(3(3F14.9,/))') kscell
 write(104,*) ''
 write(104,*) 'nkpoints', nkpnt 
 write(104,*) ''
 write(104,*) 'kpoints:'
 write(104,'(3F14.9)') kpoint
 close(104)

 
end subroutine siesta2bt

