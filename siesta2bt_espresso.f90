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

 real(dp)                   :: tcell(3,3) ! transposed cell
 real(dp)                   :: fxa(3,na)  ! fractional atomic position

 ! irreducible k-points, number of k-points, num. of *_ikp

 real(dp),      allocatable :: kpoint(:,:)
 real(dp),      allocatable :: frac_kpoint(:,:)
 integer                    :: nk, maxnk
 integer                    :: ink
 real(dp),      allocatable :: kweight(:)
 real(dp)                   :: eff_kgrid_cutoff 


 ! spglib variables
 
 logical                   :: is_time_reversal_symmetry
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

 integer                    :: bt_rcell(3,3)
 integer                    :: mesh(3)
 real(dp)                   :: bt_shift(3)

 ! module and angule of lattice vectors

 real(dp)                   :: cellm(3), celang(3)
 real(dp),        parameter ::  pi = 3.1415926d0

 ! Band energies
 real(dp),    allocatable :: ebk(:,:,:)

 !transposing cell for c-fortran wrapper used in spglib
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
     bt_shift(i) = fdf_bvalues(pline,4)
   else
     bt_shift(i) = 0._dp
   end if
 enddo
 mesh(1) = bt_rcell(1,1)
 mesh(2) = bt_rcell(2,2)
 mesh(3) = bt_rcell(3,3)
 ! TODO: Conferir contribuição do bt_shift nessa conta
! nk = abs( bt_rcell(1,1) * bt_rcell(2,2) * bt_rcell(3,3) + &
!              bt_rcell(2,1) * bt_rcell(3,2) * bt_rcell(1,3) + &
!              bt_rcell(3,1) * bt_rcell(1,2) * bt_rcell(2,3) - &
!              bt_rcell(1,1) * bt_rcell(3,2) * bt_rcell(2,3) - &
!              bt_rcell(2,1) * bt_rcell(1,2) * bt_rcell(3,3) - &
!              bt_rcell(3,1) * bt_rcell(2,2) * bt_rcell(1,3) )   

! allocate( ebk(no_u, spinor_dim, nk) )

 ! read symprec from fdf file
 symprec = fdf_single('SPG.Symprec',1.e-2)
 
call cart2frac(na, xa(1,:), xa(2,:), xa(3,:), cell, fxa)


!enddo
! TODO: calculate rotation matrices, translation vectors
 ! and lattice type.
 ! ds usage:
 ! https://atztogo.github.io/spglib/api.html?highlight=kpoint#spg-get-dataset-and-spg-get-dataset-with-hall-number
! write(*,*) 'siesta2bt: About to run spg_get_dataset'
 write(*,*) 'siesta2bt: chamando spg_get_dataset'
 ds = spg_get_dataset(tcell, fxa, isa(1:na), na, symprec)
! write(*,*) 'siesta2bt: Dataset ok'

 is_time_reversal_symmetry = .false.
 
 nk = product(mesh)
 ink = nk
 write(*,*) 'siesta2bt: allocando vetor frac_kpoint(3,',ink,')'
 allocate(frac_kpoint(3,ink), kweight(ink))
 write(*,*) 'siesta2bt: chamando find_irreducible_kpoint'
 call find_irreducible_kpoints(mesh, bt_shift, is_time_reversal_symmetry,ds%rotations, ds%n_operations, frac_kpoint, kweight, ink)
 write(*,*) 'siesta2bt: ink', ink
! nk = ink
! allocate(kpoint(3,ink), frac_kpoint(3,ink))
 allocate( ebk(no_u, spinor_dim, ink) )

 ! get reciprical cell vectors
! call get_kpoints_scale('BandLinesScale',rcell,ierr)
!
! if (ierr /= 0) return

 ! calculate ik-prints from spglib's output as in the spglib's manual
! do i=1, ink
!   kpoint(:,i) = matmul(rcell,(grid_address(:,i) + is_shift/2.)/mesh)
!   frac_kpoint(:,i) = (grid_address(:,i) + is_shift/2.)/mesh
! enddo
 
 ! calculate eigenenergies
! call bands( no_s, h_spin_dim, spinor_dim, no_u, no_l, maxnh, nk, &
!             numh, listhptr, listh, H, S, ef, xijo, indxuo, &
!             .false., ink, kpoint, ebk, occtol, .false. )

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
   write(103,'(3F14.9,I5)') frac_kpoint(:,i), spinor_dim*no_u
   
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

     write(104,'(3F14.9)') frac_kpoint(:,i)
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
 write(104,*) bt_rcell(1,1),bt_rcell(2,2),bt_rcell(3,3)
 write(104,*) ''
 write(104,*) ''
 write(104,*) 'nkpoints', ink 
 write(104,*) ''
 write(104,*) 'kpoints:'
 write(104,'(4F14.9)') (frac_kpoint(:,i), kweight(i), i=1,ink)
 close(104)

 
end subroutine siesta2bt

subroutine find_irreducible_kpoints(mesh, shift, is_time_reversal,rot, nrot, kpoint, kweight, ink)
 use precision,    only : dp
 implicit none
 ! ***************** INPUT **********************************************
 ! integer  mesh(3)    : Division of reciprocal space
 ! real(dp) shift(3)   : Supercell reciprocal of k-grid unit cell
 ! integer  nrot       : Numer of rotation matrixes
 ! integer  rot        : Rotations matrixes
 ! ***************** IN/OUTPUT *********************************************
 ! integer  ink        : In: Number of k-points
 !                     : Out: Number of irreducible k-points 
 ! ***************** OUTPUT *********************************************
 ! real(dp) kpoints(3,ink) : K-points in terms of reciprocal vectors
 ! real(dp) kweight(ink)   : weights
 !
 ! (1) calculos como no quantum espresso
 ! (2) nossos calculos
 !
 ! Input
 integer,   intent(in) :: mesh(3)
 real(dp),  intent(in) :: shift(3)
 integer,   intent(in) :: nrot      
 integer,   intent(in) :: rot(3,3,nrot) 
 logical,   intent(in) :: is_time_reversal
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
! ds2 = minval(mesh/2.d0)
 ds2 = 1.0E-8
 !
 nk = product(mesh)
 !
!inv_xyz = find_rot_mat(rot,(/-1,0,0,0,-1,0,0,0,-1/)) ! (2)
 !
 !
 do i=1, mesh(1)
   do j=1, mesh(2)
     do k=1, mesh(3)
       n = (k-1) + (j-1)*mesh(3) + (i-1)*mesh(2)*mesh(3) + 1
       kpoint(:,n) = dble((/i,j,k/)-1)/mesh + dble(shift)/mesh            ! (1)
!      print*, 'siesta2bt: i,j,k : ', i,j,k
!      print*, 'siesta2bt: dble((/i,j,k/)-1)/mesh: ', dble((/i,j,k/)-1)/mesh
!      print*, 'siesta2bt: dble(shift)/mesh: ', dble(shift)/mesh
!      kpoint(:,n) = dble( 2.d0*( (/i,j,k/) + shift ) - mesh - 1.d0)/mesh ! (2)
     enddo 
   enddo 
 enddo 
 !
 do n=1,nk
   equiv(n) = n
 enddo
 !
 kweight = 1.d0
 !
 print*,'entrando no laço de nk',nk
 do l=1, nk
   print*,'dentro do laço entrando no equiv'
   if (equiv(l) == l) then
     print*,'entrei no equiv'
     do m = 1, nrot
!      kpoint_rot = matmul(transpose(rot(:,:,m)),kpoint(:,l))
       kpoint_rot = matmul(rot(:,:,m),kpoint(:,l))
!      print('(A,I3,(3I3),/)'),'mat:  ',m,  rot(:,:,m)
       print('(A,I3,A,3I3)'),  'mat:  ',m,':', rot(1,:,m)
       print('(A,3I3)'),       '          ',    rot(2,:,m)
       print('(A,3I3)'),       '          ',    rot(3,:,m)
       print('(A,I3,3F7.3)'),'kpoint    : ',l,kpoint(:,l)
!      print('(A,3F7.3)'),   'kpoint rot: ',kpoint_rot(:)
       kpoint_rot = kpoint_rot - nint(kpoint_rot)
       print('(A,I3,3F7.3)'),'kpoint rot - nint(kpoint_rot): ',l,kpoint_rot(:)
       !
       rr = kpoint_rot*mesh - shift
       print('(A,3F10.4)'), 'rr: ',rr
       print*, 'maxval(dabs(rr - nint(rr))/2.d0): ', maxval(dabs(rr - nint(rr))/2.d0)
       in_the_list = maxval(dabs(rr - nint(rr))/2.d0) <= ds2
       print('(A,L)'), 'in_the_list: ', in_the_list
       if (in_the_list) then
         kindex_rot = mod ( nint ( kpoint_rot*mesh - shift + 2*mesh ) , mesh ) + 1
         n = (kindex_rot(3)-1) + (kindex_rot(2)-1)*mesh(3) + (kindex_rot(1)-1)*mesh(2)*mesh(3) + 1
!        print('(A,I3)'), 'kindex_rot: ', kindex_rot
!        print('(A,I3)'), 'n: ', n
         !
         print*,n,'>',l,' and ',equiv(n),' == ',n
         print*, ' '
         if (n > l .and. equiv(n) == n) then
           equiv(n) = l
           kweight(l) = kweight(l) + 1.0d0
          print*,'l: ',l
          print*,'kweight(l): ',kweight(l)
         else
           IF (equiv(n)/=l .or. n<l ) write(*,*) 'siesta2bt: Erro kpoint_grid!'
         endif
         !
         if (is_time_reversal) then
           rr = -kpoint_rot*mesh - shift
           in_the_list = maxval(dabs(rr - nint(rr))/2.d0) <= ds2
           !
           if (in_the_list) then
             kindex_rot = mod ( nint ( - kpoint_rot*mesh - shift + 2*mesh ) , mesh ) + 1
             n = (kindex_rot(3)-1) + (kindex_rot(2)-1)*mesh(3) + (kindex_rot(1)-1)*mesh(2)*mesh(3) + 1
             !
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
     print*,'n:       ', n
     print*,'ink:     ', ink
     print*,'wk(ink): ', kweight(ink)
     print*,'fact:    ', sum(kweight(1:ink))

     print*, 'before'
     print('(A,3F7.3)'), 'kpoint(ink): ',kpoint(:,n)
     kpoint(:,ink) = kpoint(:,n) - nint(kpoint(:,n))
     print*, 'after'
     print('(A,3F7.3)'), 'kpoint(ink): ',kpoint(:,ink)

   endif
 enddo
 print*, 'ink:  ',ink
 print*, 'fact: ',sum(kweight(1:ink))
 print*, 'before'
 print*,'wk(1:ink): ', kweight(1:ink)
 kweight(1:ink) = 2*kweight(1:ink)/sum(kweight(1:ink))
 print*, 'after'
 print*,'wk(1:ink): ', kweight(1:ink)

end subroutine find_irreducible_kpoints
