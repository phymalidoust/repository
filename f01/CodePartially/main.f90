PROGRAM BDG
  !use mpi
  use constants
  use field_constants
  use user_input
  use cfield_class
  use rfield_class
  use utilities
  use spin_flip_procedures
  use eigen_vector_mod
  use iso_fortran_env
  use hdf5_utilities
  
  IMPLICIT NONE

  include "mpif.h"

 !+++++
  integer(kind=4) :: ierr             ! MPI communication variables
  integer, dimension(:), allocatable   :: numstates_loc

  type(Distributed_Matrix)    :: A        ! the total matrix to be evaluated at each x momentum

 !          ***************** 2D BdG  f1f2S
  real    :: rerror
  logical :: converged
  real    :: lambda,factors1,factors2,factors3,factorsj,factorsn
  real    :: ratio1,nfe,sfe,nfree,nfree_loc,sfree,sfree_loc,freediff,avg
  integer :: i,j,it
  integer :: kx_idx ! index for a given x momentum
  integer :: F_idx  ! index for a F section within a field (F1, F2, F3)
  integer :: k,kp,p,q,qp
  integer :: k_idx, kp_idx
  integer :: nstates
  logical :: exist_already
  real, dimension(:,:), allocatable :: divx, divy, divz
  real    :: divxsum,divysum,divzsum
  real    :: const_factor2, const_factor3
  integer, dimension(2) :: c1_bounds, c2_bounds
  character(99) :: tmp_fname,tmp_fname_2 
  complex :: tmp_cmplx
  complex :: val
  complex, dimension(:), allocatable :: column14
  complex, dimension(:), allocatable :: column23
  integer :: clustersize
  real :: nphase

! rfield holds a real array that has dimensions (0:npy-1,0:npz-1)
  type(rfield)  :: diff
  type(rfield)  :: tx, ty, tz
  type(rfield)  :: hxvec, hyvec, hzvec
  type(rfield)  :: tmp_rfield
  type(rfield)  :: phase

  type(cfield)  :: old_delta, delta
  type(cfield)  :: mx, my, mz, rmx, rmy, rmz 
  type(cfield)  :: cooper0, cooper1, cooper2, cooper3
  type(cfield)  :: ret_cooper0, ret_cooper1, ret_cooper2, ret_cooper3
  type(cfield)  :: cooper0_r, cooper1_r, cooper2_r
  type(cfield)  :: ret_delta
  type(cfield)  :: nup, ndown, retnup, retndown, density
  type(cfield)  :: jyup, jydown, retjyup, retjydown
  type(cfield)  :: jytot, jztot
  type(cfield)  :: jzup, jzdown, retjzup, retjzdown
  type(cfield)  :: jyxspin, jyyspin, jyzspin
  type(cfield)  :: retjyx, retjyy, retjyz, retjzx, retjzy, retjzz
  type(cfield)  :: jzxspin, jzyspin, jzzspin

! HDF5 variables
  integer(hid_t)    :: h5_file
  real, dimension(:), allocatable :: tmp_arr
  real, dimension(:,:), allocatable :: tmp_arr_2d

! Density of states variables
  integer           :: n_DOS
  integer           :: kk,nn
  type(rfield), dimension(:,:), allocatable :: DOS
  type(rfield) :: DOS_tot
  real, dimension(:), allocatable :: DOS_energies

  ! Initializations
  ! MPI
  call MPI_Init(ierr)

  ! User input
  call ui_init_from('input.txt', i)
  if(i /= 0) then
    ! no input file, use defaults (small problem)
    call ui_set_defaults
    call ui_set_derived_values()
  endif

  h5_file = -1

  if(comm_rank == 0) then
    open(unit=76, file='Job_parameters.txt')
    call ui_write_constants(76)
    close(76)
  endif

  ! spin flip procedures
  call initialize_spin_flip_procedures()

  ! Timer
  start_time = MPI_Wtime()

  n_DOS = nint((abs(real(DOS_range(2))-real(DOS_range(1))))/DOS_inter)+1

  ! allocations
  allocate(numstates_loc    (0:nkx-1))
  allocate(column14         (qn_local_rows))
  allocate(column23         (qn_local_rows))
  allocate(DOS              (2,n_DOS))
  allocate(DOS_energies     (n_DOS))

  DOS_tot = rfield(0.0)

  DOS_energies(1) = DOS_range(1)
  do k=2,n_DOS
     DOS_energies(k) = DOS_energies(k-1) + DOS_inter
  end do

  do k=1,n_DOS
    DOS(1, k) = rfield(0.0)
    DOS(2, k) = rfield(0.0)
  end do

  diff       = rfield(0.0)
  tx         = rfield(0.0)
  ty         = rfield(0.0)
  tz         = rfield(0.0)
  hxvec      = rfield(0.0)
  hyvec      = rfield(0.0)
  hzvec      = rfield(0.0)
  tmp_rfield = rfield(0.0)
  phase      = rfield(0.0)

  old_delta   = cfield(cZero)
  delta       = cfield(cZero)
  mx          = cfield(cZero)
  my          = cfield(cZero)
  mz          = cfield(cZero)
  rmx         = cfield(cZero)
  rmy         = cfield(cZero)
  rmz         = cfield(cZero) 
  cooper0     = cfield(cZero)
  cooper1     = cfield(cZero)
  cooper2     = cfield(cZero)
  cooper3     = cfield(cZero)
  cooper0_r   = cfield(cZero)
  cooper1_r   = cfield(cZero)
  cooper2_r   = cfield(cZero)
  ret_cooper0 = cfield(cZero)
  ret_cooper1 = cfield(cZero)
  ret_cooper2 = cfield(cZero)
  ret_cooper3 = cfield(cZero)
  ret_delta   = cfield(cZero)
  nup         = cfield(cZero)
  ndown       = cfield(cZero)
  retnup      = cfield(cZero)
  retndown    = cfield(cZero)
  jyup        = cfield(cZero)
  jydown      = cfield(cZero)
  retjyup     = cfield(cZero)
  retjydown   = cfield(cZero)
  jytot       = cfield(cZero)
  jztot       = cfield(cZero)
  jzup        = cfield(cZero)
  jzdown      = cfield(cZero)
  retjzup     = cfield(cZero)
  retjzdown   = cfield(cZero)
  jyxspin     = cfield(cZero)
  jyyspin     = cfield(cZero)
  jyzspin     = cfield(cZero)
  retjyx      = cfield(cZero)
  retjyy      = cfield(cZero)
  retjyz      = cfield(cZero)
  retjzx      = cfield(cZero)
  retjzy      = cfield(cZero)
  retjzz      = cfield(cZero)
  jzxspin     = cfield(cZero)
  jzyspin     = cfield(cZero)
  jzzspin     = cfield(cZero)

  call status_message('Setting values of H')


  ! Initialize Eigen Vector module (allocates memory)
  clustersize = nint(blocking * sqrt(real(blocking)))
  if (nkx == 1 .and. kx_factor > 2.0) clustersize = blocking * blocking
  call init_eigen_vector_module(tot_N, n_proc_rows, n_proc_cols, nkx, blocking, clustersize)
  call file_status_message("Eigen Vector module initialized")

  if(mod(nkx, n_grids()) /= 0 .and. iam_root()) then
    call time_stamp()
    write(*,*) "nkx % ngrids != 0"
    write(*,*) "nkx =", nkx, ", n_grids =", n_grids()
    error stop 102
  endif

  ! Initialize HDF5 and create h5 file if it doesn't already exist
  if(comm_rank == 0) then
     call h5_init()
     inquire(file='spin_flip.h5',exist=exist_already)     
     if(.not. exist_already) h5_file = create_hdf5('spin_flip.h5', .true.)
  end if

  ! Don't execute if the user has set imp_conc to zero
  if (tot_imp > 0) then
     if(comm_rank == 0) then
        ! Initialize random number function
        if (rand_opt == 2) then
           call date_and_time(VALUES=date_info)
           rand_seed = sum(date_info)
        endif
        call srand(rand_seed)

        ! Open a file to write the impurity locations
        open(unit=67, file='scattering.txt')
        write(67,'("random seed for this run: ",i8)') rand_seed
        write(67,'(" ")')

        allocate(tmp_arr(tot_imp))


        imp_idx = 0
        do k=S1,S2
           do nn = 1,num_imp(k)
              imp_idx = imp_idx + 1 
              imp_loc_y(imp_idx) = nint(L_idx(k)+rand()*(R_idx(k)-L_idx(k)))
              imp_loc_z(imp_idx) = nint(rand()*(npz-1))
              ! Check if the new impurity is a duplicate
              do j = 1,nn
                 if ((imp_loc_y(nn) == imp_loc_y(j)) .and. (imp_loc_z(nn) == imp_loc_z(j))) then
                    imp_loc_y(imp_idx) = imp_loc_y(imp_idx) + 1
                    if (imp_loc_y(imp_idx) > R_idx(k)) &
                         imp_loc_y(imp_idx) = imp_loc_y(imp_idx) - 2
                    exit
                 endif
              enddo
              write(67,'(2(i6,2x))') imp_loc_y(imp_idx),imp_loc_z(imp_idx)
           end do
           write(67,'(" ")')
        end do

        h5_file = open_hdf5('spin_flip.h5', .true.)

        do i = 1,tot_imp
           tmp_arr(i) = real(imp_loc_y(i))*dY
        enddo
        call write_real_array_to_hdf5(h5_file, 'imp_loc_y', tmp_arr)

        do i = 1,tot_imp
           tmp_arr(i) = real(imp_loc_z(i))*dZ
        enddo
        call write_real_array_to_hdf5(h5_file, 'imp_loc_z', tmp_arr)

        call close_hdf5(h5_file)

        deallocate(tmp_arr)

        close(67)
     end if !comm_rank == 0

     ! Let all other cores know about where the scattering centers are placed
     call MPI_BCAST(imp_loc_y,num_imp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(imp_loc_z,num_imp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  endif

  ! set the row and column indices used by this rank
  local_rows = row_indices()
  local_cols = col_indices()

  ratio1  = omega/bulkdel
  lambda  = 1.0 / (LOG(ratio1 + SQRT(1.0+ratio1**2)))
  deltakx = max_kx/real(nkx)                            !zero is included so must subtract 1 to get delta

  factors1= 8.0*pi*(lambda*deltakx)/(kfw*kfl*bulkdel)   !prefactors for pair potential sum
  factors2=-2.0*deltakx*pi/(bulkdel*kfw*kfl)            !free energy
  factors3= 8.0*deltakx*pi/(kfw*kfl)                    !prefactors for mag. mom. sum	
  factorsj= 8.0*pi*deltakx/(kfw*kfl)                    !current prefactors
  factorsn= 4.0*deltakx*pi/(bulkdel*kfw*kfl)            !prefactors for particle density

  T = 2.0*temp_fact*temp*bulkdel

  do F_idx = F1,F3
    call hxvec%set_section(F_idx, H(1, F_idx))
    call hyvec%set_section(F_idx, H(2, F_idx))
    call hzvec%set_section(F_idx, H(3, F_idx))
  end do

  if(comm_rank == 0) then
    call hxvec%write_to_file('hxvec.txt')
    call hyvec%write_to_file('hyvec.txt')
    call hzvec%write_to_file('hzvec.txt')

    h5_file = open_hdf5('spin_flip.h5', .true.)

    tmp_arr = [(i*dY, i=0,npy-1)]
    call write_real_array_to_hdf5(h5_file, 'y_coord', tmp_arr)

    tmp_arr = [(i*dZ, i=0,npz-1)]
    call write_real_array_to_hdf5(h5_file, 'z_coord', tmp_arr)

    call close_hdf5(h5_file)
  endif

  ! set [c]psi{l|w}
  const_factorp = pi/real(npy-1)
  const_factorq = pi/real(npz-1)
  const_factor2 = pi/kfl
  const_factor3 = pi/kfw

  call status_message('Setting [c]psi{l|w}')
  call file_status_message('Setting [c]psi{l|w}')
  FORALL(p=1:n_p,i=0:npy-1)  psil(p,i)=                  SIN(i*p*const_factorp)
  FORALL(p=1:n_p,i=0:npy-1) cpsil(p,i)=(p*const_factor2)*COS(i*p*const_factorp)
  FORALL(q=1:n_q,i=0:npz-1)  psiw(q,i)=                  SIN(i*q*const_factorq)
  FORALL(q=1:n_q,i=0:npz-1) cpsiw(q,i)=(q*const_factor3)*COS(i*q*const_factorq)

  do q=1,n_q / n_proc_rows
    qp = (q-1)*n_proc_rows + grid_row() + 1
    psiw_local(q,:)  =  psiw(qp,:)
    cpsiw_local(q,:) = cpsiw(qp,:)
  end do

  ! tranpose them for hfunct and cfunct
  call file_status_message('Transposing [c]psi{l|w}')
  call status_message('Transposing [c]psi{l|w}')
  psil_trans  = transpose(psil)
  cpsil_trans = transpose(cpsil)
  psiw_trans  = transpose(psiw)
  cpsiw_trans = transpose(cpsiw)
  call file_status_message('Finished transposing [c]psi{l|w}')

  FORALL(p=1:n_p) hlookupl(p)=(p*const_factor2)**2
  FORALL(q=1:n_q) hlookupw(q)=(q*const_factor3)**2

  converged=.false.
  rerror=minerror+.01 !ensures that it starts higher

  ! set the (1,2) and (3,4) portions of A
  call status_message('Creating Distributed Matrix')
  call file_status_message('Creating Distributed Matrix')
  A = Distributed_Matrix(cZero)
  call status_message('Finished creating Distributed Matrix')
  call file_status_message('Finished creating Distributed Matrix')

  call file_status_message('About to set (1,2) and (3,4) values of A')
  c1_bounds = [                1,   qn_local_rows]  ! section 1
  c2_bounds = [2*qn_local_rows+1, 3*qn_local_rows]  ! section 3

  !$omp parallel private(kp, F_idx, val)
  do kp_idx=1,qn_local_cols
    kp = local_cols(kp_idx)

    !$omp sections
    !$omp section
    column14 = 0.0
    !$omp section
    column23 = 0.0
    !$omp end sections

    !$omp do
    do k_idx = 1,qn_local_rows
      k = local_rows(k_idx)
      val = 0
      do F_idx = F1,F3
        !                      X direction                  Y direction
        val = val - hfunct(1, F_idx, k, kp) + im*hfunct(2, F_idx, k, kp)
      end do

      column14(k_idx) = val               ! (1,2)
      column23(k_idx) = conjg(val)        ! (3,4)
    end do
    !$omp end do

    !$omp sections
    !$omp section
    call A%set_local_col_section(kp_idx + qn_local_cols, c1_bounds(1), c1_bounds(2), column14)
    !$omp section
    call A%set_local_col_section(kp_idx + 3*qn_local_cols, c2_bounds(1), c2_bounds(2), column23)
    !$omp end sections
  end do
  !$omp end parallel

  call file_status_message("Finished setting (1,2) and (3,4) values of A")

  nfree=0
  nfree_loc=0


  call status_message("About to call normalstate")
  call file_status_message("About to call normalstate")

  h5_file = open_hdf5('spin_flip.h5', iam_root())

  call write_eigen_values('W_initial.txt', 'W_initial_min.txt')
  call write_eigen_values_to_hdf5(h5_file, 'initial_eigen_values')

  call status_message("Finished calling normalstate")
  call file_status_message("Finished calling normalstate")

  call close_hdf5(h5_file)

 !+++++
  numstates_loc = numstates
  call MPI_Allreduce(   numstates_loc, numstates, nkx, MPI_INTEGER, &
                        MPI_SUM, rank_comm(), ierr)
 !-----	        

  ! sfree will be reduced over MPI_COMM_WORLD
  call MPI_Allreduce(   nfree_loc,nfree,1,MPI_REAL, &
                        MPI_SUM, rank_comm(), ierr)
  nfree=nfree*factors2

  if(iam_root()) then
    open(2,file='normal.txt')
    do j=0,nkx-1
      write(2, 22) numstates(j)
    end do
    close(2)

    h5_file = open_hdf5('spin_flip.h5', .true.)
    call write_int_array_to_hdf5(h5_file, 'normal', numstates)
    call close_hdf5(h5_file)

    open(123,file='nfree.txt')
    write(123,23) nfree
    close(123)
  end if

  22 format(1x,1i5)    
  23 format('n',4x,1f11.7)

  inquire(file=g_txt,exist=exist_already)
  if (restart_flag .and. exist_already) then
     !inquire(file=g_txt,exist=exist_already)
     !if (.not. exist_already) error stop "Could not find g.txt for restart"
     call file_status_message("Reading previous run data from g.txt")
     call old_delta%read_from_file(g_txt)
  else
     call file_status_message("Running without reading g.txt")
     call old_delta%set(cZero)
     call old_delta%set_section(S1, exp(im*vphi1))
     call old_delta%set_section(S2, exp(im*vphi2))
     if(comm_rank == 0) then
        h5_file = open_hdf5('spin_flip.h5', .true.)
        call old_delta%write_to_file('delta-00.txt')
        call old_delta%write_to_hdf5(h5_file, 'delta-00')
        call close_hdf5(h5_file)
     endif
  endif

  ! reset row bounds from row section 3 to row section 2, to set (2,3)
  c2_bounds = [qn_local_rows+1, 2*qn_local_rows]

  DO it=1,itmax !*******start s-c loop*******************
    if (rerror<minerror .or. it==itmax) then
      converged=.true.
    endif

    call file_status_message("Starting Iteration ", it)
    call status_message('Starting Iteration ', it)
    
    call delta%set(cZero)
    
    ! Cacluate the new value of delta (for iterations), This will go in blocks
    ! Override values of delta if kfdx is present
    if (npsx > 0) then
       do i = 0, npsx-1
          do j = 0,npz-1
             call old_delta%set(j,i,exp(im*vphi1))
          enddo
       enddo
       do i = R_idx(S2)-npsx,npy-1
          do j = 0,npz-1
             call old_delta%set(j,i,exp(im*vphi2))
          enddo
       enddo
    endif

    ! (1,4) and (2,3) of the Hamiltonian
    kp=0

    !$omp parallel private(kp_idx, kp, k_idx, k, val)
    do kp_idx = 1,qn_local_cols
      kp = local_cols(kp_idx)

      !$omp single
      column14 = cZero
      column23 = cZero
      !$omp end single
      !$omp barrier
    
      !$omp do
      do k_idx = 1,qn_local_rows
        k = local_rows(k_idx)

        val = cfunct(k,kp,old_delta)

        column14(k_idx) = val    ! row block 1, col block 4, [del](k,kp)
        column23(k_idx) = val    ! row block 2, col block 3, [del](k,kp)
      end do
      !$omp end do

      ! A is set to 0 at start of interation, so the lines below are equivilent to setting them

      !$omp sections
      !$omp section
      call A%set_local_col_section(kp_idx+3*qn_local_cols, c1_bounds(1), c1_bounds(2), column14)
      !$omp section
      call A%set_local_col_section(kp_idx+2*qn_local_cols, c2_bounds(1), c2_bounds(2), column23)
      !$omp end sections
    end do
    !$omp end parallel

    if (.not. converged) then

       call   jyup%set(cZero)
       call jydown%set(cZero)
       call cooper3%set(cZero)
       sfree = 0
       sfree_loc = 0
       avg = 0

      call file_status_message("  Calling deltafind on all kx")
      call status_message("  Calling deltafind on all kx")
      ! call deltafind on all of the nkx assigned to this rank
      do kx_idx = start_kx(),stop_kx(),kx_step()
        if (numstates(kx_idx)/=0) then
          call file_status_message("    Calling deltafind on kx_idx ", kx_idx)
          call status_message(         "    Calling deltafind on kx_idx ", kx_idx)

          call deltafind(A,kx_idx,ret_delta,retjyup,retjydown,sfe,ret_cooper3)

          delta = delta + ret_delta
          jyup = jyup +retjyup
          jydown=jydown+retjydown
          sfree_loc = sfree_loc + sfe
          cooper3 = cooper3 + ret_cooper3

          call status_message(         "    finished Calling deltafind on kx_idx ", kx_idx)
          call file_status_message("    Finished calling deltafind on kx_idx ", kx_idx)
        end if
      end do
      !+++++

      ! Every rank in this grid already has the same delta
      ! Reduce using rank_comm() to give each rank the full reduction
      call delta%all_reduce(MPI_SUM, rank_comm())
      delta = delta*factors1

      call file_status_message("  Write delta to file")
      if(iam_root())then
        ! text writes
        write(tmp_fname, '(A, I2.2, A)') 'delta-', it, '.txt'
        call delta%write_to_file(tmp_fname)
        call delta%write_to_file(g_txt)
        call delta%write_to_file_with_pos(g_coord_txt)

        ! HDF5 writes
        h5_file = open_hdf5('spin_flip.h5', .true.)
        tmp_fname(:) = ' '
        write(tmp_fname, '(A, I2.2)') 'delta-', it
        call delta%write_to_hdf5(h5_file, trim(tmp_fname))
        call delta%write_to_hdf5(h5_file, 'delta')
        call close_hdf5(h5_file)
      end if
      call file_status_message("  Finished writing delta to file")

      call   jyup%all_reduce(MPI_SUM, rank_comm())
      call jydown%all_reduce(MPI_SUM, rank_comm())

      jyup  =jyup  *(-im*factorsj)
      jydown=jydown*(-im*factorsj)
      jytot=jyup+jydown

      call MPI_Allreduce(sfree_loc,sfree,1,MPI_REAL,&
           MPI_SUM,rank_comm(),ierr)
      call cooper3%all_reduce(MPI_SUM, rank_comm())

      cooper3=cooper3*factors1

      do i = 0,npy-1
         do j = 0,npz-1
            tmp_cmplx = cooper3%get(j,i)
            call phase%set(j,i,atan2(imag(tmp_cmplx),real(tmp_cmplx))/degrees)
         enddo
      enddo

      avg = sum(abs(cooper3%get_section(S1))**2)
      avg = avg + sum(abs(cooper3%get_section(S2))**2)
      
      avg=avg*dY*dZ/(kfl*kfw*lambda)

      sfree=sfree*factors2
      freediff=sfree+avg-nfree

      if(comm_rank == 0) then 
         h5_file = open_hdf5('spin_flip.h5', .true.)         
         
         nphase = phase%get(nz0,R_idx(S2)-npsx-1) - phase%get(nz0,npsx)

         call jyup%write_to_hdf5(h5_file, 'jyup')
         call jydown%write_to_hdf5(h5_file, 'jydown')
         call jytot%write_to_hdf5(h5_file, 'jytot')
         
         tmp_arr = [  real(  jyup%get(zcenter,YC(F1))+    &
                             jydown%get(zcenter,YC(F1))), &
                      real(  jyup%get(zcenter,YC(F2))+    &
                             jydown%get(zcenter,YC(F2))), &
                      real(  jyup%get(zcenter,YC(F3))+    &
                             jydown%get(zcenter,YC(F3)))]

         open(12,position="append",file=avgjcharge_txt)
         write(12,'(i3,2x,3(F10.7,1X),F15.8)') it, tmp_arr, nphase
         close(12)

         call append_avgjcharge_to_hdf5(h5_file, it, tmp_arr, nphase)

         tmp_arr = [nfree, sfree, freediff]

         call append_free_to_hdf5(h5_file, it, tmp_arr)
         call close_hdf5(h5_file)
         
         open(15,position="append",file=free0_txt)
         WRITE (15,'(i3,2x,3(F12.6,1x))') it,nfree,sfree,freediff
         close(15)
      endif

    end if  ! .not. converged

    if (converged) then
      call cooper0%set(cZero)
      call cooper1%set(cZero)
      call cooper2%set(cZero)
      call cooper3%set(cZero)
      call delta%set(cZero)

      call mx%set(cZero)
      call my%set(cZero)
      call mz%set(cZero)

      call   jyup%set(cZero)
      call jydown%set(cZero)
      call   jzup%set(cZero)
      call jzdown%set(cZero)

      call   nup%set(cZero)
      call ndown%set(cZero)

      avg=0
      sfree=0
      sfree_loc=0

      call file_status_message('  Starting calls to finalfind')
      call status_message('  Starting calls to finalfind')
      do kx_idx = start_kx(),stop_kx(),kx_step()
        call file_status_message('    Starting call to finalfind for kx_idx = ', kx_idx)
        call status_message('    Starting call to finalfind for kx_idx = ', kx_idx)
        call finalfind( A, kx_idx, DOS_energies,                &
                        rmx,rmy,rmz,                            &
                        sfe,                                    &
                        ret_cooper0,ret_cooper1,                &
                        ret_cooper2,ret_cooper3,                &
                        retnup,retndown,                        &
                        retjyup,retjydown,                      &
                        retjzup,retjzdown,                      &
                        retjyx,retjyy,retjyz,                   &
                        retjzx,retjzy,retjzz,                   &
                        DOS, ret_delta)

        call file_status_message('    Finished call to finalfind for kx_idx = ', kx_idx)
        cooper0 = cooper0 + ret_cooper0
        cooper1 = cooper1 + ret_cooper1
        cooper2 = cooper2 + ret_cooper2
        cooper3 = cooper3 + ret_cooper3

        delta = delta + ret_delta

        sfree_loc = sfree_loc + sfe
                     
        mx = mx + rmx
        my = my + rmy
        mz = mz + rmz

        jyup = jyup +retjyup
        jydown=jydown+retjydown
        jzup = jzup +retjzup
        jzdown=jzdown+retjzdown

        nup = nup +retnup
        ndown=ndown+retndown

        jyxspin=jyxspin+retjyx
        jyyspin=jyyspin+retjyy
        jyzspin=jyzspin+retjyz
        jzxspin=jzxspin+retjzx
        jzyspin=jzyspin+retjzy
        jzzspin=jzzspin+retjzz
      end do

      !+++++
      call write_eigen_values('W_final.txt', 'W_final_min.txt')

      h5_file = open_hdf5('spin_flip.h5', iam_root())
      call write_eigen_values_to_hdf5(h5_file, 'final_eigen_values')
      call close_hdf5(h5_file)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call file_status_message("  Finished all calls to finalfind")


      call MPI_Allreduce(sfree_loc,sfree,1,MPI_REAL,&
                       MPI_SUM,MPI_COMM_WORLD,ierr)
      sfree = sfree / n_proc_rows

      call cooper0%all_reduce(MPI_SUM, rank_comm())
      
      call jyxspin%all_reduce(MPI_SUM, rank_comm())
      call jyyspin%all_reduce(MPI_SUM, rank_comm())
      call jyzspin%all_reduce(MPI_SUM, rank_comm())

      call jzxspin%all_reduce(MPI_SUM, rank_comm())
      call jzyspin%all_reduce(MPI_SUM, rank_comm())
      call jzzspin%all_reduce(MPI_SUM, rank_comm())

      do k = 1,n_DOS
         call DOS(1,k)%all_reduce(MPI_SUM, rank_comm())
         DOS(1,k) = DOS(1,k) * factors3
         call DOS(2,k)%all_reduce(MPI_SUM, rank_comm())
         DOS(2,k) = DOS(2,k) * factors3
      end do
      !-----	



      do i = 0,npy-1
         do j = 0,npz-1
            tmp_cmplx = cooper3%get(j,i)
            call phase%set(j,i,atan2(imag(tmp_cmplx),real(tmp_cmplx))/degrees)
         enddo
      enddo
      
      ! Find the Cooper pairs in the spin-rotated frame
      ! for the three ferromagnetic regions
      do F_idx = F1,F3
         do i = L_idx(F_idx),R_idx(F_idx)
            do j = 0,npz-1
               ! Find rotated cooper0
               tmp_cmplx = cooper0%get(j,i)*cos(phi(F_idx))*sin(theta(F_idx))     &
                                - cooper1%get(j,i)*cos(theta(F_idx))              &
                                - im*cooper2%get(j,i)*sin(theta(F_idx))*sin(phi(F_idx))

               call cooper0_r%set_at_to_complex(j,i,tmp_cmplx)

               ! Find rotated cooper1
               tmp_cmplx = cooper0%get(j,i)*cos(theta(F_idx))*cos(phi(F_idx))     &
                                + cooper1%get(j,i)*sin(theta(F_idx))              &
                                - im*cooper2%get(j,i)*cos(theta(F_idx))*sin(phi(F_idx))

               call cooper1_r%set_at_to_complex(j,i,tmp_cmplx)

               ! Find rotated cooper2
               tmp_cmplx = cooper2%get(j,i)*cos(phi(F_idx)) &
                           - im*cooper0%get(j,i)*sin(phi(F_idx))

               call cooper2_r%set_at_to_complex(j,i,tmp_cmplx)
            end do
         end do
      end do

      delta = delta*factors1
                    
      nup  =  nup*factorsn
      ndown=ndown*factorsn

      jyxspin=jyxspin*( factorsj)
      jyyspin=jyyspin*(-factorsj)
      jyzspin=jyzspin*( factorsj)
      jzxspin=jzxspin*( factorsj)
      jzyspin=jzyspin*(-factorsj)
      jzzspin=jzzspin*( factorsj)

      mx=mx*(-factors3) !recall normalization divides by  minus(Bohr magneton)
      my=my*(-factors3*im) !recall normalization divides by  minus(Bohr magneton)
      mz=mz*(-factors3) !recall normalization divides by  minus(Bohr magneton)

      tx = ((my%re() * hzvec) - (mz%re() * hyvec)) * (-2.0)
      ty = ((mx%re() * hzvec) - (mz%re() * hxvec)) * ( 2.0)
      tz = ((mx%re() * hyvec) - (my%re() * hxvec)) * (-2.0)
           
      if(comm_rank == 0) then 
        call status_message('Starting final data dump')
        
        h5_file = open_hdf5('spin_flip.h5', .true.)
        write(tmp_fname, '(A, I2.2, A)') 'delta-', it, '.txt'
        call delta%write_to_file(tmp_fname)

        ! HDF5 writes
        write(tmp_fname, '(A, I2.2)') 'delta-', it


        allocate(tmp_arr_2d(n_DOS,4))

        ! Write two forms of DOS output
        do k=1,3           
           do nn=1,n_DOS
              if (abs(DOS_energies(nn) - DOS_energies_out(k)) <= 1e-4) kk = nn
           end do

           ! Write whole domain for selected energies
           ! To DOS
           write(tmp_fname, '(A, f6.3, A)') 'DOS', DOS_energies(kk), '-up'
           tmp_fname = remove_blanks(tmp_fname)
           call DOS(1,kk)%write_to_hdf5(h5_file, tmp_fname)
           write(tmp_fname, '(A, f6.3, A)') 'DOS', DOS_energies(kk), '-down'
           tmp_fname = remove_blanks(tmp_fname)
           call DOS(2,kk)%write_to_hdf5(h5_file, tmp_fname)
           write(tmp_fname, '(A, f6.3, A)') 'DOS', DOS_energies(kk), '-total'
           DOS_tot = DOS(1,kk) + DOS(2,kk)
           tmp_fname = remove_blanks(tmp_fname)
           call DOS_tot%write_to_hdf5(h5_file, tmp_fname)

           ! To .txt
           write(tmp_fname, '(A, f6.3, A)') 'DOS', DOS_energies(kk), '.txt'
           tmp_fname = remove_blanks(tmp_fname)
           open(12,action="write",file=tmp_fname)
           do i=0,npy-1
              do j=0,npz-1
                 DOS_tot = DOS(1,kk) + DOS(2,kk)
                 write(12,'(2F10.4,1X,3(es15.6,1x))') i*dY,j*dZ, &
                      DOS(1,kk)%get(j,i), DOS(2,kk)%get(j,i),DOS_tot%get(j,i)
              end do
           end do
           close(12)

           ! Write energy spectrum for specified point
           write(tmp_fname_2, '(A,f7.2,A,f7.2,A)') 'DOS_y',DOS_loc_y_out(k), &
                "_z",DOS_loc_z_out(k),'.txt'
           tmp_fname_2 = remove_blanks(tmp_fname_2)
           open(13,action='write',file=tmp_fname_2)
           i = nint(grid_int*DOS_loc_y_out(k))
           j = nint(grid_int*DOS_loc_z_out(k))
           do nn=1,n_DOS
              DOS_tot = DOS(1,nn) + DOS(2,nn)
              write(13,'(f6.3,1X,3(es15.6,1x))') DOS_energies(nn), &
                   DOS(1,nn)%get(j,i), DOS(2,nn)%get(j,i),DOS_tot%get(j,i)
              tmp_arr_2d(nn,1) = DOS_energies(nn)
              tmp_arr_2d(nn,2) = DOS(1,nn)%get(j,i)
              tmp_arr_2d(nn,3) = DOS(2,nn)%get(j,i)
              tmp_arr_2d(nn,4) = DOS_tot%get(j,i)
           end do
           close(13)
           write(tmp_fname_2, '(A,f7.2,A,f7.2)') 'DOS_y',DOS_loc_y_out(k), &
                "_z",DOS_loc_z_out(k)
           call write_real_mat_to_hdf5(h5_file,tmp_fname_2,tmp_arr_2d)
        end do !k=1,3

        deallocate(tmp_arr_2d)

        call status_message('Finished writing DOS output')
            
        density = nup + ndown
        call nup%write_to_hdf5(h5_file, 'nup')
        call ndown%write_to_hdf5(h5_file, 'ndown')
        call density%write_to_hdf5(h5_file, 'density')

        ! ------- Statistics (handled in write_to_hdf5
        call cooper0%write_stats(avgf0_txt, maxf0_txt)
        call cooper1%write_stats(avgf1_txt, maxf1_txt)
        call cooper2%write_stats(avgf2_txt, maxf2_txt)
        call cooper3%write_stats(avgf3_txt, maxf3_txt)

        call cooper0_r%write_stats(avgf0rot_txt, maxf0rot_txt)
        call cooper1_r%write_stats(avgf1rot_txt, maxf1rot_txt)
        call cooper2_r%write_stats(avgf2rot_txt, maxf2rot_txt)

        call mx%write_stats(avgmx_txt, maxmx_txt)
        call my%write_stats(avgmy_txt, maxmy_txt)
        call mz%write_stats(avgmz_txt, maxmz_txt)

        ! ------- Area Integrals (handled in write_to_hdf5
        call jyxspin%write_rarea_integrals(syxavg_txt)
        call jyyspin%write_rarea_integrals(syyavg_txt)
        call jyzspin%write_rarea_integrals(syzavg_txt)

        call jzxspin%write_rarea_integrals(szxavg_txt)
        call jzyspin%write_rarea_integrals(szyavg_txt)
        call jzzspin%write_rarea_integrals(szzavg_txt)

        ! ------- Torques
        call jyxspin%write_torque_file(tx_tot_txt, tx)
        call jyyspin%write_torque_file(ty_tot_txt, ty)
        call jyzspin%write_torque_file(tz_tot_txt, tz)

        call jyxspin%write_torque_to_hdf5(h5_file, 'tx_tot', tx)
        call jyyspin%write_torque_to_hdf5(h5_file, 'ty_tot', ty)
        call jyzspin%write_torque_to_hdf5(h5_file, 'tz_tot', tz)

        call status_message('Finished writing Torques to HDF5')

        call jytot%write_cline_integrals(divjc)

        avg = sum(abs(cooper3%get_section(S1))**2)
        avg = avg + sum(abs(cooper3%get_section(S2))**2)

        avg=avg*dY*dZ/(kfl*kfw*lambda)

        sfree=sfree*factors2
        freediff=sfree+avg-nfree

        tmp_arr = [nfree, sfree, freediff]

        call append_free_to_hdf5(h5_file, it, tmp_arr)
        call status_message('Finished writing free0 to HDF5')

        open(15,position="append",file=free0_txt)
        WRITE (15,'(2x,i3,2x,3(F12.6,1x))') it,nfree,sfree,freediff
        close(15)

        open(12,position="append",file=mm_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,'(2F10.4,1X,3(F11.8,1x))') i*dY,j*dZ, &
                    real(mx%get(j,i)),real(my%get(j,i)),real(mz%get(j,i))
          end do
        end do
        write(12, *) ''
        close(12)

        open(12,position="append",file=mm_y)
        do i=0,npy-1
          write(12,73) i*dY, &
                real(mx%get(nz0,i)),real(my%get(nz0,i)),real(mz%get(nz0,i))
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=mm_z)
        do j=0,npz-1
          write(12,73) j*dZ, &
                real(mx%get(j,ny0)),real(my%get(j,ny0)),real(mz%get(j,ny0))
        end do
        write(12,'(A)') ''
        close(12)

        73 format(F10.4,3(1X,f11.8))

        open(12,position="append",file=torque_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,46) i*dY,j*dZ,tx%get(j,i),ty%get(j,i),tz%get(j,i)
          end do
        end do
        write(12,'(A)') ''
        close(12)
        46 format(2(F10.4,1X),3(f13.10,1X))

        open(12,position="append",file=pa_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,78) i*dY,j*dZ, cooper0%get(j,i), &
                        cooper1%get(j,i), cooper2%get(j,i), &
                        cooper3%get(j,i)
          end do
        end do
        write(12,'(A)') ''
        close(12)
        78 format(2(F10.4,1X),8(F10.7,1X))
        788 format(2(F10.4,1X),6(F10.7,1X))

        open(12,position="append",file=pa1D_y)
        do i=0,npy-1
          write(12,799) i*dY, cooper0%get(nz0,i), &
                        cooper1%get(nz0,i), cooper2%get(nz0,i), &
                        cooper3%get(nz0,i)
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=pa1D_z)
        do j=0,npz-1
          write(12,799) j*dZ, cooper0%get(j,ny0), &
                        cooper1%get(j,ny0), cooper2%get(j,ny0), &
                        cooper3%get(j,ny0)
        end do
        write(12,'(A)') ''
        close(12)

        ! Spin_rotated Cooper pairs
        open(12,position="append",file=parot_txt)
        do i=0,npy-1
          do j=0,npz-1
             if (i < L_idx(F1) .or. i > R_idx(F3))then
                ! Write the normal Cooper pairs for the semiconductors
                write(12,788) i*dY,j*dZ, cooper0%get(j,i), &
                        cooper1%get(j,i), cooper2%get(j,i)
             else
                ! Write the rotated pairs for the ferromagnets
                write(12,788) i*dY,j*dZ, cooper0_r%get(j,i), &
                        cooper1_r%get(j,i), cooper2_r%get(j,i)
             end if
          end do
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=pa1Drot_y)
        do i=0,npy-1
           if (i < L_idx(F1) .or. i > R_idx(F3))then
              ! Write the normal Cooper pairs for the semiconductors
              write(12,79) i*dY, cooper0%get(nz0,i), &
                        cooper1%get(nz0,i), cooper2%get(nz0,i)
           else
              ! Write the rotated pairs for the ferromagnets
              write(12,79) i*dY, cooper0_r%get(nz0,i), &
                        cooper1_r%get(nz0,i), cooper2_r%get(nz0,i)
           end if
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=pa1Drot_z)
        do j=0,npz-1
          write(12,79) j*dZ, cooper0_r%get(j,ny0), &
                        cooper1_r%get(j,ny0), cooper2_r%get(j,ny0)
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=jspin1D_y)
        do i=0,npy-1
          write(12,79)  i*dY,                                               &
                        real(jyxspin%get(nz0,i)), real(jyyspin%get(nz0,i)), &
                        real(jyzspin%get(nz0,i)), real(jzxspin%get(nz0,i)), &
                        real(jzyspin%get(nz0,i)), real(jzzspin%get(nz0,i))
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=jspin1D_z)
        do j=0,npz-1
          write(12,79)  j*dZ,                                               &
                        real(jyxspin%get(j,ny0)),real(jyyspin%get(j,ny0)),  &
                        real(jyzspin%get(j,ny0)),real(jzxspin%get(j,ny0)),  &
                        real(jzyspin%get(j,ny0)),real(jzzspin%get(j,ny0))
        end do
        write(12,'(A)') ''
        close(12)

        79 format(F10.4,1X,6(F10.7,1X))
        799 format(F10.4,1X,8(F10.7,1X))

        open(12,position="append",file=dens_z)
        do j=0,npz-1
          write(12,52) j*dZ, real(nup%get(j,ny0) + ndown%get(j,ny0))
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=dens_y)
        do i=0,npy-1
          write(12,52) i*dY,real(nup%get(nz0,i) + ndown%get(nz0,i)) 
        end do
        write(12,'(A)') ''
        close(12)

        52 format(F10.4,1X,F10.7)
        522 format(F10.4,1X,F14.7)
        open(12,position="append",file=jc_z)
        do j=0,npz-1
          write(12,51) j*dZ, real(jytot%get(j,ny0)), real(jztot%get(j,ny0))
        end do
        write(12,'(A)') ''
        close(12)

        open(12,position="append",file=jc_y)
        do i=0,npy-1
          write(12,51) i*dY, real(jytot%get(nz0,i)), real(jztot%get(nz0,i))
        end do
        write(12,'(A)') ''
        close(12)

        51 format(F10.4,1X,F10.7,1X,F10.7)

        open(12,position="append",file=density_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,50) i*dY, j*dZ,                                &
                         real(nup%get(j,i) + ndown%get(j,i)),       &
                         real(nup%get(j,i)), real(ndown%get(j,i))
          end do
        end do
        write(12,'(A)') ''
        close(12)
        50 format(2(F10.4,1X),3(F10.7,1X))
        
        open(12,position="append",file=jcharge_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,58) i*dY, j*dZ, real( jytot%get(j,i)), &
                real(  jyup%get(j,i)), real(jydown%get(j,i)), &
                real( jztot%get(j,i)), &
                real(  jzup%get(j,i)), real(jzdown%get(j,i))
          end do
        end do
        write(12,'(A)') ''
        close(12)
        58 format(2(F10.4,1X),6(es15.6,1X))

        open(12,position="append",file=djspin_txt)
        open(13,position="append",file=tx_line_txt)
        divxsum=0
        divysum=0
        divzsum=0

        call status_message('About to calculate divx, divy, and divz')

        allocate(divx(0:npz-2,0:npy-2))
        allocate(divy(0:npz-2,0:npy-2))
        allocate(divz(0:npz-2,0:npy-2))
        do i=0,npy-2
          do j=0,npz-2

            write(12,18)  i*dY,j*dZ,divx(j,i),divy(j,i),divz(j,i)
            if (i<= nf1-2) then
              divxsum=divxsum+divx(j,i)
              divysum=divysum+divy(j,i)
              divzsum=divzsum+divz(j,i)
            end if
            if (j==npz/2) then
              write(13,20) i*dY,tx%get(j,i),ty%get(j,i),tz%get(j,i), &
                                divx(j,i),  divy(j,i),  divz(j,i)
            end if
          end do
        end do

        call status_message('Writing divx, divy, and divz to HDF5')
        call write_real_mat_to_hdf5(h5_file, 'divx', divx)
        call write_real_mat_to_hdf5(h5_file, 'divy', divy)
        call write_real_mat_to_hdf5(h5_file, 'divz', divz)
        call status_message('Finished writing divx, divy, and divz to HDF5')

        tmp_arr = [divxsum*dY*dZ, divysum*dY*dZ, divzsum*dY*dZ]
        call write_real_array_to_hdf5(h5_file, 'dtavg', tmp_arr)

        close(12)
        close(13)
        18 format(2(F10.4,1X),3(F13.10,1X))
        20 format(F10.4,6(1X,F13.10))

        open(15,position="append",file=dtavg_txt)
        WRITE (15,'(3(F13.10,1x))') tmp_arr
        close(15)

        open(12,position="append",file=jspin_txt)
        do i=0,npy-1
          do j=0,npz-1
            write(12,17) i*dY, j*dZ,                                    &
                         real(jyxspin%get(j,i)),real(jyyspin%get(j,i)), &
                         real(jyzspin%get(j,i)),real(jzxspin%get(j,i)), &
                         real(jzyspin%get(j,i)),real(jzzspin%get(j,i))
          end do
        end do
        write(12,'(A)') ''
        close(12)

        nphase = phase%get(nz0,R_idx(S2)-npsx-1) - phase%get(nz0,npsx)

        tmp_arr = [  real(  jyup%get(zcenter,YC(F1))+    &
                            jydown%get(zcenter,YC(F1))), &
                     real(  jyup%get(zcenter,YC(F2))+    &
                            jydown%get(zcenter,YC(F2))), &
                     real(  jyup%get(zcenter,YC(F3))+    &
                            jydown%get(zcenter,YC(F3)))]

        open(12,position="append",file=avgjcharge_txt)
        write(12,599) it,tmp_arr,nphase
        close(12)

        call append_avgjcharge_to_hdf5(h5_file, it, tmp_arr, nphase)

        open(12,position='append', file="phase_y.txt")
        do i = 0,npy-1
           write(12,522) i*dY, phase%get(nz0,i)
        enddo
        close(12)

        599 format (i3,2x,3(F10.7,1X),F15.8)
        17 format(2(F10.4,1X),6(F13.10,1X))

        call close_hdf5(h5_file)
      end if

      call MPI_barrier(MPI_COMM_WORLD,ierr)
      call file_status_message('Performed all final find file writes')
      call status_message('Performed all final find file writes')

      exit  ! break from sc-loop

    end if ! converged == true condition
    
    tmp_rfield = delta%mag()
    diff   = old_delta%abs_diff(delta)
    
    if (npsx > 0) then
       rerror = diff%maximum_w_kfdx() / tmp_rfield%maximum_w_kfdx()
    else
       rerror = diff%maximum() / tmp_rfield%maximum()
    endif
    old_delta = delta

    if (iam_root()) then
      open (45,position="append",file=iter_txt)
      write (45,84) it,rerror !record iteration history
      write (*,*) "  iteration error =", rerror
      close(45)
      84 format(I5,1X,F12.8)

      h5_file = open_hdf5('spin_flip.h5', iam_root())
      call append_iter_to_hdf5(h5_file, it, rerror)
      call close_hdf5(h5_file)
    end if

    if ( tmp_rfield%maximum() < delsmall ) then
      if (iam_root()) then
        write (*,*) "  delta maximum magnitude =", tmp_rfield%maximum()
      end if
      exit  ! break from sc-loop
    endif
  end do ! sc-loop

  call file_status_message("Ending program")
  close(local_log)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call close_eigen_vector_module()

  call MPI_Finalize(ierr)
END PROGRAM BDG
