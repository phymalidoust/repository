module eigen_vector_mod
    use utilities
    use user_input

    implicit none

    private     ! everything private by default
    public      :: local_comm, rank_comm, iam_root, n_grids, grid_idx
    public      :: local_col_comm, grid_row, col_indices
    public      :: local_row_comm, grid_col, row_indices
    public      :: init_eigen_vector_module, close_eigen_vector_module
    public      :: get_local_rows, rank_with_col, grid_root, get_col_idx

    ! work variables
    complex, allocatable   :: pr_work(:)
    integer                :: pr_lwork
    real, allocatable      :: pr_rwork(:)
    integer                :: pr_lrwork
    integer, allocatable   :: pr_iwork(:)
    integer                :: pr_liwork
    integer                :: pr_blacs_context
    integer, parameter     :: pr_dlen = 9
    integer                :: pr_descA(pr_dlen)
    integer                :: pr_descZ(pr_dlen)
    integer                :: pr_iam
    integer                :: pr_block_factor
    logical                :: pr_initialized = .false.

    ! variabes used that user may want access to
    integer     :: pr_nprocs                    ! total number of processors (ranks) in program
    integer     :: pr_n_grids                   ! number of grids to create
    integer     :: pr_grid_size                 ! the total number of ranks assigned to 1 matrix, n_proc_rows * n_proc_cols
    integer     :: pr_grid_N                    ! rank of grid matrix, input
    integer     :: pr_grid_idx                  ! The grid that this rank is assigned to
    integer(kind=4)     :: pr_local_comm        ! MPI communicator created for this grid
    integer(kind=4)     :: pr_local_col_comm    ! MPI communicator created for a grid column
    integer(kind=4)     :: pr_local_row_comm    ! MPI communicator created for a grid column
    integer(kind=4)     :: pr_rank_comm         ! MPI communicator created for each rank r accross grids

    integer     :: pr_n_proc_rows               ! number of processor rows to use on each grid
    integer     :: pr_n_proc_cols               ! number of processor columns to use on each grid
    integer     :: pr_row                       ! grid row this rank is assigned to
    integer     :: pr_col                       ! grid column this rank is assigned to
    integer     :: pr_local_rows                ! number of rows in the local matrix
    integer     :: pr_local_cols                ! number of columsn in the local matrix

    real        :: pr_tolerance               ! input/output
    real        :: pr_orfac                   ! input/output

    logical     :: pr_iam_root = .false.      ! output
    logical     :: pr_grid_root = .false.     ! output

! This is the primary type that procedures outside of the module will deal with
! This utilizes the ScaLAPACK BLACS library. The NxN Cyclic matrix will be broken
! up over (n_proc_row x n_proc_col) ranks in the following manner
! npc => number of processors cols, npr => number of processor rows
! Each cell in the grid below holds a block_factor x block_factor chunk of the matrix
!
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |  1,1  |  1,2  | ... |  1,npc-1  |   1,npc   |  1,1  |  1,2  | ... |  1,npc-1  |  1,npc  |
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |  2,1  |  2,2  | ... |  2,npc-1  |   2,npc   |  2,1  |  2,2  | ... |  2,npc-1  |  2,npc  |
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!       .       .               .          .         .       .               .           . 
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |npr-1,1|npr-1,2| ... |npr-1,npc-1|npr-1,npc|npr-1,1|npr-1,2| ... |npr-1,npc-1|npr-1,npc|
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   | npr,1 | npr,2 | ... | npr,npc-1 | npr,npc | npr,1 | npr,2 | ... | npr,npc-1 | npr ,npc|
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |  1,1  |  1,2  | ... |  1,npc-1  |  1,npc  |  1,1  |  1,2  | ... |  1,npc-1  |  1,npc  |
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |  2,1  |  2,2  | ... |  2,npc-1  |  2,npc  |  2,1  |  2,2  | ... |  2,npc-1  |  2,npc  |
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!       .       .               .          .         .       .               .           . 
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   |npr-1,1|npr-1,2| ... |npr-1,npc-1|npr-1,npc|npr-1,1|npr-1,2| ... |npr-1,npc-1|npr-1,npc|
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!   | npr,1 | npr,2 | ... | npr,npc-1 | npr,npc | npr,1 | npr,2 | ... | npr,npc-1 | npr,npc |
!   +-------+-------+     +-----------+---------+-------+-------+     +-----------+---------+
!
    public Distributed_Matrix
    type Distributed_Matrix
        private
        ! the actual values for the block cyclic matrix
        complex, dimension(:,:), allocatable  :: m_values         ! (pr_local_rows, pr_local_cols)
        integer, dimension(:), allocatable    :: m_failed_vec_idx ! (pr_grid_N)
        integer, dimension(:), allocatable    :: m_icluster       ! (2*pr_grid_size)
        integer, dimension(:), allocatable    :: m_gap            ! (pr_grid_size)
        real   , dimension(:), allocatable    :: m_eigen_values   ! (pr_grid_N)

        integer   :: m_values_found = 0
        integer   :: m_values_computed = 0

        integer   :: m_vectors_found = 0
        integer   :: m_vectors_computed = 0

        integer   :: m_info = -1

        logical   :: m_values_up_to_date = .false.
    contains
        procedure :: n_values_found
        procedure :: n_values_computed

        procedure :: set_to
        procedure :: set_at
        generic   :: set => set_to, set_at

        procedure :: set_local_col_section

        procedure :: get
        procedure :: get_col
        
        procedure :: check_allocated

        ! calculate and return eigen values, has side effects
        procedure :: eigen_values

        ! calculate and return eigen vectors, has side effects
        procedure :: eigen_vectors

        ! return the number of rows in the matrix
        procedure :: rows

        ! write out the array to a txt file
        procedure :: write_to_txt
        procedure :: write_local_to_txt

        ! gather/scatter
        procedure :: distribute
        procedure :: gather
    end type Distributed_Matrix

    interface Distributed_Matrix
        module procedure value_constructor
    end interface

    public  :: have_col
    interface have_col
        module procedure have_col_scalar
        ! module procedure have_col_array
    end interface

    public  :: have_row
    interface have_row
        module procedure have_row_scalar
        module procedure have_row_array
    end interface

    type(Distributed_Matrix)     :: pr_work_matrix

contains
    subroutine pr_check_initialized()
        if(.not. pr_initialized) then
            write(*,*) 'Attemted to use Eigen_vectors module prior to initialization'
            error stop 1
        endif
    end subroutine

    subroutine init_eigen_vector_module(N, n_proc_rows, n_proc_cols, n_calculations, block_factor, clustersize)
        !use mpi

        implicit none
        include "mpif.h"
        ! N is the rank of the matrix to evaluate
        integer, intent(in) :: N
        ! n_proc_rows is the number of grid rows to use to evaluate a matrix
        integer, intent(in) :: n_proc_rows
        ! n_proc_cols is the number of grid columns to use to evaluate a matrix
        integer, intent(in) :: n_proc_cols
        ! n_calculations is the number of different orientations to evaulate (e.g. nkx)
        integer, intent(in):: n_calculations
        ! block_factor is the block_factor factor to use for the 2d block cyclical matrix
        ! distribution
        integer, intent(in), optional :: block_factor
        integer, intent(in), optional :: clustersize

        integer, dimension(:,:) ,allocatable :: usermap

        integer(kind=4) :: rank
        integer :: local_start
        integer :: i, j
        integer :: NP0, MQ0, NN
        integer :: sqnpc
        integer :: NPS
        integer :: ANB
        integer :: INFO
        integer :: nhetrd_lwork 
        integer :: clustersz
        integer :: blocking
        integer(kind=4) :: ierr
        integer(kind=4) :: i4

        ! external function names
        integer :: numroc
        integer :: iceil
        integer :: pjlaenv
        real    :: pslamch

        external blacs_pinfo
        external blacs_setup
        external blacs_get
        external blacs_gridmap
        external blacs_gridinfo
        external blacs_abort
        external numroc
        external iceil
        external pslamch
        external pjlaenv

        ! if already iniialized, do nothing
        if (pr_initialized) then
            return
        end if

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

        blocking = 8
        if(present(block_factor)) blocking = max(1, block_factor)
        clustersz = 0
        if(present(clustersize)) clustersz = max(clustersize - 1, 0)

        pr_grid_N = N
        pr_n_proc_rows = n_proc_rows
        pr_n_proc_cols = n_proc_cols
        pr_grid_size = pr_n_proc_rows * pr_n_proc_cols
        pr_block_factor = blocking

        pr_nprocs = pr_n_proc_rows * pr_n_proc_cols
        call blacs_pinfo( pr_iam, pr_nprocs )
        if( ( pr_nprocs < 1 ) ) then
            pr_nprocs = pr_n_proc_rows * pr_n_proc_cols
            call blacs_setup( pr_iam,  pr_nprocs)
        end if

        pr_iam_root = (pr_iam == 0)

        if(pr_iam_root) then
            if(modulo(N, (n_proc_cols * blocking)) /= 0) then
                write(*,*) "Please ensure that the number of process columns * the block factor is a factor of N"
                write(*,*) "N = ", N, " block factor = ", blocking, " # process columns = ", n_proc_cols
                error stop 2
            endif

            if(modulo(N, (n_proc_rows * blocking)) /= 0) then
                write(*,*) "Please ensure that the number of process rows * the block factor is a factor of N"
                write(*,*) "N = ", N, " block factor = ", blocking, " # process rows = ", n_proc_rows
                error stop 3
            endif

            if(N < 1) then
                write(*,*) "Please provide a positive value for N"
                error stop 4
            endif

            if(mod(pr_nprocs, pr_grid_size) /= 0) then
                write(*,*) "Grid does not tile onto provided number of processors"
                error stop 5
            endif
        endif

        pr_local_rows = N / n_proc_rows
        pr_local_cols = N / n_proc_cols

        ! get the BLACS system context
        call blacs_get(-1, 0, pr_blacs_context)

        ! determine the rank where the local grid starts
        ! local_start is the starting local nkx.
        local_start = pr_iam - mod(pr_iam, pr_grid_size)
        pr_grid_root = (pr_iam - local_start == 0)      ! the root rank for this grid
        pr_grid_idx = (local_start / pr_grid_size)      ! which grid this rank belongs to
        pr_n_grids = pr_nprocs / pr_grid_size

        call file_status_message('grid_idx', pr_grid_idx)
        call file_status_message('nprocs', pr_nprocs)
        call file_status_message('n_grids', pr_n_grids)
        call file_status_message('local_rows', pr_local_rows)
        call file_status_message('local_cols', pr_local_cols)
        call file_status_message('N', pr_grid_N)
        call file_status_message('local_start', local_start)
        call file_status_message('grid_size', pr_grid_size)

        if ((n_calculations == 1) .and. (pr_n_grids /= 1)) error stop "2D case trying to run with nkx > 1"

        ! comm for the local grid
        i4 = pr_grid_idx
        call MPI_Comm_split(MPI_COMM_WORLD, i4, rank, pr_local_comm, ierr)

        ! comm to communicate between grids
        i4 = pr_iam - local_start
        call MPI_Comm_split(MPI_COMM_WORLD, i4, rank, pr_rank_comm, ierr)

        allocate(usermap(0:n_proc_rows-1, 0:n_proc_cols-1),stat=INFO)
        if(INFO /= 0) then
            call file_status_message('Failed to allocate usermap')
            error stop 6
        endif

        ! Create the map for this grid, use row major natural ordering
        do j = 0, n_proc_cols-1
            do i = 0, n_proc_rows-1
                usermap(i, j) = local_start + i * n_proc_cols + j
            end do
        end do

        ! create the grid, pr_blacs_context will be changed from the system context
        ! to a Matrix specific context associated with a grid
        call blacs_gridmap(pr_blacs_context, usermap, pr_n_proc_rows, pr_n_proc_rows, pr_n_proc_cols)

        ! get grid info 
        call blacs_gridinfo(pr_blacs_context, i, j, pr_row, pr_col )

        if(i /= n_proc_rows .or. j /= n_proc_cols) then
            if(pr_grid_root) then
                write(*,*) "BLACS_GRIDINFO did not return what was expected"
                write(*,*) "[i          , j          ] = [", i, ",", j, "]"
                write(*,*) "[n_proc_rows, n_proc_cols] = [", n_proc_rows, ",",  &
                                            n_proc_cols, "]"
            end if
            error stop 7
        end if

        ! Comm to communicate along a given column of a particualr grid
        i4 = local_start + pr_col
        call MPI_Comm_split(MPI_COMM_WORLD, i4, rank, pr_local_col_comm, ierr)

        ! Comm to communicate along a given row of a particualr grid
        i4 = local_start + pr_row
        call MPI_Comm_split(MPI_COMM_WORLD, i4, rank, pr_local_row_comm, ierr)

        ! Array descriptors for use in ScaLAPACK
        call descinit( pr_descA, N, N, pr_block_factor, pr_block_factor,         &
                        0, 0, pr_blacs_context, pr_local_rows, INFO )
        if(INFO /= 0) then
            if(pr_grid_root) then
                write(*,*) "DESCINIT did not successfully exit"
                write(*,*) "INFO = ", INFO
            end if
            error stop 8
        end if

        call descinit( pr_descZ, N, N, pr_block_factor, pr_block_factor,         &
                        0, 0, pr_blacs_context, pr_local_rows, INFO )
        if(INFO /= 0) then
            if(pr_grid_root) then
                write(*,*) "DESCINIT did not successfully exit"
                write(*,*) "INFO = ", INFO
            end if
            error stop 9
        end if

        ! workspace allocation
        !complex
        pr_lwork = N + (pr_local_rows + pr_local_cols + pr_block_factor) * pr_block_factor
        ANB = pjlaenv(pr_blacs_context, 3, 'PCHETTRD', 'L', 0, 0, 0, 0)
        sqnpc = int(SQRT(DBLE(pr_grid_size)))
        NPS = max(numroc(pr_grid_N, 1, 0, 0, sqnpc), 2 * ANB)
        nhetrd_lwork = pr_grid_N + 2 * (ANB + 1) * (4 * NPS + 2) + (NPS + 2) * NPS
        pr_lwork = max(pr_lwork, nhetrd_lwork)
        allocate(pr_work(pr_lwork),stat=INFO)
        if(INFO /= 0) then
            call file_status_message('Failed to allocate pr_work')
            error stop 20
        endif

        !real, NB = block_factor
        NN = max(pr_grid_N, block_factor, 2)
        ! IROFFA = MOD( IA-1, MB_A )
        NP0 = numroc(pr_grid_N, block_factor, 0, 0, n_proc_rows)
        MQ0 = numroc(max(pr_grid_N, block_factor, 2), block_factor, 0, 0, n_proc_cols)
        pr_lrwork = 4 * pr_grid_N + max(5 * NN, NP0 * MQ0) +           &
                    iceil(pr_grid_N, pr_grid_size) * NN   +            &
                    pr_grid_N * clustersz
        allocate(pr_rwork(pr_lrwork),stat=INFO)
        if(INFO /= 0) then
            call file_status_message('Failed to allocate pr_rwork')
            error stop 21
        endif

        !integer
        NP0 = max(pr_grid_N, pr_grid_size + 1, 4)
        pr_liwork = 6 * NP0
        allocate(pr_iwork(pr_liwork),stat=INFO)
        if(INFO /= 0) then
            call file_status_message('Failed to allocate pr_iwork')
            error stop 22
        endif

        pr_tolerance = 2*pslamch(pr_blacs_context, 'U')
        !pr_orfac = 1.0e-4
        pr_orfac = 1.0e-3

        ! create the local MPI intra communicator associated with this particular grid

        pr_initialized = .true.
    end subroutine

    pure logical function have_col_scalar(col_idx)
        integer, intent(in) :: col_idx
        integer :: tmp

        ! map the column to the first set of blocks
        tmp = modulo(col_idx-1, pr_n_proc_cols*pr_block_factor)
        have_col_scalar = (pr_col == (tmp / pr_block_factor))
    end function

    pure logical function have_row_scalar(row_idx)
        integer, intent(in) :: row_idx
        integer :: tmp

        ! map the row to the first set of blocks
        tmp = modulo(row_idx, pr_n_proc_rows*pr_block_factor)
        have_row_scalar = (pr_row == (tmp / pr_block_factor))
    end function

    pure logical function have_row_array(row_idx)
        integer, dimension(:), intent(in)  :: row_idx
        integer, dimension(:), allocatable :: itmp
        logical, dimension(:), allocatable :: btmp
        integer :: i, gl_block_sz, arr_sz

        arr_sz = size(row_idx)
        allocate(itmp(arr_sz), btmp(arr_sz))

        ! map the column to the first set of blocks
        have_row_array = .false.

        gl_block_sz = pr_n_proc_cols * pr_block_factor
        do i = 1,arr_sz
            itmp(i) = modulo(row_idx(i), gl_block_sz)
        end do

        btmp = ((itmp / pr_block_factor) == pr_row)
        have_row_array = any(btmp)
    end function

    ! returns the total number of grids that are in use by BLACS
    function num_of_grids()
        integer :: num_of_grids

        call pr_check_initialized()

        num_of_grids = pr_nprocs / pr_grid_size
    end function

    subroutine close_eigen_vector_module()
        implicit none

        external BLACS_EXIT
        external BLACS_BARRIER

        call pr_check_initialized()

        ! all processes must call
        call BLACS_BARRIER(pr_blacs_context, 'A')
        
        ! close out BLACS framework
        call BLACS_EXIT(1)

        ! mark module as uninitialized
        pr_initialized = .false.
    end subroutine

! Distributed_Matrix procedures
    function value_constructor(value) result(new_matrix)
        complex, intent(in) :: value
        type(Distributed_Matrix) :: new_matrix

        call new_matrix%check_allocated()

        new_matrix%m_values = value
        new_matrix%m_values_up_to_date = .false.
    end function

    subroutine check_allocated(this)
        class(Distributed_Matrix), intent(inout) :: this

        integer :: ierr

        call pr_check_initialized()

        if(.not. allocated(this%m_values)) then
            allocate(this%m_values(pr_local_rows, pr_local_cols), stat=ierr)
            if(ierr /= 0) then
                call file_status_message('Failed to allocate m_values')
                error stop 31
            endif

            allocate(this%m_failed_vec_idx (pr_grid_N),stat=ierr)
            if(ierr /= 0) then
                call file_status_message('Failed to allocate m_failed_vec_idx')
                error stop 32
            endif

            allocate(this%m_icluster(2*pr_grid_size),stat=ierr)
            if(ierr /= 0) then
                call file_status_message('Failed to allocate m_icluster, size =', 2*pr_grid_size)
                error stop 33
            endif

            allocate(this%m_gap(pr_grid_size),stat=ierr)
            if(ierr /= 0) then
                call file_status_message('Failed to allocate m_gap, size =', pr_grid_size)
                error stop 34
            endif

            allocate(this%m_eigen_values(pr_grid_N),stat=ierr)
            if(ierr /= 0) then
                call file_status_message('Failed to allocate m_eigen_values')
                error stop 35
            endif
        endif
    end subroutine

    function eigen_vectors(this, IL, IU)
        !use mpi

        implicit none
        include "mpif.h"
        class(Distributed_Matrix), intent(inout)  :: this   ! 2d block cyclic input matrix
                                                            ! m_values will not be overwritten
        integer, intent(in),optional :: IL, IU              ! if present represent lower and
                                                            ! upper indices of eigen{values|vectors}
                                                            ! to calculate
        integer                      :: IL_, IU_            ! varaibles that can be changed
        type(Distributed_Matrix) :: eigen_vectors           ! 2d block cyclic output matrix
                                                            ! contains eigen vectors of A

        integer                 :: N_A                      ! variable to pass into unused
                                                            ! subroutine parameters

        character               :: I_or_A

        external BLACS_BARRIER

        ! make a copy of this, so that when PCHEEVX gets called, it doesn't destroy what
        ! was held within this

        N_A = pr_grid_N

        I_or_A = 'A'
        IL_ = 1
        IU_ = pr_grid_N
        if(present(IL) .or. present(IU)) then
            if(present(IL)) then
                IL_ = IL
            endif
            if(present(IU)) then
                IU_ = IU
            endif
            if(IU_ < IL_) then
                write(*,*) 'IU must not be less than IL', IL_, IU_
                error stop 101
            endif
            if(IL_ < 1 .or. IL_ > pr_grid_N) then
                write(*,*) 'IL is out of bounds'
                error stop 102
            endif
            if(IU_ < 1 .or. IU_ > pr_grid_N) then
                write(*,*) 'IL is out of bounds'
                error stop 102
            endif
            I_or_A = 'I'
        endif

        call this%check_allocated()
        call pr_work_matrix%check_allocated()

        pr_work_matrix%m_values = this%m_values

        eigen_vectors = Distributed_Matrix((0.0, 0.0))

        call BLACS_BARRIER(pr_blacs_context, 'A')

        this%m_values_found = 0
        this%m_vectors_computed = 0
        this%m_eigen_values = 0
        this%m_failed_vec_idx = 0
        this%m_icluster = 0
        this%m_gap = 0

        call file_status_message("    Calling PCHEEVX")
        call PCHEEVX(   'V', I_or_A, 'U',                                           &
                        pr_grid_N, pr_work_matrix%m_values, 1, 1, pr_descA,         &
                        IL_, IU_, N_A, N_A, pr_tolerance,                           &
                        this%m_values_found, eigen_vectors%m_vectors_computed,      &
                        this%m_eigen_values, pr_orfac, eigen_vectors%m_values,      &
                        1, 1, pr_descZ,                                             &
                        pr_work, pr_lwork, pr_rwork, pr_lrwork, pr_iwork, pr_liwork,&
                        this%m_failed_vec_idx, this%m_icluster, this%m_gap,         &
                        this%m_info)
        call file_status_message("    Finished Calling PCHEEVX")

        if(this%m_info /= 0) then
            call file_status_message("Call to PCHEEVX resulted in info =", this%m_info)
            call file_status_message("icluster")
            write(local_log,*)  this%m_icluster
            flush(local_log)
        endif

        eigen_vectors%m_eigen_values    = this%m_eigen_values
        eigen_vectors%m_values_found    = this%m_values_found
        eigen_vectors%m_gap             = this%m_gap
        eigen_vectors%m_icluster        = this%m_icluster
        eigen_vectors%m_failed_vec_idx  = this%m_failed_vec_idx
        eigen_vectors%m_info            = this%m_info

        this%m_values_up_to_date        = .true.
        eigen_vectors%m_values_up_to_date = .false.
    end function

    function eigen_values(this)
        implicit none
        class(Distributed_Matrix), intent(inout)    :: this   ! 2d block cyclic input matrix
                                                              ! m_values will not be overwritten
        integer                 :: N_A = 0      ! variable to pass into unused
                                                ! subroutine parameters
        real, dimension(:), allocatable  :: eigen_values

        complex,dimension(pr_local_rows, pr_local_cols) :: Z

        external BLACS_BARRIER

        call this%check_allocated()

        ! if the eigen values are all up to date, don't have to call PCHEEVX
        if(this%m_values_up_to_date) then
            eigen_values = this%m_eigen_values
            return
        endif

        ! make a copy of this, so that when PCHEEVX gets called, it doesn't destroy what
        ! was held within this
        pr_work_matrix%m_values = this%m_values

        call BLACS_BARRIER(pr_blacs_context, 'A')
        Z = 0

        call file_status_message('    Calling PCHEEVX')
        call PCHEEVX(   'N', 'A', 'U',                                              &
                        pr_grid_N, pr_work_matrix%m_values, 1, 1, pr_descA,         &
                        N_A, N_A, N_A, N_A,                                         &
                        pr_tolerance, this%m_values_found, N_A,                     &
                        this%m_eigen_values, pr_orfac, pr_work_matrix%m_values,     &
                        1, 1, pr_descZ,                                             &
                        pr_work, pr_lwork, pr_rwork, pr_lrwork, pr_iwork, pr_liwork,&
                        this%m_failed_vec_idx, this%m_icluster, this%m_gap, this%m_info)
        call file_status_message('    finished Calling PCHEEVX')

        eigen_values = this%m_eigen_values
        this%m_values_up_to_date = .true.
    end function

    pure integer function rows(this)
        class(Distributed_Matrix), intent(in)    :: this   ! 2d block cyclic input matrix

        rows = pr_grid_N
    end function

    pure integer function n_values_found(this)
        class(Distributed_Matrix), intent(in)  :: this

        n_values_found = this%m_values_found
    end function

    pure integer function n_values_computed(this)
        class(Distributed_Matrix), intent(in)  :: this

        n_values_computed = this%m_values_computed
    end function

    subroutine set_to(this, value)
        class(Distributed_Matrix), intent(inout)   :: this
        complex, intent(in) :: value

        call this%check_allocated()

        this%m_values_up_to_date = .false.
        this%m_values = value
    end subroutine

    subroutine set_at(this, global_i, global_j, value)
        class(Distributed_Matrix), intent(inout)   :: this
        integer, intent(in) :: global_i, global_j
        complex, intent(in) :: value
        integer, dimension(2) :: local_idx

        if(.not. have_row(global_i) .or. .not. have_col(global_j)) then
            return
        endif

        local_idx = get_local_idx(global_i, global_j)

        this%m_values_up_to_date = .false.
        this%m_values(local_idx(1), local_idx(2)) = value
    end subroutine

    pure complex function get(this, global_i, global_j)
        class(Distributed_Matrix), intent(in)  :: this
        integer, intent(in) :: global_i, global_j
        integer, dimension(2) :: local_idx

        ! return 0.0 if we don't have this row or column
        if(.not. have_row(global_i) .or. .not. have_col(global_j)) then
            get = (0.0, 0.0)
            return
        endif

        local_idx = get_local_idx(global_i, global_j)
        get = this%m_values(local_idx(1), local_idx(2))
    end function

    function get_col(this, idx)
        class(Distributed_Matrix), intent(in), target :: this
        integer, intent(in)     :: idx
        complex, dimension(:), pointer :: get_col
        integer :: col_idx(2)
        
        col_idx = get_col_idx(idx)

        get_col => this%m_values(:,col_idx(2))
    end function

    pure function get_local_idx(global_i, global_j) result(ans)
        integer, intent(in) :: global_i, global_j
        integer, dimension(2) :: ans
        integer :: i, j             ! 0-based array idx
        integer :: block_number, block_idx

        i = global_i - 1
        j = global_j - 1
        block_number = i / (pr_n_proc_rows * pr_block_factor)
        block_idx = modulo(i, pr_block_factor)
        ans(1) = block_number * pr_block_factor + block_idx + 1

        block_number = j / (pr_n_proc_cols * pr_block_factor)
        block_idx = modulo(j, pr_block_factor)
        ans(2) = block_number * pr_block_factor + block_idx + 1
    end function

    pure function get_col_idx(global_j) result(ans)
        integer, intent(in) :: global_j
        integer :: ans
        integer :: j             ! 0-based array idx
        integer :: block_number, block_idx

        j = global_j - 1
        block_number = j / (pr_n_proc_cols * pr_block_factor)
        block_idx = modulo(j, pr_block_factor)
        ans = block_number * pr_block_factor + block_idx + 1
    end function


    pure function get_local_rows(global_rows) result(ans)
        integer, dimension(:), intent(in) :: global_rows
        integer, dimension(:), allocatable :: ans
        integer, dimension(:), allocatable :: block_number, block_idx

        block_number = global_rows / (pr_n_proc_rows * pr_block_factor)
        block_idx    = global_rows - (block_number * (pr_n_proc_rows * pr_block_factor))
        ans = (block_number * pr_block_factor) + block_idx
    end function

    pure integer(kind=4) function local_comm()
        local_comm = pr_local_comm
    end function

    pure integer(kind=4) function local_row_comm()
        local_row_comm = pr_local_row_comm
    end function

    pure integer(kind=4) function local_col_comm()
        local_col_comm = pr_local_col_comm
    end function

    pure integer(kind=4) function rank_comm()
        rank_comm = pr_rank_comm
    end function

    pure logical function iam_root()
        iam_root = pr_iam_root
    end function

    pure integer function n_grids()
        n_grids = pr_n_grids
    end function

    pure integer function grid_idx()
        grid_idx = pr_grid_idx
    end function

    pure logical function grid_root()
        grid_root = pr_grid_root
    end function

    pure integer function grid_col()
        grid_col = pr_col
    end function

    pure integer function grid_row()
        grid_row = pr_row
    end function

    function col_indices()
        integer,dimension(pr_local_cols) :: col_indices
        integer :: n_blocks, i_block, i_idx, k, offset

        n_blocks = pr_local_cols / pr_block_factor
        k = 1
        do i_block = 0,n_blocks-1       ! 0-indexed for easier math
            do i_idx = 1,pr_block_factor
                offset = i_block * pr_n_proc_cols * pr_block_factor + pr_col * pr_block_factor
                col_indices(k) = offset + i_idx
                k = k+1
            end do
        end do
    end function

    function row_indices()
        integer,dimension(pr_local_rows) :: row_indices
        integer :: n_blocks, i_block, i_idx, k, offset

        n_blocks = pr_local_rows / pr_block_factor
        k = 1
        do i_block = 0,n_blocks-1       ! 0-indexed for easier math
            do i_idx = 1,pr_block_factor
                offset = i_block * pr_n_proc_rows * pr_block_factor + pr_row * pr_block_factor
                row_indices(k) = offset + i_idx
                k = k+1
            end do
        end do
    end function


    pure function rank_with_col(n)
        integer, intent(in) :: n
        integer :: rank_with_col

        ! shift down to first set of blocks
        rank_with_col = modulo(n, pr_block_factor * pr_n_proc_cols)

        ! determine which block (i.e. rank) it's in
        rank_with_col = rank_with_col / pr_block_factor
    end function

    subroutine distribute(this, A, root)
        !use mpi
        include "mpif.h"
        class(Distributed_Matrix), intent(inout)  :: this
        complex, dimension(:,:), intent(in)       :: A
        integer, intent(in)                       :: root
        integer :: global_k, local_k, ierr

        call this%check_allocated()

        global_k = 1
        local_k = 1

        do while(global_k < pr_grid_N)
            call MPI_Scatter(A(1,global_k),                             &
                             pr_grid_N*pr_block_factor, MPI_COMPLEX,    &
                             this%m_values(1,local_k),                  &
                             pr_grid_N*pr_block_factor, MPI_COMPLEX,    &
                             root, pr_local_comm, ierr)
            if(ierr /= 0) then
                write(*,*) 'Error with MPI_Scatter'
                call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
            endif

            global_k = global_k + (pr_block_factor*pr_n_proc_cols)
            local_k = local_k + pr_block_factor
        enddo

        this%m_values_up_to_date = .false.
    end subroutine

    function gather(this, root)
        !use mpi
        include "mpif.h"
        class(Distributed_Matrix), intent(inout)    :: this
        integer, intent(in)                         :: root
        complex, dimension(:,:), allocatable        :: gather

        integer     :: global_k, local_k
        integer     :: local_rank, ierr

        call MPI_Comm_rank(pr_local_comm, local_rank, ierr)
        if(local_rank == root) then
            allocate(gather(pr_grid_N, pr_grid_N))
        else
            allocate(gather(1,1))
        endif

        global_k = 1
        local_k = 1

        do while(global_k < pr_grid_N)
            if(local_rank == root) then
              call MPI_Gather(this%m_values(1,local_k),                     &
                              pr_grid_N*pr_block_factor, MPI_COMPLEX,       &
                              gather(1,global_k),                           &
                              pr_grid_N*pr_block_factor, MPI_COMPLEX,       &
                              0, pr_local_comm, ierr)
            else
              call MPI_Gather(this%m_values(1,local_k),                     &
                              pr_grid_N*pr_block_factor, MPI_COMPLEX,       &
                              gather(1,1),                                  &
                              pr_grid_N*pr_block_factor, MPI_COMPLEX,       &
                              0, pr_local_comm, ierr)
            endif

            global_k = global_k + (pr_block_factor*pr_n_proc_cols)
            local_k = local_k + pr_block_factor
        end do
    end function

    subroutine write_to_txt(this, filename)
        class(Distributed_Matrix), intent(inout)    :: this
        character(len=*), intent(in)                :: filename
        complex, dimension(:,:), allocatable        :: A
        integer :: k
        character(99)       :: fmt_

        ! gather the distrubted matrix to the comm root
        A = this%gather(0)

        if(pr_grid_root) then
            open(123, file=filename)

            write(fmt_, '(A,i5,A)') '(', pr_grid_N, "('[',e12.5,';',e12.5,'],'))"
        
            do k=1,pr_grid_N
                write(123,fmt_) A(k,:)
            end do
            close(123)
            write(*,*) 'A(1,3073) =', A(1,3073), this%get(1,3073)
        endif
    end subroutine

    subroutine write_local_to_txt(this, filenum)
        class(Distributed_Matrix), intent(in)   :: this
        integer, intent(in) :: filenum
        integer             :: i,j

        do i=1,pr_local_cols
            do j=1,pr_local_rows
                write(filenum, '(e12.5, 1x, e12.5)') this%m_values(j,i)
            end do
        end do
    end subroutine

!    subroutine write_to_h5(this, filename)
!        use hdf5
!        use mpi
!        implicit none
!        
!        class(Distributed_Matrix), intent(in)   :: this
!        character(len=*),intent(in)             :: filename
!
!        integer(HID_T) :: file_id       ! File identifier
!        integer(HID_T) :: dset_id       ! Dataset identifier
!        integer(HID_T) :: filespace     ! Dataspace identifier in file
!        integer(HID_T) :: plist_id      ! Property list identifier
!
!        integer(HSIZE_T), DIMENSION(2) :: dimsf = ([5,8]) ! Dataset dimensions.
!
!        integer         :: ierr
!        integer         :: comm, rank
!        integer         :: info
!
!        call h5open_f(ierr)
!
!        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
!        call h5fcreate_f(filename
!
!    end subroutine
end module eigen_vector_mod
