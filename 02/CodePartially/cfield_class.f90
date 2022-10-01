

module cfield_class
  use field_constants
  use user_input
  use rfield_class

  implicit none
  private                           ! hide everything by default

  public :: cfield
  type cfield
    private
    complex, dimension(:,:), allocatable :: m_values
  contains
    procedure   :: set_to_real
    procedure   :: set_to_complex
    procedure   :: set_at_to_real
    procedure   :: set_at_to_complex
    generic     :: set => set_to_real, set_to_complex, set_at_to_real, set_at_to_complex

    procedure   :: set_section_to_real      ! set all of a section to a given real value
    procedure   :: set_section_to_complex   ! set all of a section to a given complex value
    procedure   :: set_section_to_rarray    ! set all of a section to the real values in an array
    procedure   :: set_section_to_carray    ! set all of a section to the complex values in an array
    generic     :: set_section =>   set_section_to_real, set_section_to_complex, &
                                    set_section_to_rarray, set_section_to_carray

    procedure   :: get                      ! returns an individual element
    procedure   :: get_section              ! returns an array with the values in that section

    procedure   :: inc_at
    procedure   :: inc_by_real
    procedure   :: inc_by_complex
    procedure   :: inc_by_rfield
    procedure   :: inc_by_cfield
    generic     :: operator(+) => inc_by_real,inc_by_complex,inc_by_rfield,inc_by_cfield

    procedure   :: dec_by_rfield
    procedure   :: dec_by_cfield
    generic     :: operator(-) => dec_by_rfield,dec_by_cfield

    procedure   :: mult_by_real
    procedure   :: mult_by_complex
    procedure   :: mult_by_rfield
    procedure   :: mult_by_cfield
    generic     :: operator(*) => mult_by_real, mult_by_complex, mult_by_rfield, mult_by_cfield

    procedure   :: mag                      ! absolute value, different name to prevent collision
    procedure   :: abs_diff                 ! ||A| - |B||

    ! get portions of the answer
    procedure   :: con                      ! conjugate, returns a cfield
    procedure   :: re                       ! real part, returns an rfield
    procedure   :: i_                       ! imaginary part, returns an rfield

    ! I/O subroutines
    procedure   :: write_to_file
    procedure   :: write_to_file_with_pos
    procedure   :: read_from_file

    procedure :: write_stats
    procedure :: write_rarea_integrals
    procedure :: write_cline_integrals
    procedure :: write_torque_file

    procedure   :: write_to_hdf5
    procedure   :: write_torque_to_hdf5

    procedure :: rline_integrals

    ! MPI procedures
    procedure :: all_reduce
    procedure :: broadcast
  end type cfield

  interface cfield
    procedure complex_value_constructor
    procedure real_value_constructor
    procedure complex_array_constructor
    procedure copy_constructor
    procedure copy_rfield_constructor
  end interface

contains
! constructor sets the values elements to the values of the array
  function complex_array_constructor(values) result(new_field)
    complex, intent(in), dimension(0:npz-1,L_idx(S1):R_idx(S2)) :: values
    type(cfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = values
  end function

! constructor sets the values of all elements to the intial value
  function complex_value_constructor(initial_value) result(new_field)
    complex, intent(in) :: initial_value
    type(cfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = initial_value
  end function

! constructor sets the values of all elements to the intial value
  function real_value_constructor(initial_value) result(new_field)
    real, intent(in) :: initial_value
    type(cfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = initial_value
  end function

! constructor creates a copy of the provided cfield
  function copy_constructor(old_field) result(new_field)
    type(cfield), intent(in) :: old_field
    type(cfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = old_field%m_values
  end function

! constructor creates a copy of the provided rfield
  function copy_rfield_constructor(old_field) result(new_field)
    type(rfield), intent(in) :: old_field
    type(cfield)             :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = old_field%to_array()
  end function

! sets the value of a specific element
  subroutine set_at_to_real(this, z_idx, y_idx, value)
    class(cfield), intent(inout) :: this
    integer, intent(in) :: z_idx, y_idx
    real, intent(in) :: value

    this%m_values(z_idx, y_idx) = value
  end subroutine

! sets the value of a specific element
  subroutine set_at_to_complex(this, z_idx, y_idx, value)
    class(cfield), intent(inout) :: this
    integer, intent(in) :: z_idx, y_idx
    complex, intent(in) :: value

    this%m_values(z_idx, y_idx) = value
  end subroutine

! sets all of the elemnts to the given value
  subroutine set_to_real(this, value)
    class(cfield), intent(inout) :: this
    real, intent(in) :: value

    this%m_values = value
  end subroutine

! sets all of the elemnts to the given value
  subroutine set_to_complex(this, value)
    class(cfield), intent(inout) :: this
    complex, intent(in) :: value

    this%m_values = value
  end subroutine

! sets all of the elements in a section to a real value
  subroutine set_section_to_real(this, s_idx, value)
    class(cfield), intent(inout)    :: this
    integer, intent(in)             :: s_idx
    real, intent(in)                :: value

    this%m_values(:,L_idx(s_idx):R_idx(s_idx)) = value
  end subroutine

! sets all of the elements in a section to a complex value
  subroutine set_section_to_complex(this, s_idx, value)
    class(cfield), intent(inout)    :: this
    integer, intent(in)             :: s_idx
    complex, intent(in)             :: value

    this%m_values(:,L_idx(s_idx):R_idx(s_idx)) = value
  end subroutine

! sets a section to the corresponding values in a real array
  subroutine set_section_to_rarray(this, s_idx, arr)
    class(cfield), intent(inout)        :: this
    integer, intent(in)                 :: s_idx
    real, dimension(:,:), intent(in)    :: arr

    this%m_values(:,L_idx(s_idx):R_idx(s_idx)) = arr
  end subroutine

! sets a section to the corresponding values in a complex array
  subroutine set_section_to_carray(this, s_idx, arr)
    class(cfield), intent(inout)        :: this
    integer, intent(in)                 :: s_idx
    complex, dimension(:,:), intent(in) :: arr

    this%m_values(:,L_idx(s_idx):R_idx(s_idx)) = arr
  end subroutine

! returns the value at a specfic element
  function get(this, z_idx, y_idx) result(value)
    class(cfield), intent(in) :: this
    integer, intent(in) :: z_idx, y_idx
    complex :: value

    value = this%m_values(z_idx, y_idx)
  end function

! get an array that for a specific section
  function get_section(this, s_idx)
    class(cfield), target, intent(in) :: this
    integer, intent(in) :: s_idx
    complex, dimension(:,:), pointer :: get_section

    get_section => this%m_values(:,L_idx(s_idx):R_idx(s_idx))
  end function

! inc the element at z_idx,y_idx by the provided amount
  subroutine inc_at(this, z_idx, y_idx, value)
    class(cfield), intent(inout) :: this
    complex, intent(in) :: value
    integer, intent(in) :: z_idx, y_idx

    this%m_values(z_idx, y_idx) = this%m_values(z_idx, y_idx) + value
  end subroutine

! inc all elements by the given scalar value, return answer
  function inc_by_real(this, value) result(ans)
    class(cfield), intent(in) :: this
    real, intent(in) :: value
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values + value
  end function

! inc all elements by the given scalar value, return answer
  function inc_by_complex(this, value) result(ans)
    class(cfield), intent(in) :: this
    complex, intent(in) :: value
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values + value
  end function

! inc each element by the corresponding element in the addend
  function inc_by_rfield(this, addend) result(ans)
    class(cfield), intent(in) :: this
    class(rfield), intent(in) :: addend
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values + addend%to_array()
  end function

! inc each element by the corresponding element in the addend
  function inc_by_cfield(this, addend) result(ans)
    class(cfield), intent(in) :: this
    class(cfield), intent(in) :: addend
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values + addend%m_values
  end function

! decrement each element by the corresponding element in the argument
  function dec_by_rfield(this, subtrahend) result(ans)
    class(cfield), intent(in) :: this
    class(rfield), intent(in) :: subtrahend
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values - subtrahend%to_array()
  end function

! decrement each element by the corresponding element in the argument
  function dec_by_cfield(this, subtrahend) result(ans)
    class(cfield), intent(in) :: this
    class(cfield), intent(in) :: subtrahend
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values - subtrahend%m_values
  end function

! multiply all elements by the given scalar value, return answer
  function mult_by_real(this, value) result(ans)
    class(cfield), intent(in) :: this
    real, intent(in) :: value
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values * value
  end function

! multiply all elements by the given scalar value, return answer
  function mult_by_complex(this, value) result(ans)
    class(cfield), intent(in) :: this
    complex, intent(in) :: value
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values * value
  end function

! multiply two fields, element-by-element
  function mult_by_rfield(this, factor) result(ans)
    class(cfield), intent(in) :: this
    class(rfield), intent(in) :: factor
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values * factor%to_array()
  end function

! multiply two fields, element-by-element
  function mult_by_cfield(this, factor) result(ans)
    class(cfield), intent(in) :: this
    class(cfield), intent(in) :: factor
    type(cfield)  :: ans

    allocate(ans%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    ans%m_values = this%m_values * factor%m_values
  end function

! the magnitude (absolute value) of each of the elements
  function mag(this)
    class(cfield), intent(in)   :: this
    type(rfield)                :: mag
    
    mag = rfield(abs(this%m_values))
  end function

! ||this| - |that||, where |X| => abs(X)
  function abs_diff(this, that)
    class(cfield), intent(in)           :: this
    class(cfield), intent(in)           :: that
    type(rfield)                       :: abs_diff

    abs_diff = this%mag() - that%mag()
    abs_diff = abs_diff%mag()
  end function

! return the conjugate of the values
  function con(this) result(ans)
    class(cfield), intent(in) :: this
    type(cfield)              :: ans

    ans = cfield(conjg(this%m_values))
  end function

! return the real portion of the values
  function re(this) result(ans)
    class(cfield), intent(in) :: this
    type(rfield)              :: ans

    ans = rfield(real(this%m_values))
  end function

! return the imaginary portion of the values
  function i_(this) result(ans)
    class(cfield), intent(in) :: this
    type(rfield)              :: ans

    ans = rfield(imag(this%m_values))
  end function

! write out the cfield
  subroutine write_to_file(this, filename)
    class(cfield), intent(in)       :: this
    character(len=*), intent(in)    :: filename
    integer :: z_idx, y_idx

    open(88, file=filename, action="write")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        write(88,'(2(es15.7,1x))') this%m_values(z_idx, y_idx)
      end do
    end do
    close(88)
  end subroutine

  subroutine write_torque_to_hdf5(this, loc_id, field_name, torque_field)
    use hdf5_utilities
    use hdf5

    class(cfield), intent(in)    :: this
    class(rfield), intent(in),target    :: torque_field
    character(len=*), intent(in) :: field_name
    integer(hid_t), intent(in)   :: loc_id

    real, dimension(F1:F3+3)  :: torques
    integer                 :: s_idx

    do s_idx = F1,F3
      torques(s_idx) = sum(torque_field%get_section(s_idx)) * dA
    end do

    torques(F3+1:F3+3) = this%rline_integrals()

    call write_real_array_to_hdf5(loc_id, field_name, torques)
  end subroutine

  subroutine write_to_hdf5(this, loc_id, field_name)
    use hdf5_utilities
    use hdf5

    class(cfield), intent(in)       :: this
    character(len=*), intent(in)    :: field_name
    integer(hid_t), intent(in)      :: loc_id
    integer(hsize_t), dimension(3)  :: dims
    integer                         :: ierr, s_idx, i
    integer(hid_t)                  :: dset_id, filespace, grp_id
    real, dimension(npz, npy, 2)    :: scratch

    logical                         :: already_exists
    real, dimension(5)              :: averages, maxes, integrals
    complex, dimension(F1:F3)       :: clines
    complex, dimension(:,:), pointer :: section
    complex, dimension(0:npy-1)     :: y_scan
    complex, dimension(0:npz-1)     :: z_scan
    complex                         :: left, right

    ! the actual values
    dims = [npz, npy, 2]
    scratch(:,:,1) = real(this%m_values)
    scratch(:,:,2) = imag(this%m_values)
    call H5Lexists_f(loc_id, field_name, already_exists, ierr)

    if(already_exists) then
        call H5Gopen_f(loc_id, field_name, grp_id, ierr)

        call H5Dopen_f(grp_id, 'values', dset_id, ierr)
        call H5Dget_space_f(dset_id, filespace, ierr)
    else
        call H5Gcreate_f(loc_id, field_name, grp_id, ierr)
        call H5Screate_simple_f(3, dims, filespace, ierr, dims)
        call H5Dcreate_f(grp_id, 'values', H5T_NATIVE_REAL, filespace, &
                        dset_id, ierr)
    endif

    call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)
    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(filespace, ierr)

    ! Section based values
    do s_idx = S1,S2
      section => this%get_section(s_idx)
      averages(s_idx) = sum(abs(section))*dA/Area(s_idx)
      maxes(s_idx)    = maxval(abs(section))
      integrals(s_idx) = sum(real(section))*dA   ! S1
    end do

    call write_real_array_to_hdf5(grp_id, 'averages', averages)
    call write_real_array_to_hdf5(grp_id, 'maxes', maxes)
    call write_real_array_to_hdf5(grp_id, 'area_integrals', integrals)

    ! line integrals
!   rlines = this%rline_integrals()
!   call write_real_array_to_hdf5(grp_id, 'r_line_integrals', rlines)

    do s_idx=F1,F3
      left  = sum(this%m_values(:,L_idx(s_idx)))*dZ
      right = sum(this%m_values(:,R_idx(s_idx)))*dZ
      clines(s_idx) = right - left
    end do

    call write_complex_array_to_hdf5(grp_id, 'c_line_integrals', clines)

    ! centerline scans
    do i = 0,npy-1
        y_scan(i) = this%get(nz0,i)
    enddo
    call write_complex_array_to_hdf5(grp_id, 'y_scan', y_scan)

    do i = 0,npz-1
        z_scan(i) = this%get(i,ny0)
    enddo
    call write_complex_array_to_hdf5(grp_id, 'z_scan', z_scan)
  end subroutine

! write out the cfield, include the position
  subroutine write_to_file_with_pos(this, filename)
    class(cfield), intent(in)       :: this
    character(len=*), intent(in)    :: filename
    integer :: z_idx, y_idx

    open(88, file=filename, action="write")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        write(88,'(2(f10.5,1x), 2(f10.7,1x))') y_idx*dY, z_idx*dZ,        &
                real(this%m_values(z_idx, y_idx)), imag(this%m_values(z_idx, y_idx))
      end do
    end do
    close(88)
  end subroutine

! read in the cfield from a file
  subroutine read_from_file(this, filename)
    class(cfield), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: z_idx, y_idx
    real    :: re_, im_

    open(88, file=filename, action="READ")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        read(88,'(2(e15.7,1x))') re_, im_
        this%m_values(z_idx, y_idx) = cmplx(re_, im_)
      end do
    end do
    close(88)
  end subroutine

! write out some stats about the field
  subroutine write_stats(this, avg_filename, max_filename)
    class(cfield), intent(in), target :: this
    real, dimension(5)  :: averages
    real, dimension(5)  :: maxes
    integer             :: s_idx
    character(len=*), intent(in) :: avg_filename
    character(len=*), intent(in) :: max_filename
    complex, dimension(:,:), pointer :: section

    do s_idx = S1,S2
      section => this%get_section(s_idx)
      averages(s_idx) = sum(abs(section))*dA/Area(s_idx)
      maxes(s_idx)    = maxval(abs(section))
    end do

    open(15, position="append", file=avg_filename)
    write(15, '(5(F10.7,1x))') averages(F1:F3), averages(S1), averages(S2)
    close(15)

    open(15, position="append", file=max_filename)
    write(15, '(5(F10.7,1x))')  maxes(F1:F3), maxes(S1), maxes(S2)
    close(15)
  end subroutine

! write out the area integrals of all of the sections
  subroutine write_rarea_integrals(this, filename)
    class(cfield), intent(in), target :: this
    character(len=*), intent(in) :: filename
    real, dimension(5)  :: integrals
    integer             :: s_idx
    complex, dimension(:,:), pointer :: section

    do s_idx = S1,S2
      section => this%get_section(s_idx)
      integrals(s_idx) = sum(real(section))*dA   ! S1
    end do

    open(15, position="append", file=filename)
    write(15, '(5(F10.7,1x))') integrals(F1:F3), integrals(S1), integrals(S2)
    close(15)
  end subroutine

  subroutine write_cline_integrals(this, filename)
    class(cfield), intent(in) :: this
    character(len=*), intent(in) :: filename
    complex :: left, right
    complex, dimension(F1:F3) :: cline_integrals
    integer :: s_idx

    do s_idx=F1,F3
      left  = sum(this%m_values(:,L_idx(s_idx)))*dZ
      right = sum(this%m_values(:,R_idx(s_idx)))*dZ
      cline_integrals(s_idx) = right - left
    end do

    open(15, position="append", file=filename)
    write(15, '(6(F10.7,1x))') real(cline_integrals(:)), imag(cline_integrals(:))
    close(15)
  end subroutine

  function rline_integrals(this)
    class(cfield), intent(in)   :: this
    real, dimension(F1:F3)      :: rline_integrals
    real                        :: left, right
    integer                     :: s_idx

    do s_idx=F1,F3
      left =sum(real(this%m_values(:,L_idx(s_idx))))*dZ  ! L_F1
      right=sum(real(this%m_values(:,R_idx(s_idx))))*dZ  ! R_F1
      rline_integrals(s_idx) = right - left
    end do
  end function

 subroutine write_torque_file(this, filename, torque_field)
    class(cfield), intent(in)    :: this
    class(rfield), intent(in),target    :: torque_field
    character(len=*), intent(in) :: filename

    real, dimension(F1:F3)  :: torques
    integer                 :: s_idx

    do s_idx = F1,F3
      torques(s_idx) = sum(torque_field%get_section(s_idx)) * dA
    end do

    open(15, position="append", file=filename)
    write(15, '(6(F12.7,1x))') torques(:), this%rline_integrals()
    close(15)
  end subroutine

  subroutine all_reduce(this, operation, comm)
    !use mpi
    include "mpif.h"
    class(cfield), intent(inout) :: this
    integer(kind=4), intent(in) :: operation
    integer(kind=4), intent(in), optional :: comm
    complex, dimension(0:npz-1,L_idx(S1):R_idx(S2)) :: values
    integer(kind=4) :: comm_, ierr, y_idx
    integer :: cutoff = 1048576  ! 2^20
    integer :: cols_at_once

    cols_at_once = cutoff / npz

    comm_ = mpi_comm_world
    if(present(comm)) then
        comm_ = comm
    endif

    ! Perform Reduce on ~2^20 elements at a time (max number of columns <= 2^20)
    y_idx = L_idx(S1)
    do while (y_idx + cols_at_once < R_idx(S2))
        call MPI_Allreduce(this%m_values(:,y_idx:y_idx+cols_at_once - 1),   &
                            values(:,y_idx:y_idx+cols_at_once - 1),         &
                            npz*cols_at_once, MPI_COMPLEX, operation, comm_, ierr)
        y_idx = y_idx + cols_at_once
    enddo

    ! get final set of columns
    if(y_idx <= R_idx(S2)) then
        call MPI_Allreduce(this%m_values(:,y_idx:R_idx(S2)),                &
                            values(:,y_idx:R_idx(S2)),                      &
                            npz*(R_idx(S2)-y_idx + 1),                      &
                            MPI_COMPLEX, operation, comm_, ierr)
    endif
    this%m_values = values
  end subroutine

  subroutine broadcast(this, root, comm)
    !use mpi
    include "mpif.h"
    class(cfield), intent(inout) :: this
    integer(kind=4), intent(in)  :: root
    integer(kind=4), intent(in), optional :: comm
    integer(kind=4) :: comm_, ierr, y_idx
    integer, parameter :: cutoff = 1048576  ! 2^20
    integer :: cols_at_once

    cols_at_once = cutoff / npz

    comm_ = mpi_comm_world
    if(present(comm)) then
        comm_ = comm
    endif

    ! Perform Reduce on ~2^20 elements at a time (max number of columns <= 2^20)
    y_idx = L_idx(S1)
    do while (y_idx + cols_at_once < R_idx(S2))
        call MPI_Bcast(this%m_values(:,y_idx:y_idx+cols_at_once - 1),   &
                       npz*cols_at_once, MPI_COMPLEX, root, comm_, ierr)
        y_idx = y_idx + cols_at_once
    enddo

    ! get final set of columns
    if(y_idx <= R_idx(S2)) then
        call MPI_Bcast(this%m_values(:,y_idx:R_idx(S2)),                &
                       npz*(R_idx(S2)-y_idx + 1),                       &
                       MPI_COMPLEX, root, comm_, ierr)
    endif
  end subroutine
end module cfield_class
