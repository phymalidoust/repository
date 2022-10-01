

module rfield_class
  use field_constants
  use user_input

  implicit none
  private                           ! hide everything by default

  public :: rfield
  type rfield
    private
    real, dimension(:,:), allocatable :: m_values
  contains
    procedure   :: set_to
    procedure   :: set_at
    generic     :: set => set_to, set_at

    procedure   :: set_section_to_real      ! set all of a section to a given real value
    procedure   :: set_section_to_array     ! set all of a section to the values in an array
    generic     :: set_section => set_section_to_real, set_section_to_array

    procedure   :: get
    procedure   :: get_section              ! returns an allocated array with the values in that section

    procedure   :: inc_at
    procedure   :: inc_by_real
    procedure   :: inc_by_rfield
    generic     :: operator(+) => inc_by_real,inc_by_rfield

    procedure   :: dec_by_rfield
    generic     :: operator(-) => dec_by_rfield

    procedure   :: mult_by_real
    procedure   :: mult_by_rfield
    generic     :: operator(*) => mult_by_real,mult_by_rfield

    procedure   :: maximum
    procedure   :: mag
    procedure   :: abs_diff
    procedure   :: maximum_w_kfdx

    procedure   :: to_array

    ! I/O subroutines
    procedure   :: write_to_file
    procedure   :: write_to_hdf5
    procedure   :: write_to_file_with_pos
    procedure   :: read_from_file

    procedure :: all_reduce
  end type rfield

  interface rfield
    procedure array_constructor
    procedure value_constructor
    procedure copy_constructor
  end interface

contains
! constructor sets the values of all elements to the values of the array
  function array_constructor(values) result(new_field)
    real, dimension(0:npz-1,L_idx(S1):R_idx(S2)), intent(in) :: values
    type(rfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = values
  end function

! constructor sets the values of all elements to the intial value
  function value_constructor(value) result(new_field)
    real, intent(in) :: value
    type(rfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = value
  end function

! constructor creates a copy of the provided constructor
  function copy_constructor(old_field) result(new_field)
    type(rfield), intent(in) :: old_field
    type(rfield) :: new_field

    allocate(new_field%m_values(0:npz-1,L_idx(S1):R_idx(S2)))
    new_field%m_values = old_field%m_values
  end function

! returns the value at a specfic element
  pure function get(this, z_idx, y_idx) result(value)
    class(rfield), intent(in) :: this
    integer, intent(in) :: z_idx, y_idx
    real :: value
    value = this%m_values(z_idx, y_idx)
  end function

! sets the value of a specific element
  subroutine set_at(this, z_idx, y_idx, value)
    class(rfield), intent(inout) :: this
    integer, intent(in) :: z_idx, y_idx
    real, intent(in) :: value

    this%m_values(z_idx, y_idx) = value
  end subroutine

! sets all of the elemnts to the given value
  subroutine set_to(this, value)
    class(rfield), intent(inout) :: this
    real, intent(in) :: value

    this%m_values = value
  end subroutine

! inc the element at z_idx,y_idx by the provided amount
  subroutine inc_at(this, z_idx, y_idx, value)
    class(rfield), intent(inout) :: this
    real, intent(in) :: value
    integer, intent(in) :: z_idx, y_idx

    this%m_values(z_idx, y_idx) = this%m_values(z_idx, y_idx) + value
  end subroutine

! inc all elements by the given scalar value, return answer
  function inc_by_real(this, value) result(ans)
    class(rfield), intent(in) :: this
    type(rfield) :: ans
    real, intent(in) :: value

    ans = rfield(0.0)
    ans%m_values = this%m_values + value
  end function

! inc each element by the corresponding element in the addend
  function inc_by_rfield(this, addend) result(ans)
    class(rfield), intent(in) :: this
    class(rfield), intent(in) :: addend
    type(rfield) :: ans

    ans = rfield(0.0)
    ans%m_values = this%m_values + addend%m_values
  end function

! decrement each element by the corresponding element in the argument
  function dec_by_rfield(this, subtrahend) result(ans)
    class(rfield), intent(in) :: this
    class(rfield), intent(in) :: subtrahend
    type(rfield) :: ans

    ans = rfield(0.0)
    ans%m_values = this%m_values - subtrahend%m_values
  end function

! multiply all elements by the given scalar value, return answer
  function mult_by_real(this, value) result(ans)
    class(rfield), intent(in) :: this
    real, intent(in) :: value
    type(rfield) :: ans

    ans = rfield(0.0)
    ans%m_values = this%m_values * value
  end function

! multiply all elements by the given scalar value, return answer
  function mult_by_rfield(this, factor) result(ans)
    class(rfield), intent(in) :: this
    class(rfield), intent(in) :: factor
    type(rfield) :: ans

    ans = rfield(0.0)
    ans%m_values = this%m_values * factor%m_values
  end function

! write out the rfield
  subroutine write_to_file(this, filename)
    class(rfield), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer :: z_idx, y_idx

    open(88, file=filename, action="write")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        write(88,'(f10.7)') this%m_values(z_idx, y_idx)
      end do
    end do
    close(88)
  end subroutine

  subroutine write_to_hdf5(this, loc_id, dset_name)
    use hdf5

    class(rfield), intent(in)       :: this
    character(len=*), intent(in)    :: dset_name
    integer(hid_t), intent(in)      :: loc_id
    integer(hsize_t), dimension(2)  :: dims
    integer                         :: ierr
    integer(hid_t)                  :: dset_id, filespace

    dims = [npz, npy]

    call H5Screate_simple_f(2, dims, filespace, ierr, dims)
    call H5Dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, filespace, &
                     dset_id, ierr)

    call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, this%m_values, dims, ierr)

    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(filespace, ierr)
  end subroutine

! write out the rfield, include the position
  subroutine write_to_file_with_pos(this, filename)
    class(rfield), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer :: z_idx, y_idx

    open(88, file=filename, action="write")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        write(88,'(2(f10.51x), f10.7)') z_idx*dY, y_idx*dZ, this%m_values(z_idx, y_idx)
      end do
    end do
    close(88)
  end subroutine

! read in the rfield from a file
  subroutine read_from_file(this, filename)
    class(rfield), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: z_idx, y_idx

    open(88, file=filename, action="READ")
    do y_idx = 0,npy-1
      do z_idx = 0,npz-1
        read(88,'(f10.7)') this%m_values(z_idx,y_idx)
      end do
    end do
    close(88)
  end subroutine

  function get_section(this, idx)
    class(rfield), intent(in), target :: this
    integer, intent(in) :: idx
    real, dimension(:,:), pointer :: get_section

    get_section => this%m_values(:,L_idx(idx):R_idx(idx))
  end function

  subroutine set_section_to_real(this, section_idx, value)
    class(rfield), intent(inout)    :: this
    integer, intent(in)             :: section_idx
    real, intent(in)             :: value

    this%m_values(:,L_idx(section_idx):R_idx(section_idx)) = value
  end subroutine

  subroutine set_section_to_array(this, section_idx, arr)
    class(rfield), intent(inout)     :: this
    integer, intent(in)              :: section_idx
    real, dimension(:,:), intent(in) :: arr

    this%m_values(:,L_idx(section_idx):R_idx(section_idx)) = arr
  end subroutine

  function mag(this)    result(ans)
    class(rfield), intent(in)           :: this
    type(rfield)                        :: ans

    ans = rfield(0.0)
    ans%m_values = abs(this%m_values)
  end function

  function abs_diff(this, that)
    class(rfield), intent(in)           :: this
    class(rfield), intent(in)           :: that
    type(rfield)                        :: abs_diff

    abs_diff = this%mag() - that%mag()
    abs_diff = abs_diff%mag()
  end function

  pure function to_array(this)
    class(rfield), intent(in)           :: this
    real, dimension(:,:),allocatable    :: to_array

    allocate(to_array(0:npz-1, L_idx(S1):R_idx(S2)))
    to_array = this%m_values
  end function

  pure real function maximum(this)
    class(rfield), intent(in)           :: this

    maximum = maxval(this%m_values)
  end function

  pure real function maximum_w_kfdx(this)
    class(rfield), intent(in)           :: this
    integer :: left,right

    left = npsx+1
    right = R_idx(S2) - npsx-1

    maximum_w_kfdx = maxval(this%m_values(0:npz-1,left:right))
  end function

  subroutine all_reduce(this, operation, comm)
    !use mpi
    include "mpif.h"
    class(rfield), intent(inout) :: this
    integer(kind=4), intent(in) :: operation
    integer(kind=4), intent(in), optional :: comm
    real, dimension(0:npz-1,L_idx(S1):R_idx(S2)) :: values
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
                            npz*cols_at_once, MPI_REAL, operation, comm_, ierr)
        y_idx = y_idx + cols_at_once
    enddo

    ! get final set of columns
    if(y_idx <= R_idx(S2)) then
        call MPI_Allreduce(this%m_values(:,y_idx:R_idx(S2)),                &
                            values(:,y_idx:R_idx(S2)),                      &
                            npz*(R_idx(S2)-y_idx + 1),                      &
                            MPI_REAL, operation, comm_, ierr)
    endif
    this%m_values = values
  end subroutine
end module rfield_class
