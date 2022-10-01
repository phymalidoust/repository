module hdf5_utilities
    use hdf5

    implicit none

    integer(hid_t) :: h5_complex
contains
    subroutine h5_init()
        implicit none
        integer :: ierr

        call h5open_f(ierr)
    end subroutine

    subroutine write_slab(dset_id, values, layer)
        implicit none

        integer(hid_t), intent(in)          :: dset_id
        real, dimension(:,:), intent(in)    :: values
        integer, intent(in)                 :: layer

        integer(hid_t)                  :: memspace, filespace
        integer(hsize_t), dimension(3)  :: start, count
        integer, dimension(2)           :: rank
        integer                         :: ierr

        call H5Dget_space_f(dset_id, filespace, ierr)

        rank = shape(values)

        start = [1, 1, layer]
        count = [rank(1), rank(2), 1]

        call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                                    start, count, ierr)

        call H5Screate_simple_f(3, count, memspace, ierr)

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, values, count, ierr, &
                memspace, filespace)

        call H5Sclose_f(memspace, ierr)
        call H5Sclose_f(filespace, ierr)
        
    end subroutine

    subroutine write_real_array_to_hdf5(loc_id, dset_name, values)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        character(len=*), intent(in)    :: dset_name
        real, dimension(:), intent(in)  :: values

        integer(hsize_t), dimension(1)  :: dims
        integer(hid_t)                  :: filespace, dset_id
        integer                         :: ierr
        logical                         :: already_exists

        dims(1) = size(values)

        call H5Lexists_f(loc_id, dset_name, already_exists, ierr)

        if(already_exists) then
            call H5Dopen_f(loc_id, dset_name, dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
        else
            call H5Screate_simple_f(1, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, filespace, &
                             dset_id, ierr)
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, values, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    subroutine write_real_mat_to_hdf5(loc_id, dset_name, values)
        implicit none

        integer(hid_t), intent(in)          :: loc_id
        character(len=*), intent(in)        :: dset_name
        real, dimension(:,:), intent(in)    :: values

        integer(hsize_t), dimension(2)      :: dims
        integer(hid_t)                      :: filespace, dset_id
        integer                             :: ierr
        logical                             :: already_exists

        dims = shape(values)

        call H5Lexists_f(loc_id, dset_name, already_exists, ierr)

        if(already_exists) then
            call H5Dopen_f(loc_id, dset_name, dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
        else
            call H5Screate_simple_f(2, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, filespace, &
                             dset_id, ierr)
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, values, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    subroutine write_complex_array_to_hdf5(loc_id, dset_name, values)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        character(len=*), intent(in)    :: dset_name
        complex, dimension(:), intent(in)  :: values

        integer(hsize_t), dimension(2)  :: dims
        integer(hid_t)                  :: filespace, dset_id
        integer                         :: ierr
        logical                         :: already_exists
        real, dimension(:,:), allocatable :: scratch

        dims(1) = size(values)
        dims(2) = 2

        allocate(scratch(dims(1), 2))

        scratch(:,1) = real(values)
        scratch(:,2) = imag(values)

        call H5Lexists_f(loc_id, dset_name, already_exists, ierr)

        if(already_exists) then
            call H5Dopen_f(loc_id, dset_name, dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
        else
            call H5Screate_simple_f(2, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, dset_name, H5T_NATIVE_REAL, filespace, &
                             dset_id, ierr)
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    subroutine write_int_array_to_hdf5(loc_id, dset_name, values)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        character(len=*), intent(in)    :: dset_name
        integer, dimension(:), intent(in)  :: values

        integer(hsize_t), dimension(1)  :: dims
        integer(hid_t)                  :: filespace, dset_id
        integer                         :: ierr
        logical                         :: already_exists

        dims(1) = size(values)

        call H5Lexists_f(loc_id, dset_name, already_exists, ierr)

        if(already_exists) then
            call H5Dopen_f(loc_id, dset_name, dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
        else
            call H5Screate_simple_f(1, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, &
                             filespace, dset_id, ierr)
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    subroutine append_iter_to_hdf5(loc_id, it, error)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        integer, intent(in)             :: it
        real, intent(in)                :: error

        integer                         :: ierr, i
        logical                         :: already_exists
        integer(hid_t)                  :: dset_id, filespace
        integer(hsize_t), dimension(2)  :: dims

        real, dimension(2,100)          :: scratch

        scratch = -1
        dims = [2, 100]

        call H5Lexists_f(loc_id, 'iter', already_exists, ierr)

        if(already_exists) then
            ! find first "empty" row, and replace it
            call H5Dopen_f(loc_id, 'iter', dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)

            do i=1,100
                if(scratch(1,i) == -1) exit
            end do

            scratch(1,i) = real(it)
            scratch(2,i) = error
        else
            call H5Screate_simple_f(2, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, 'iter', H5T_NATIVE_REAL, &
                             filespace, dset_id, ierr)

            scratch(1,1) = it
            scratch(2,1) = error
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    subroutine append_avgjcharge_to_hdf5(loc_id, it, charge_arr, phase)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        integer, intent(in)             :: it
        real,dimension(3), intent(in)   :: charge_arr
        real, intent(in)                :: phase

        integer                         :: ierr, i
        logical                         :: already_exists
        integer(hid_t)                  :: dset_id, filespace
        integer(hsize_t), dimension(2)  :: dims

        real, dimension(5,100)          :: scratch

        scratch = -1
        dims = [5, 100]

        call H5Lexists_f(loc_id, 'avg_j_charge', already_exists, ierr)

        if(already_exists) then
            ! find first "empty" row, and replace it
            call H5Dopen_f(loc_id, 'avg_j_charge', dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)

            do i=1,100
                if(scratch(1,i) == -1) exit
            end do

            scratch(1,i) = real(it)
            scratch(2,i) = charge_arr(1)
            scratch(3,i) = charge_arr(2)
            scratch(4,i) = charge_arr(3)
            scratch(5,i) = phase
            
        else
            call H5Screate_simple_f(2, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, 'avg_j_charge', H5T_NATIVE_REAL, &
                             filespace, dset_id, ierr)

            scratch(1,1) = real(it)
            scratch(2,1) = charge_arr(1)
            scratch(3,1) = charge_arr(2)
            scratch(4,1) = charge_arr(3)
            scratch(5,1) = phase
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine


    subroutine append_free_to_hdf5(loc_id, it, free_arr)
        implicit none

        integer(hid_t), intent(in)      :: loc_id
        integer, intent(in)             :: it
        real, dimension(3), intent(in)  :: free_arr

        integer                         :: ierr, i
        logical                         :: already_exists
        integer(hid_t)                  :: dset_id, filespace
        integer(hsize_t), dimension(2)  :: dims

        real, dimension(4,100)          :: scratch

        scratch = -1
        dims = [4, 100]

        call H5Lexists_f(loc_id, 'free0', already_exists, ierr)

        if(already_exists) then
            ! find first "empty" row, and replace it
            call H5Dopen_f(loc_id, 'free0', dset_id, ierr)
            call H5Dget_space_f(dset_id, filespace, ierr)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)

            do i=1,100
                if(scratch(1,i) == -1) exit
            end do

            scratch(1,i) = real(it)
            scratch(2,i) = free_arr(1)
            scratch(3,i) = free_arr(2)
            scratch(4,i) = free_arr(3)
        else
            call H5Screate_simple_f(2, dims, filespace, ierr, dims)
            call H5Dcreate_f(loc_id, 'free0', H5T_NATIVE_REAL, &
                             filespace, dset_id, ierr)

            scratch(1,1) = real(it)
            scratch(2,1) = free_arr(1)
            scratch(3,1) = free_arr(2)
            scratch(4,1) = free_arr(3)
        endif

        call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, scratch, dims, ierr)
        call H5Dclose_f(dset_id, ierr)
        call H5Sclose_f(filespace, ierr)
    end subroutine

    integer(hid_t) function create_hdf5(filename, is_root)
        implicit none

        character(len=*), intent(in)    :: filename
        logical, intent(in)             :: is_root

        integer                         :: h5_err

        if(is_root) then
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, create_hdf5, h5_err)
        else
            create_hdf5 = -1
        endif
    end function

    integer(hid_t) function open_hdf5(filename, is_root)
        implicit none

        character(len=*), intent(in)    :: filename
        logical, intent(in)             :: is_root

        integer                         :: h5_err

        if(is_root) then
            call h5fopen_f(filename, H5F_ACC_RDWR_F, open_hdf5, h5_err)
        else
            open_hdf5 = -1
        endif
    end function

    subroutine close_hdf5(file_id)
        integer(hid_t), intent(inout) :: file_id

        integer :: h5_err
        if(file_id /= -1) then
            call h5fclose_f(file_id, h5_err)
            file_id = -1
        endif
    end subroutine
end module
