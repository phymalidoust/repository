module utilities
  use field_constants
  use constants

  implicit none

  double precision  :: start_time

  contains
    subroutine time_stamp(fid)
      include "mpif.h"
        integer, intent(in), optional :: fid
        double precision :: current_time, delta_time
        integer          :: hours, minutes, seconds

        current_time = MPI_Wtime()
        delta_time = current_time - start_time
        hours = floor(delta_time / 3600)
        delta_time = delta_time - (hours * 3600)
        minutes = floor(delta_time / 60)
        seconds = nint(delta_time - minutes * 60)

        if(present(fid)) then
            write(fid,'(I0.4,A,I0.2,A,I0.2,A)', advance="no") &
                hours, ':', minutes, ':', seconds, ' '
        else
            write(*,'(I0.4,A,I0.2,A,I0.2,A)', advance="no") &
                hours, ':', minutes, ':', seconds, ' '
        endif
    end subroutine

    subroutine status_message(txt, num)
      !use mpi
      include "mpif.h"
        character(len=*), intent(in)    :: txt
        integer, intent(in), optional   :: num
        integer(kind=4) :: rank, ierr

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        if (rank == 0) then
            call time_stamp()
            if(present(num)) then
                write(*,*) txt, num
            else
                write(*,*) txt
            endif
        endif
    end subroutine

    subroutine file_status_message(txt, num)
        use user_input

        character(len=*), intent(in)    :: txt
        integer, intent(in), optional   :: num

        call time_stamp(local_log)
        if(present(num)) then
            write(local_log,*) txt, num
        else
            write(local_log,*) txt
        endif
        flush(local_log)
    end subroutine

    character(99) function remove_blanks(in_str)
    
    character(*),intent(in) :: in_str
    character(99) :: out_str
    character :: ch
    integer :: j

      out_str = " "
      do j=1,len_trim(in_str)
         ch = in_str(j:j)
         if (ch /= " ") then
            out_str = trim(out_str) // ch
         endif
         remove_blanks = out_str
      end do
    end function remove_blanks

end module utilities
