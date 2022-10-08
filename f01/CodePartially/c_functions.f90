module c_functions
    use iso_c_binding
    implicit none

    interface
        subroutine my_mkdir(string) bind(C, name='my_mkdir')
            use iso_c_binding, only: c_char
            character(kind=c_char) :: string(*)
        end subroutine my_mkdir
    end interface

contains
    subroutine mkdir(string)
        use iso_c_binding, only:c_char, c_null_char
        implicit none

        character(len=*), intent(in) :: string
        character(kind=c_char,len=1024) :: c_string
        integer :: n

        n = len(string)

        c_string = ''

        ! need to append a null char to the end of the string
        write(c_string, '(a,a)') string, c_null_char

        call my_mkdir(c_string)
    end subroutine
end module

