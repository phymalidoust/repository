program real_size
    implicit none
    real :: a = 1.0

    open(unit=10, file='tmp.bin', form='unformatted', access='stream')
    write(10) a
    write(10) a
    close(10)
    
end program real_size
