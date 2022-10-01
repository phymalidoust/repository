! this module sets up the various constants used for describing 

module field_constants
    implicit none

    ! index for each of the sections of the junction
    integer, parameter :: S1 = 1    ! Superconductor 1
    integer, parameter :: F1 = 2    ! Ferromagnet 1
    integer, parameter :: F2 = 3    ! Ferromagnet 2
    integer, parameter :: F3 = 4    ! Ferromagnet 3
    integer, parameter :: S2 = 5    ! Superconductor 2

    ! index for references Y or Z coordinates
    integer, parameter :: Z = 1
    integer, parameter :: Y = 2

end module field_constants
