module constants
    implicit none

    ! parameter variables that are used by the program
    ! and external subroutines
    real   , parameter  :: pi = 4*atan(1.0)     ! use Fortran intrinsic to define PI
    integer, parameter  :: false=0
    integer, parameter  :: true=1
    complex, parameter  :: cZero = (0.0, 0.0) ! complex zero
    complex, parameter  :: im = (0.0,1.0) ! convenient way to write 'i'
    real   , parameter  :: degrees=pi/180.0

    character(99),parameter :: g_txt="g.txt" 
    character(99),parameter :: g_coord_txt="g_coord.txt"
    character(99),parameter :: jcharge_txt="jcharge.txt" 
    character(99),parameter :: density_txt="density.txt"
    character(99),parameter :: avgjcharge_txt="avgjcharge.txt" 
    character(99),parameter :: free0_txt="free0.txt"
    character(99),parameter :: avgf3_txt="avgf3.txt"
    character(99),parameter :: avgf0_txt="avgf0.txt"
    character(99),parameter :: avgf1_txt="avgf1.txt"
    character(99),parameter :: avgf2_txt="avgf2.txt"
    character(99),parameter :: maxf3_txt="maxf3.txt"
    character(99),parameter :: maxf0_txt="maxf0.txt"
    character(99),parameter :: maxf1_txt="maxf1.txt"
    character(99),parameter :: maxf2_txt="maxf2.txt"
    character(99),parameter :: avgf0rot_txt="avgf0_rot.txt"
    character(99),parameter :: avgf1rot_txt="avgf1_rot.txt"
    character(99),parameter :: avgf2rot_txt="avgf2_rot.txt"
    character(99),parameter :: maxf0rot_txt="maxf0_rot.txt"
    character(99),parameter :: maxf1rot_txt="maxf1_rot.txt"
    character(99),parameter :: maxf2rot_txt="maxf2_rot.txt"
    character(99),parameter :: mm_txt="mm.txt"
    character(99),parameter :: mm_y="mm_y.txt"
    character(99),parameter :: mm_z="mm_z.txt"
    character(99),parameter :: avgmx_txt="avgmx.txt"
    character(99),parameter :: avgmy_txt="avgmy.txt"
    character(99),parameter :: avgmz_txt="avgmz.txt"   
    character(99),parameter :: maxmx_txt="maxmx.txt"
    character(99),parameter :: maxmy_txt="maxmy.txt"
    character(99),parameter :: maxmz_txt="maxmz.txt"   
    character(99),parameter :: pa_txt="pa.txt"
    character(99),parameter :: pa1D_y="pa1D_y.txt"
    character(99),parameter :: pa1D_z="pa1D_z.txt" 
    character(99),parameter :: parot_txt="pa_rot.txt"
    character(99),parameter :: pa1Drot_y="pa1D_y_rot.txt"
    character(99),parameter :: pa1Drot_z="pa1D_z_rot.txt" 
    character(99),parameter :: dens_y="density_y.txt"
    character(99),parameter :: dens_z="density_z.txt"
    character(99),parameter :: jc_y="jcharge_y.txt"
    character(99),parameter :: jc_z="jcharge_z.txt"
    character(99),parameter :: divjc="divjc.txt"
    character(99),parameter :: iter_txt="iter.txt"
    character(99),parameter :: syxavg_txt="Syxavg.txt"
    character(99),parameter :: syyavg_txt="Syyavg.txt"
    character(99),parameter :: syzavg_txt="Syzavg.txt"
    character(99),parameter :: szxavg_txt="Szxavg.txt"
    character(99),parameter :: szyavg_txt="Szyavg.txt"
    character(99),parameter :: szzavg_txt="Szzavg.txt"
    character(99),parameter :: jspin1D_y="jspin_y.txt" 
    character(99),parameter :: jspin1D_z="jspin_z.txt"
    character(99),parameter :: tx_tot_txt="tx_tot.txt"
    character(99),parameter :: ty_tot_txt="ty_tot.txt"
    character(99),parameter :: tz_tot_txt="tz_tot.txt"
    character(99),parameter :: djspin_txt="djspin.txt"
    character(99),parameter :: jspin_txt="jspin.txt"
    character(99),parameter :: torque_txt="torque.txt"
    character(99),parameter :: tx_line_txt="tx_line.txt"
    character(99),parameter :: dtavg_txt="dtavg_tot.txt"
    character(99),parameter :: w_txt="W.txt"
    character(99),parameter :: w_min_txt="W_min.txt"

end module constants
