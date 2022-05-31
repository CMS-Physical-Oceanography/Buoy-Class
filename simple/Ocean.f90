module Ocean
    use math
    implicit none

    real, parameter :: swdens = 1025.0
    real,parameter :: OmegaE = (2*pi)/86400 

    public :: swdens,OmegaE,coriolis
    private
    contains 
        !==============================================
        !FUNCTION: coriolis(latt)result(cor)
        real function coriolis(latt) result(cor)
            implicit none 

            real, intent(in) :: latt ! lattitude in deg

            cor = 2*OmegaE*sin(deg2rad(latt))
        end function coriolis
end module Ocean