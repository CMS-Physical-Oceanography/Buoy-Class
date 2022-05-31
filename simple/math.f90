module math

    implicit none

    real,parameter :: pi = 3.141592
    real,parameter :: e = 2.718281

    public :: pi,e,deg2rad

    contains 

    real function deg2rad(deg) result(rad)
        implicit none 

        real, intent(in) :: deg ! angle in degrees to convert to radians
        rad = (pi*deg)/180
    end function deg2rad
    
end module math 