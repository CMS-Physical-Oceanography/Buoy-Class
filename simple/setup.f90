module setup
    use  arrayin
    implicit none

    character(99),parameter :: udafile = 'uda.dat' ! u comp of depth averaged flow
    character(99),parameter :: vdafile = 'vda.dat' ! v comp of depth averaged flow
    character(99),parameter :: cswsfile = 'csws.dat' ! cross-shelf comp of wind stress
    character(99),parameter :: aswsfile = 'asws.dat' ! along-shelf comp of wind stress
    character(99),parameter :: Avfile = 'Av.dat' ! Eddy Viscosities
    character(99),parameter :: dzfile = 'dz.dat' ! vertical spacial steps
    character(99),parameter :: dtfile = 'dt.dat' ! time steps
    character(99),parameter :: scalefile = 'scale.dat' ! time steps

    real,parameter :: lat = 33
    complex,parameter :: j = (0,1) ! sqrt(-1)
    ! integer,parameter :: Nz = 39

    public :: PressureGradient,WindStress,EddyViscosity,Discritization,lat,j
    private

    contains 
    
    subroutine EddyViscosity(Av,lendata,Nz)
        use arrayin
        implicit none
        real,allocatable, intent(out) :: Av(:,:)
        integer,intent(out) :: Nz
        integer,intent(in) :: lendata

        Av = Av_in(Avfile,lendata) ! vertical Eddy Viscosity
        Nz = size(Av(:,1)) ! vertical gridpoints

    end subroutine EddyViscosity
! ===============================================================
    subroutine PressureGradient(pg,Nz) 
        ! This subroutine initializes a
        ! pressure gradient term assuming 
        ! depth averaged geostrophically balanced 
        ! flow. Input files are read in as 1D column
        ! vectors with real_in1d(fil) from the 
        ! array in module.
        ! ===============================================
        ! OUTPUTS DEPTH AVERAGED CORIOLIS TERM:
        !    pg = f(-vda+i*uda)
        ! Where f is the coriolis parameter and (uda,vda)
        ! are the depth averaged flow components. 
        ! coriolis(latt) is contained in Ocean.f90
        ! module.
        ! ===============================================
        ! OUTPUT ARRAY IS Nz x size(da):
        !    pg = [[pg(0,0) ... pg(size(uda))]
        !            :                  :       
        !          [pg(Nz,0) ...pg(size(uda))]]
        ! Where Nz is the number of vertical gridpoints
        ! which should be gotten with discritization(discfile)
        ! function if I wrote that (5/10/22).
        ! ===============================================
        use arrayin;use Ocean
        implicit none 

        real,allocatable:: uda(:),vda(:)
        complex,allocatable :: pg(:,:)
        integer :: i
        integer,intent(in) :: Nz
    
        uda = real_in1d(udafile) ! read uda.dat
        vda = real_in1d(vdafile) ! read vda.dat

        allocate(pg(Nz,size(uda))) ! allocate memory for pg

        do i = 1,size(uda) ! wrap into complex vectors
            pg(:,i) = complex(-1*vda(i),uda(i))
        end do
        pg(:,:) = pg(:,:) * coriolis(lat) ! scale by coriolis parameter
        deallocate(uda)
        deallocate(vda)
    end subroutine PressureGradient
! ===============================================================
    subroutine WindStress(ws,lendata) 

        ! This subroutine initializes a 
        ! complex windstress vector of the 
        ! same form as pg in PressureGradient(pg)
        ! subroutine above.
        ! ===============================================
        ! OUTPUTS WIND STRESS:
        !    ws = (csws+i*asws)
        ! Where csws is the across-shelf component
        ! of the wind stress and asws is the 
        ! along-shelf comp.
        ! ===============================================
        use arrayin
        implicit none 

        real,allocatable:: csws(:),asws(:)
        complex,allocatable :: ws(:)
        integer :: i,lendata
    
        csws = real_in1d(cswsfile) ! read csws.dat
        asws = real_in1d(aswsfile) ! read asws.dat

        lendata = size(csws)
        allocate(ws(lendata)) ! allocate memory for ws

        do i = 1,lendata ! wrap into complex vector
            ws(i) = complex(csws(i),asws(i))
        end do

        ! dealocate memory
        deallocate(csws)
        deallocate(asws)
    end subroutine WindStress
! ===============================================================
    subroutine Discritization(dz,dt,logscale,Nz,lendata)
        use arrayin
        implicit none

        real,allocatable ::dz(:,:),dt(:,:)
        real, allocatable :: logscale(:)
        integer, intent(in) :: Nz,lendata

        allocate(dz(Nz,lendata))
        allocate(dt(Nz,lendata))

        dz = dz_in(dzfile,Nz,lendata)
        dt = dt_in(dtfile,Nz,lendata)
        logscale = real_in1d(scalefile)

    end subroutine Discritization


end module setup