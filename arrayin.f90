module arrayin
    implicit none
    contains
    ! ===============================================================
    function real_in(fil) result(data_)
        ! This function inputs file name fil
        ! containing an array and outputs a real
        ! array data_. Dimensions must be specified 
        ! as real or integer as the sole values in 
        ! the first column of the input array. 
        ! Input arrays can be created in python with 
        ! pandas.DataFrame(data).to_csv('fil.dat',header=False,index=False)
        ! Header and index must be passed as False for 
        ! this function to work.
        ! INPUT FORMAT:
        ! INPUT = [[Nrows,Ncols],
        !          [real(0,0) ......... real(0,Ncols)],
        !          [real(Nrows,0)..real(Nrows,Ncols)]]
        ! THIS WILL THROW AND ERROR IF MEMORY FOR DATA_ IS
        ! DEALLOCATED IN FUNCTION.
        implicit none
        !                  VARIABLES
        character(99),intent(in) :: fil ! input file name
        real,allocatable :: data_(:,:) ! output
        real :: cols,rows ! input array dimensions
        integer :: row,col ! index
        !                   OPEN FILE  
        open(unit=1,file=fil)
        !              CREATE OUTPUT ARRAY
        read(1,*) cols ! N columns
        read(1,*) rows ! N rows
        allocate(data_(int(rows),int(cols)))! allocate memory

        !              WRITE TO OUTPUT ARRAY
        do row = 1,int(rows) ! loop through file 
            do col = 1,int(cols)
                read(1,*) data_(row,col)
            end do
        end do
    end function real_in
    ! =================================================================================
    function real_in1d(fil) result(data_)
        
        ! SAME AS ABOVE JUST 1D
        implicit none
        !                  VARIABLES
        character(99),intent(in) :: fil ! input file name
        real,allocatable :: data_(:) ! output
        real :: len,null ! input array dimensions
        integer :: idx ! index
        !                   OPEN FILE  
        open(unit=1,file=fil)
        !              CREATE OUTPUT ARRAY
        read(1,*) len ! input array length (cols)
        allocate(data_(int(len)))! allocate memory
        len = len
        !              WRITE TO OUTPUT ARRAY
        do idx = 1,int(len) ! loop through file 
            read(1,*) data_(idx) ! data_(row) = fil(row)
        end do
    end function real_in1d
    ! ===============================================================
    function complex_in1d(fil) result(data_)
        
        ! SAME AS ABOVE JUST 1D
        implicit none
        !                  VARIABLES
        character(99),intent(in) :: fil ! input file name
        complex,allocatable :: data_(:) ! output
        complex :: len ! input array dimensions
        integer :: idx ! index
        !                   OPEN FILE  
        open(unit=1,file=fil)
        !              CREATE OUTPUT ARRAY
        read(1,*) len ! input array length (cols)
        allocate(data_(int(len)))! allocate memory

        !              WRITE TO OUTPUT ARRAY
        do idx = 1,int(len) ! loop through file 
            read(1,*) data_(idx) ! data_(row) = fil(row)
        end do
    end function complex_in1d
    ! ===============================================================
    function Av_in(fil,lendata) result(Av)
        implicit none 

        character(99),intent(in) :: fil
        integer :: Nz,lendata,row,col
        real :: Nz_
        real,allocatable :: Av(:,:), tempAv(:)

        open(unit=1,file=fil) 

        read(1,*) Nz_
        Nz = int(Nz_)

        allocate(Av(Nz,lendata))
        allocate(tempAv(Nz))

        do  col = 1,lendata
            do row = 1,Nz
                read(1,*) tempAv(row)
            end do 
            Av(:,col) = tempAv(:)
        end do
        deallocate(tempAv)
    end function Av_in
! ===============================================================
    function dz_in(fil,Nz,lendata) result(dz)
        implicit none 

        character(99),intent(in) :: fil
        integer,intent(in) :: Nz,lendata
        integer :: row
        real,dimension(Nz,lendata) :: dz
        real ::dz_

        open(unit=1,file=fil) 

        do  row = 1,Nz
            read(1,*) dz_
            dz(row,:) = dz_
        end do
    end function dz_in
! ===============================================================   
    function dt_in(fil,Nz,lendata) result(dt)
        implicit none 
        character(99),intent(in) :: fil
        integer,intent(in) :: Nz,lendata
        integer :: col
        real,dimension(Nz,lendata) :: dt
        real ::dt_
        open(unit=1,file=fil) 
        do  col = 1,lendata
            read(1,*) dt_
            dt(:,col) = dt_
        end do
    end function dt_in
end module arrayin