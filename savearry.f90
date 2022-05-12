module savearry
    implicit none

    contains
    subroutine saveNdreal(arr,fil,Nrows,Ncols)
        ! This function saves an Nd array to a
        ! a .txt file. Filename is inputted as fil
        ! array. 
        implicit none

        character(len=*), intent(in) :: fil
        real,dimension(:,:),intent(in) ::  arr(:,:)
        integer :: Nrows,Ncols,i

        open(unit=1,file=fil,action='write')

        do i = 1,Nrows ! loop through rows
            write(1,*) arr(i,:) ! save row
        end do 
    end subroutine saveNdreal
end module savearry