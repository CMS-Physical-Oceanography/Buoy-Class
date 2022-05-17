program simple 
    use setup;use Ocean;use savearry
    implicit none 
    ! ==========================================
    !               DATA ARRAYS
    ! Depth Unifom Geostrophic Pressure Gradient
    complex,allocatable :: pg(:,:),W(:,:),W_copy(:,:)
    ! Windstress
    complex,allocatable :: ws(:),wndstress(:)
    ! ==========================================
    !               PARAMETERS
    real,allocatable :: Av(:,:),dz(:,:),dt(:,:),logscale(:)
    ! ===================================x=======
    !                 INPUTS
    integer :: lendata,Nz,t
    print*,'loading Wind Stress...'
    call WindStress(ws,lendata) ! (csws,asws)
    print*,'Loading Eddy Viscosity...'
    call EddyViscosity(Av,lendata,Nz) ! Nz x lendata
    print*,'Loading Pressure Gradient...'
    call PressureGradient(pg,Nz) ! f(-vda,uda)
    print*,'Creating Grid...'
    call Discritization(dz,dt,logscale,Nz,lendata)
    ! ==========================================
    !                  MAIN
    ! ALLOCATE MEMORY
    allocate(W(Nz,lendata))
    allocate(W_copy(Nz,lendata))
    allocate(wndstress(size(ws)))

    W(:,:) = complex(0,0)

    dt(:,:) = 12
    print*, 'Crunching Numbers...'
    do t = 1,size(logscale) ! loop through time
        W_copy = W
        wndstress = ws*logscale(t)

        W(1,:) = W_copy(1,:) + dt(1,:)*PG(1,:) & ! Pressure Gradient
        - j*dt(1,:)*coriolis(lat)*W_copy(1,:)  & ! Coriolis Term
        + ((dt(1,:)*Av(1,:))/(dz(1,:)**2))     &
        *(2*W_copy(2,:) - 2*W_copy(1,:) + ((2*dz(1,:)*wndstress)/(swdens*Av(1,:)))) &
        +(dt(1,:)/(dz(1,:)**2))*(W_copy(1,:)-W_copy(2,:))*(Av(1,:)-Av(2,:)) 
        !size(W(:,1))
        W(2:Nz-1,:) = W_copy(2:Nz-1,:) + dt(2:Nz-1,:)*PG(2:Nz-1,:) & 
        - j*dt(2:Nz-1,:)*coriolis(lat)*W_copy(2:Nz-1,:) & ! Coriolis Term
        + ((dt(2:Nz-1,:)*Av(2:Nz-1,:))/(dz(2:Nz-1,:)**2))  &
        *(W_copy(1:Nz-2,:) - 2*W_copy(2:Nz-1,:) + W_copy(3:,:)) &
        +(dt(2:Nz-1,:)/(dz(2:Nz-1,:)**2)) &
        *(W_copy(2:Nz-1,:)-W_copy(3:,:))*(Av(2:Nz-1,:)-Av(3:,:)) 
    end do
    print*,coriolis(lat)
    print*,'N time-steps=',size(logscale)
    print*,'Saving...'
    
    call saveNdreal(real(W),'CSflow.bin',Nz)
    call saveNdreal(aimag(W),'ASflow.bin',Nz)

    deallocate(logscale)
    deallocate(dz)
    deallocate(dt)
    deallocate(Av)
    deallocate(pg)
    deallocate(ws)
end program simple