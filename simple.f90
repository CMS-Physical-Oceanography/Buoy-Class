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

    dt(:,:) = 7
    print*, 'Crunching Numbers...'
    do t = 1,size(logscale) ! loop through time
        W_copy = W
        wndstress = ws*logscale(t)

        W(0,:) = W_copy(0,:) + dt(0,:)*PG(0,:) & ! Pressure Gradient
        - j*dt(0,:)*coriolis(lat)*W_copy(0,:)  & ! Coriolis Term
        + ((dt(0,:)*Av(0,:))/(dz(0,:)**2))    &
        *(2*W_copy(1,:) - 2*W_copy(0,:) + ((2*dz(0,:)*wndstress)/(swdens*Av(0,:)))) &
        +(dt(0,:)/(dz(0,:)**2))*(W_copy(0,:)-W_copy(1,:))*(Av(0,:)-Av(1,:)) 

        W(2:size(W(:,0))-1,:) = W_copy(2:size(W(:,0))-1,:) + dt(2:size(W(:,0))-1,:)*PG(2:size(W(:,0))-1,:) & 
        - j*dt(2:size(W(:,0))-1,:)*coriolis(lat)*W_copy(2:size(W(:,0))-1,:) & ! Coriolis Term
        + ((dt(2:size(W(:,0))-1,:)*Av(2:size(W(:,0))-1,:))/(dz(2:size(W(:,0))-1,:)**2))  &
        *(W_copy(1:size(W(:,0))-2,:) - 2*W_copy(2:size(W(:,0))-1,:) + W_copy(3:,:)) &
        +(dt(2:size(W(:,0))-1,:)/(dz(2:size(W(:,0))-1,:)**2)) &
        *(W_copy(2:size(W(:,0))-1,:)-W_copy(3:,:))*(Av(2:size(W(:,0))-1,:)-Av(3:,:)) 
    end do
    print*,'N time-steps=',size(logscale)
    print*,'Saving...'
    
    call saveNdreal(real(W),'CSflow.bin',Nz,lendata)
    call saveNdreal(aimag(W),'ASflow.bin',Nz,lendata)

    deallocate(logscale)
    deallocate(dz)
    deallocate(dt)
    deallocate(Av)
    deallocate(pg)
    deallocate(ws)
end program simple