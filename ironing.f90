!....................................................
!..         Subroutine ironing                    ...
!....................................................
subroutine ironing(f,nx,ny,nz)
implicit none

integer (kind=4),intent(in)   :: nx,ny,nz
real    (kind=8)              :: f(nx,ny,nz)
real    (kind=8), allocatable :: aux(:,:,:)

integer (kind=4)              :: ix,iy,iz
real    (kind=8), parameter   :: a0=8.333333333333333d-2   ! = 1/12
real    (kind=8), parameter   :: a1=0.5d0

allocate(aux(nx,ny,nz))

aux = f

do iz=2,nz-1
  do iy=2,ny-1
    do ix=2,nx-1
      f(ix,iy,iz) =  a0*( aux(ix-1,iy  ,iz  )+aux(ix+1,iy  ,iz)       &
                         +aux(ix  ,iy-1,iz  )+aux(ix  ,iy+1,iz)       &
                         +aux(ix  ,iy  ,iz-1)+aux(ix  ,iy  ,iz+1) )   &
                  +  a1*aux(ix,iy,iz)
    end do
  end do
end do

deallocate(aux)

return
end
