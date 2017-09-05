Subroutine sorteo(denx,nx,ny,nz,x,y,z,rimp,N_par,xi,yi,zi,xf,yf,zf)
implicit none


integer (kind=4) :: N_par, nx, ny, nz, i
real    (kind=8) :: denx(nx,ny,nz), x(nx), y(ny), z(nz) ! Densidad del 4He del sampling
real    (kind=8) :: rimp(N_par,3), dum, xmin, ymin, zmin, dx, dy, dz, denmax
real    (kind=8) :: xi, xf, yi, yf, zi, zf, xyz(3), f
xmin = x(1)
ymin = y(1)
zmin = z(1)
dx = x(2) -x(1); dy = y(2) -y(1); dz = z(2) -z(1)
denmax = MaxVal(denx)

!Write(6,'("From sorteo: denmax...:",1p,E15.6)')denmax
!Write(6,'("From sorteo: xi,yi,zi...:",1p,3E15.6)')xi,yi,zi
!Write(6,'("From sorteo: xf,yf,zf...:",1p,3E15.6)')xf,yf,zf

do i=1,N_par

  202 call random_number(xyz)
      xyz(1)=xyz(1)*(xf-xi)+xi
      xyz(2)=xyz(2)*(yf-yi)+yi
      xyz(3)=xyz(3)*(zf-zi)+zi
      call random_number(dum)
      dum=dum*denmax
      call interpola3denx(nx,ny,nz,xmin,ymin,zmin,dx,dy,dz,xyz,denx,f)
      if(dum.ge.f) goto 202
      Rimp(i,:)=xyz(:)
end do

return

end


subroutine interpola3denx(nx,ny,nz,x0,y0,z0,Dx,Dy,Dz,xyz,denx,f)

implicit none

real (kind=8)    :: x0,y0,z0,Dx,Dy,Dz,f,  r,s,t, ax,ay,az
integer (kind=4) :: nx,ny,nz,ix,iy,iz
real (kind=8)    :: xyz(3), denx(nx,ny,nz)


ax = (xyz(1)-x0)/Dx + 1.
ay = (xyz(2)-y0)/Dy + 1.
az = (xyz(3)-z0)/Dz + 1.

ix = int(ax)
iy = int(ay)
iz = int(az)

!write(*,*) ax,ay,az,ix,iy,iz

r = ax - ix 
s = ay - iy 
t = az - iz
!write(*,*) r,s,t
if(ix.ge.nx.or.iy.ge.ny.or.iz.ge.nz.or.ix.le.0.or.iy.le.0.or.iz.le.0) then
 f=0.d0
!write(*,*) ix,iy,iz,nx,ny,nz
else
 call blend_103(r,s,t,denx(ix,iy,iz),denx(ix,iy,iz+1),denx(ix,iy+1,iz),denx(ix,iy+1,iz+1),&
                denx(ix+1,iy,iz),denx(ix+1,iy,iz+1),denx(ix+1,iy+1,iz),denx(ix+1,iy+1,iz+1),f)
!write(*,*) f,denx(ix,iy,iz)
end if

return

end


