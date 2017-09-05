!............................................................
!...                      Subroutine derden               ...
!............................................................

subroutine denpart()

use Para_derivnD
use deriva     ! (npd,dxden,dyden,dzden,icon)
use grid       ! (nx,ny,nz,hx,hy,hz,nxyz,dxyz)
use rho        ! (den)
use util1      ! (nn,mmx,iw)
use work1      ! (nn,mmx,iw)
use impur      
use field


implicit none
Integer  (Kind=4) :: i,ix,iy,iz
Real (kind=8)  :: xmin, ymin, zmin,hxyz

xmin = x(1); ymin = y(1); zmin = z(1)
hxyz = hx*hy*hz

denx = 0.0d0
do i = 1, N_par
  ix = 1.5d0 + (rimp(i,1)-xmin)/hx
  iy = 1.5d0 + (rimp(i,2)-ymin)/hy
  iz = 1.5d0 + (rimp(i,3)-zmin)/hz
  if(iz.Ge.1.And.iz.Le.nz.And.iy.Ge.1.And.iy.Le.ny.And.ix.Ge.1.And.ix.Le.nx)Then
    denx(ix,iy,iz) = denx(ix,iy,iz) + 1.0d0
  EnDif
EndDo
denx = denx/N_Par/hxyz
Do i = 1, ironx

  Call ironing(denx,nx,ny,nz)

EndDo

psix = Dsqrt(denx)


  Call derivnD(2,nn,hx,1,psix,sto1,Icon)
  Call derivnD(2,nn,hy,2,psix,sto2,Icon)
  Call derivnD(2,nn,hz,3,psix,sto3,Icon)
  hpsix = -h2o2mx*(sto1 + sto2 + sto3)

return
end
