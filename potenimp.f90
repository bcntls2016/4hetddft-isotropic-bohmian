subroutine potenimp()
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use classicimp 
use impur
use field
use grid
use work1
implicit none

!...................................
!... First, update uimp (for He) ...
!...................................

call denpart()
call fftfw_denx()
wk1 = wk1*Vq
call fftbk_uimp()
end subroutine potenimp

!---------------------------------------------------------------------------!

subroutine forceimp()
! Esta rutina calcula el potencial sentido por
! el He debido a la impureza. Sabiendo su posicion,
! evalua el potencial de interacion en la malla de trabajo.
use Para_derivnD
use classicimp 
use deriva
 use field
use work1
use util1
use grid
use rho
use impur
implicit none
Real (kind=8) :: xmin, ymin, zmin
Integer (Kind=4) :: ix,iy,iz
xmin=x(1); ymin=y(1); zmin=z(1)

  call fftfw_x()
  call fftfw_y()
  call fftfw_z()
  Wk1 = -wk1*vq
  Wk2 = -wk2*vq
  Wk3 = -wk3*vq 
  call fftbk_x()
  call fftbk_y()
  call fftbk_z()

  Call derivnD(1,nn,hx,1,psix,sto1,Icon)
  Call derivnD(1,nn,hy,2,psix,sto2,Icon)
  Call derivnD(1,nn,hz,3,psix,sto3,Icon)

  Call derivnD(1,nn,hx,1,hpsix,sto4,Icon)
  Call derivnD(1,nn,hy,2,hpsix,sto5,Icon)
  Call derivnD(1,nn,hz,3,hpsix,sto6,Icon)
  sto7 =0.0d0
  Do iz=1, nz
    Do iy=1, ny
      Do ix=1, nx
        If(denx(ix,iy,iz).Gt.0.0d0)sto7(ix,iy,iz)= 1.0d0/denx(ix,iy,iz)
      EndDo
    EndDo
  EndDo
  fx = (fx + sto7*(sto1*hpsix - psix*sto4))/mAg_u
  fy = (fy + sto7*(sto2*hpsix - psix*sto5))/mAg_u
  fz = (fz + sto7*(sto3*hpsix - psix*sto6))/mAg_u

Call interpola3(nx,ny,nz,xmin,ymin,zmin,hx,hy,hz,N_par,rimp,fx,fy,fz,aimp)

end subroutine forceimp


subroutine interpola3(nx,ny,nz,x0,y0,z0,Dx,Dy,Dz,N_par,rimp,fx,fy,fz,aimp)

!use constants
implicit none

real (kind=8)    :: x0,y0,z0,Dx,Dy,Dz,  r,s,t, ax,ay,az,Ddx, Ddz, xx, zz, aax,aay, aaz
integer (kind=4) :: nx,ny,nz,ix,iy,iz,N_par,i,j,Nf=100
real (kind=8)    :: rimp(N_par,3), aimp(N_par,3)
real (kind=8)    :: fx(nx,ny,nx), fy(nx,ny,nz), fz(nx,ny,nz)

Real (Kind=8) :: aux1,aux2,aux3,aux4




Do i=1, N_par

  ax = (rimp(i,1)-x0)/Dx+1
  ay = (rimp(i,2)-y0)/Dy+1
  az = (rimp(i,3)-z0)/Dz+1

  ix = int(ax)
  iy = int(ay)
  iz = int(az)

  r = ax - ix
  s = ay - iy
  t = az - iz
  if(ix.ge.nx.or.iy.ge.ny.or.iz.ge.nz.or.ix.le.0.or.iy.le.0.or.iz.le.0) then
    aimp(i,:)=0.d0
  else
    call blend_103(r,s,t,fx(ix,iy,iz),fx(ix,iy,iz+1),fx(ix,iy+1,iz),fx(ix,iy+1,iz+1),&
                fx(ix+1,iy,iz),fx(ix+1,iy,iz+1),fx(ix+1,iy+1,iz),fx(ix+1,iy+1,iz+1),aimp(i,1))
    call blend_103(r,s,t,fy(ix,iy,iz),fy(ix,iy,iz+1),fy(ix,iy+1,iz),fy(ix,iy+1,iz+1),&
                fy(ix+1,iy,iz),fy(ix+1,iy,iz+1),fy(ix+1,iy+1,iz),fy(ix+1,iy+1,iz+1),aimp(i,2))
    call blend_103(r,s,t,fz(ix,iy,iz),fz(ix,iy,iz+1),fz(ix,iy+1,iz),fz(ix,iy+1,iz+1),&
                fz(ix+1,iy,iz),fz(ix+1,iy,iz+1),fz(ix+1,iy+1,iz),fz(ix+1,iy+1,iz+1),aimp(i,3))
  end if
EndDo


return
aux1 = Minval(rimp(:,1))
aux2 = Maxval(rimp(:,1))
aux3 = Minval(rimp(:,3))
aux4 = Maxval(rimp(:,3))

Ddx = (aux2-aux1)/(nf-1)
Ddz = (aux4-aux3)/(nf-1)



Open(Unit=51,File='fuerzas.dat')
Write(51,'("#x,z,aax,aay,aaz")')
do i = 1, Nf
  xx = aux1 + (i-1)*Ddx
  do j = 1, Nf
    zz = aux3 + (j-1)*Ddz

  ax = (xx-x0)/Dx+1
  ay = -y0/Dy+1
  az = (zz-z0)/Dz+1

  ix = int(ax)
  iy = int(ay)
  iz = int(az)

  r = ax - ix
  s = ay - iy
  t = az - iz
!  if(ix.ge.nx.or.iy.ge.ny.or.iz.ge.nz.or.ix.le.0.or.iy.le.0.or.iz.le.0) then
!    aimp(i,:)=0.d0
!  else
    call blend_103(r,s,t,fx(ix,iy,iz),fx(ix,iy,iz+1),fx(ix,iy+1,iz),fx(ix,iy+1,iz+1),&
                fx(ix+1,iy,iz),fx(ix+1,iy,iz+1),fx(ix+1,iy+1,iz),fx(ix+1,iy+1,iz+1),aax)
    call blend_103(r,s,t,fy(ix,iy,iz),fy(ix,iy,iz+1),fy(ix,iy+1,iz),fy(ix,iy+1,iz+1),&
                fy(ix+1,iy,iz),fy(ix+1,iy,iz+1),fy(ix+1,iy+1,iz),fy(ix+1,iy+1,iz+1),aay)
    call blend_103(r,s,t,fz(ix,iy,iz),fz(ix,iy,iz+1),fz(ix,iy+1,iz),fz(ix,iy+1,iz+1),&
                fz(ix+1,iy,iz),fz(ix+1,iy,iz+1),fz(ix+1,iy+1,iz),fz(ix+1,iy+1,iz+1),aaz)
!  end if

    Write(51,'(1p,5E15.6)')xx,zz,aax,aay,aaz

  Enddo
  Write(51,*)
Enddo
Write(6,'("Aceleraciones maximas..:",1p,3E15.6)')MaxVal(aimp(:,1)), MaxVal(aimp(:,2)), MaxVal(aimp(:,3))
Write(6,'("Aceleraciones minimas..:",1p,3E15.6)')MinVal(aimp(:,1)), MinVal(aimp(:,2)), MinVal(aimp(:,3))
Stop

end

subroutine blend_103 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, x )
!
!*******************************************************************************
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
  implicit none
!
  real(kind=8)    :: r
  real(kind=8)    :: s
  real(kind=8)    :: t
  real(kind=8)    :: x
  real(kind=8)    :: x000
  real(kind=8)    :: x001
  real(kind=8)    :: x010
  real(kind=8)    :: x011
  real(kind=8)    :: x100
  real(kind=8)    :: x101
  real(kind=8)    :: x110
  real(kind=8)    :: x111
!
!  Interpolate the interior point.
!
  x = &
    1.0E+00     * ( + x000 ) &
    + r         * ( - x000 + x100 ) &
    +     s     * ( - x000        + x010 ) &
    +         t * ( - x000               + x001 ) &
    + r * s     * ( + x000 - x100 - x010                      + x110 ) &
    + r     * t * ( + x000 - x100        - x001        + x101 ) &
    +     s * t * ( + x000        - x010 - x001 + x011 ) &
    + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
! write(*,*) x,x000
  return
end
