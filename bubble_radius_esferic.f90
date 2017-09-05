Subroutine Bubble_Radius_esferic(den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xms,yms,zms,nr,rr,denr,ir,rb,Lcontrol)

Implicit none
Integer (Kind=4), Intent(in)  :: nr, nx,ny,nz,ipx,ipy,ipz
Real    (Kind=8), Intent(in)  :: den(nx,ny,nz), x(nx),y(ny),z(nz)
Real    (Kind=8), Intent(out) :: denr(nr),rr(nr)
Real    (Kind=8), Intent(in)  :: xms,yms,zms
Real    (Kind=8), Intent(out) :: rb
Logical         , Intent(out) :: Lcontrol(2)
Logical                       :: Lprint=.false.
Integer (Kind=4), Intent(out) :: ir
Integer (Kind=4), Allocatable :: norm(:)
Integer (Kind=4)              :: naux,niter=100,ix,iy,iz,it,i,j
Real    (Kind=8)              :: rhomax,aux1,aux2,aux3,den_max,r_max,Eps=1.d-6,Step,xnorm
Real    (Kind=8)              :: Pi,rmax,dr,dxyz,hx,hy,hz,x2,y2,z2,r,r2,s1=0.d0,sn=0.d0
Real    (Kind=8)              :: r_min=2.d0, denmin=1.1d-2   !  Half maximum at t=0
Real    (Kind=8), Allocatable :: xSpline(:),ySpline(:),ASpline(:),BSpline(:), CSpline(:), Dspline(:)
Data Pi/3.141592653589793238462643D0/
Lcontrol=.false.
hx = x(2) - x(1)
hy = y(2) - y(1)
hz = z(2) - z(1)
rmax=Dsqrt((-x(1)-x(ipx))**2+(-y(1)-y(ipy))**2+(-z(1)-z(ipz))**2)
dr = rmax/(nr-1)
If(Lprint)Then
  Write(6,'(" Pas de integracio en r: ",1p,E15.6)')dr
  Write(6,'(" rmax..................: ",1p,E15.6)')rmax
EndIf  
dxyz=hx*hy*hz
xnorm=sum(den)*dxyz
If(Lprint)Write(6,'(" Nombre de partícules.....:",F10.3)')xnorm
Allocate(norm(nr))
Allocate(Aspline(nr)); Allocate(Bspline(nr)); Allocate(Cspline(nr)); Allocate(Dspline(nr));

denr = 0.d0
norm = 0
If(Lprint)Then
  Write(6,'(" ipx, x(ipx)...:",I5,1p,E15.6)')ipx,x(ipx)      
  Write(6,'(" ipy, y(ipy)...:",I5,1p,E15.6)')ipy,x(ipy)      
  Write(6,'(" ipz, z(ipz)...:",I5,1p,E15.6)')ipz,z(ipz)      
EndIF        
Do iz=1,nz
  z2=(z(iz)-z(ipz))**2
  Do iy=1,ny
    y2=(y(iy)-y(ipy))**2
    Do ix=1,nx
      x2=(x(ix)-x(ipx))**2
      r=Dsqrt(x2+y2+z2)
      ir= (r/dr + 1.5)
      If(ir.Gt.nr.Or.ir.Lt.1)Then
        Write(6,'("El ir se ha sortit de mare...:",I5)')ir
        Write(6,'("r............................:",1p,E15.6)')r
        Stop '0001'
      Endif
      denr(ir)=denr(ir) + den(ix,iy,iz)
      norm(ir)=norm(ir)+1
    Enddo
  Enddo
Enddo
Do ir=1,nr
  if(norm(ir).Ne.0)denr(ir)=denr(ir)/norm(ir)
Enddo
Deallocate(norm)
denmin = MaxVal(denr)*0.75d0
aux1=0.d0
Do ir=1,nr
  r=(ir-1)*dr
  rr(ir) = r
  r2=r**2
  aux1 = aux1 +r2*denr(ir)
Enddo
aux1=aux1*dr*4.0d0*pi
denr = denr*(xnorm/aux1)
If(Lprint)Write(6,'(" Nombre de partícules(mitjana esferica: Força bruta).....:",F10.3)')aux1
Call Spline(rr,denr,ASpline,BSpline,CSpline,DSpline,S1,Sn,nr)
aux1 = Max(xms,yms,zms)
i = (aux1/dr + 1.5)
Do j = i, nr-1
  If(denr(j).Gt. denr(j-1).And.denr(j).Gt.denr(j+1))Then
    If(rr(j).Lt.r_min.Or.denr(j).Lt.denmin)Cycle
    ir = j
    Lcontrol(1) = .true.
    Exit
  EndIf        
EndDo  
If(j.Eq.(nr-1).Or.j.Gt.(nr-1))Then
  Lcontrol(1)=.false.      
  If(Lprint)Write(6,'("We cannot found the maximum density")')
  DeAllocate(Aspline); DeAllocate(Bspline); DeAllocate(Cspline); DeAllocate(Dspline);
  Return
EndIf        
j=ir
r_max = rr(j)
den_max = denr(j)
If(Lprint)Then
  Write(6,'(" r(j-1), r(j), r(j+1)...............:",1p,3E15.6)')rr(j-1), rr(j), rr(j+1)
  Write(6,'(" denrr(j-1), denr(j), denr(j+1)....:",1p,3E15.6)')denr(j-1), denr(j), denr(j+1)
EndIf
aux1 = Max(xms,yms,zms)
aux2 = den_max*0.5d0
it = (aux1/dr + 1.5)
Do j = ir-1, it+1, -1
  If((denr(j-1)-aux2)*(denr(j+1)-aux2).Lt.0.d0)Then
    i = j
    Lcontrol(2) = .true.
    Exit
  EndIf        
EndDo  
!If(j.Eq.ir.Or.j.Gt.(ir+1))Then
If(j.Eq.it.Or.j.Eq.(it-1))Then
  Lcontrol(2)=.false.      
  If(Lprint)Write(6,'("We cannot found the density=den_max/2(1)")')
  DeAllocate(Aspline); DeAllocate(Bspline); DeAllocate(Cspline); DeAllocate(Dspline);
  Return
EndIf        
rb = rr(i)
Do it = 1, Niter
  aux1 = ((DSpline(i)*rb + CSpline(i))*rb + BSpline(i))*rb + Aspline(i)
  aux3 = (3.*DSpline(i)*rb + 2.*CSpline(i))*rb + BSpline(i) 
  Step = -(aux1 - aux2)/aux3
  rb = rb + Step
  If(Lprint)Write(6,'("rb, step...",1p,2E15.6)')rb, step
  If(Abs(Step).Lt.Eps)Exit
EndDo
If(it.Eq.Niter.Or.it.Eq.(Niter+1))Then
  Lcontrol(2) = .false.      
  If(Lprint)Write(6,'("We cannot found the density=den_max/2(2)")')
  DeAllocate(Aspline); DeAllocate(Bspline); DeAllocate(Cspline); DeAllocate(Dspline);
  Return
EndIf  
If(Lprint)Write(6,'("Bubbe radius....:",1p,E15.6)')rb
DeAllocate(Aspline); DeAllocate(Bspline); DeAllocate(Cspline); DeAllocate(Dspline);
Return
End
