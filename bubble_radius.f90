Subroutine Bubble_Radius(iaxis, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xms,yms,zms,il,ir,rl,rr,Lcontrol)

Implicit none
Integer (Kind=4), Intent(in)  :: iaxis, nx,ny,nz,ipx,ipy,ipz
Real    (Kind=8), Intent(in)  :: den(nx,ny,nz), x(nx),y(ny),z(nz)
Real    (Kind=8), Intent(in)  :: xms,yms,zms
Real    (Kind=8), Intent(out) :: rl,rr
Logical         , Intent(out) :: Lcontrol(2)
Logical                       :: Lprint=.false.
Integer (Kind=4), Intent(out) :: il,ir
Integer (Kind=4)              :: naux,niter=100,ix,iy,iz,it,j
Real    (Kind=8)              :: rhomax,aux1,aux2,aux3,step,Eps=1.d-6, Pas=0.1d0,Step_max=0.5d0,h,h1
Real    (Kind=8)              :: denmin=1.5d-2   !  Half maximum at t=0
Real    (Kind=8), Allocatable :: xSpls3(:),ySpls3(:),ASpls3(:),BSpls3(:)
Lcontrol=.false.
If(iaxis.Eq.1)Then
  h1=(x(2)-x(1))      
  h=h1*2.      
  il=0
  Do ix = ipx-1,2,-1
    If(den(ix,ipy,ipz).Gt.den(ix+1,ipy,ipz).And.    &
       den(ix,ipy,ipz).Gt.den(ix-1,ipy,ipz))Then
      If(Abs(Abs(x(ix))-Abs(x(ipx))).Lt.xms.Or.den(ix,ipy,ipz).Lt.denmin)Cycle
      il = ix
      Lcontrol(1)=.true.
      Exit
    EndIf  
  EndDo
  ir = 0
  Do ix = ipx+1,nx-1
    If(den(ix,ipy,ipz).Gt.den(ix+1,ipy,ipz).And.    &
       den(ix,ipy,ipz).Gt.den(ix-1,ipy,ipz))Then
      If(Abs(Abs(x(ix))-Abs(x(ipx))).Lt.xms.Or.den(ix,ipy,ipz).Lt.denmin)Cycle
      ir = ix
      Lcontrol(2)=.true.
      Exit
    EndIf  
  EndDo
  If(.Not.Lcontrol(1).Or..Not.Lcontrol(2))Return
  naux = Max(nx,ny,nz)
  Allocate(xSpls3(naux)); Allocate(ySpls3(naux));Allocate(ASpls3(naux));Allocate(BSpls3(naux));
  rhomax = den(il,ipy,ipz)
  naux = ipx - il +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rl....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf
  Do ix=il, ipx
    j = ix - il + 1
    xSpls3(j) = x(ix)
    ySpls3(j) = den(ix,ipy,ipz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rl = x(il) + h
!  rl = (x(ipx)+x(il))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rl,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rl,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rl, f(rl), df(rl)..:",1p,3E15.6)')rl, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rl = rl +Step
    If(Lprint)Write(6,'("Iaxis, rl, f(rl), Step...:",I5,1p,3E15.6)')iaxis,rl, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  If(it.Ge.Niter)Then
    Lcontrol(1)=.false.
    DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
    Return
  EndIf  
  rhomax = den(ir,ipy,ipz)
  naux = ir - ipx +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rr....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf
  Do ix=ipx, ir
    j = ix - ipx + 1
    xSpls3(j) = x(ix)
    ySpls3(j) = den(ix,ipy,ipz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rr = x(ir) - h
!  rr = (x(ir)+x(ipx))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rr,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rr,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rr, f(rr), df(rr)..:",1p,3E15.6)')rr, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rr = rr +Step
    If(Lprint)Write(6,'("Iaxis, rr, f(rr), Step...:",I5,1p,3E15.6)')iaxis,rr, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
  If(it.Ge.Niter)Then
    Lcontrol(2)=.false.
    Return
  EndIf  
ElseIf(Iaxis.Eq.2)Then  
  h1=(y(2)-y(1))      
  h=h1*2.      
  il=0
  Do iy = ipy-1,2,-1
    If(den(ipx,iy,ipz).Gt.den(ipx,iy+1,ipz).And.    &
       den(ipx,iy,ipz).Gt.den(ipx,iy-1,ipz))Then
      If(Abs(Abs(y(iy))-Abs(y(ipy))).Lt.yms.Or.den(ipx,iy,ipz).Lt.denmin)Cycle
      il = iy
      Lcontrol(1)=.true.
      Exit
    EndIf  
  EndDo
  ir = 0
  Do iy = ipy+1,ny-1
    If(den(ipx,iy,ipz).Gt.den(ipx,iy+1,ipz).And.    &
       den(ipx,iy,ipz).Gt.den(ipx,iy-1,ipz))Then
      If(Abs(Abs(y(iy))-Abs(y(ipy))).Lt.yms.Or.den(ipx,iy,ipz).Lt.denmin)Cycle
      ir = iy
      Lcontrol(2)=.true.
      Exit
    EndIf  
  EndDo
  If(.Not.Lcontrol(1).Or..Not.Lcontrol(2))Return
  naux = Max(nx,ny,nz)
  Allocate(xSpls3(naux)); Allocate(ySpls3(naux));Allocate(ASpls3(naux));Allocate(BSpls3(naux));
  rhomax = den(ipx,il,ipz)
  naux = ipy - il +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rl....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf  
  Do iy=il, ipy
    j = iy - il + 1
    xSpls3(j) = y(iy)
    ySpls3(j) = den(ipx,iy,ipz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rl =y(il) + h
!  rl = (y(ipy)+y(il))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rl,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rl,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rl, f(rl), df(rl)..:",1p,3E15.6)')rl, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rl = rl +Step
    If(Lprint)Write(6,'("Iaxis, rl, f(rl), Step...:",I5,1p,3E15.6)')iaxis,rl, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  If(it.Ge.Niter)Then
    Lcontrol(1)=.false.
    DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
    Return
  EndIf  
  rhomax = den(ipx,ir,ipz)
  naux = ir - ipy +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rr....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf  
  Do iy=ipy, ir
    j = iy - ipy + 1
    xSpls3(j) = y(iy)
    ySpls3(j) = den(ipx,iy,ipz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rr = y(ir) - h
!  rr = (y(ir)+y(ipy))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rr,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rr,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rr, f(rr), df(rr)..:",1p,3E15.6)')rr, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rr = rr +Step
    If(Lprint)Write(6,'("Iaxis, rr, f(rr), Step...:",I5,1p,3E15.6)')iaxis,rr, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
  If(it.Ge.Niter)Then
    Lcontrol(2)=.false.
    Return
  EndIf
ElseIf(Iaxis.Eq.3)Then
  h1=(z(2)-z(1))      
  h=h1*2.      
  il=0
  Do iz = ipz-1,2,-1
    If(den(ipx,ipy,iz).Gt.den(ipx,ipy,iz+1).And.    &
       den(ipx,ipy,iz).Gt.den(ipx,ipy,iz-1))Then
      If(Abs(Abs(z(iz))-Abs(z(ipz))).Lt.zms.Or.den(ipx,ipy,iz).Lt.denmin)Cycle
      il = iz
      Lcontrol(1)=.true.
      Exit
    EndIf  
  EndDo
  ir = 0
  Do iz = ipz+1,nz-1
    If(den(ipx,ipy,iz).Gt.den(ipx,ipy,iz+1).And.    &
       den(ipx,ipy,iz).Gt.den(ipx,ipy,iz-1))Then
      If(Abs(Abs(z(iz))-Abs(z(ipz))).Lt.zms.Or.den(ipx,ipy,iz).Lt.denmin)Cycle
      ir = iz
      Lcontrol(2)=.true.
      Exit
    EndIf  
  EndDo
  If(.Not.Lcontrol(1).Or..Not.Lcontrol(2))Return
  naux = Max(nx,ny,nz)
  Allocate(xSpls3(naux)); Allocate(ySpls3(naux));Allocate(ASpls3(naux));Allocate(BSpls3(naux));
  rhomax = den(ipx,ipy,il)
  naux = ipz - il +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rl....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf  
  Do iz=il, ipz
    j = iz - il + 1
    xSpls3(j) = z(iz)
    ySpls3(j) = den(ipx,ipy,iz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rl =z(il) + h
!  rl = (z(ipz)+z(il))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rl,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rl,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rl, f(rl), df(rl)..:",1p,3E15.6)')rl, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rl = rl +Step
    If(Lprint)Write(6,'("Iaxis, rl, f(rl), Step...:",I5,1p,3E15.6)')iaxis,rl, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  If(it.Ge.Niter)Then
    Lcontrol(1)=.false.
    DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
    Return
  EndIf  
  rhomax = den(ipx,ipy,ir)
  naux = ir - ipz +1
  If(Lprint)Then
    Write(6,'(" Determinacio de rr....")')
    Write(6,'("naux, rhomax...:,",I10,1p,E15.6)')naux,rhomax
    Write(6,'("x, y(Spls3)....:,")')
  EndIf  
  Do iz=ipz, ir
    j = iz - ipz + 1
    xSpls3(j) = z(iz)
    ySpls3(j) = den(ipx,ipy,iz)
    If(Lprint)Write(6,'(I5,1p,2E15.6)')j,xSpls3(j), ySpls3(j)
  EndDo
  aux3=rhomax*0.5d0
  Call Spls3(xSpls3,ySpls3,naux,aux1,aux2,1,ASpls3,BSpls3,0,0)
  rr = z(ir) - h
!  rr = (z(ir)+z(ipz))*0.5
  Do it = 1, Niter
    Call Spls3(xSpls3,ySpls3,naux,rr,aux1,1,ASpls3,BSpls3,1,1)
    Call Spls3(xSpls3,ySpls3,naux,rr,aux2,1,ASpls3,BSpls3,2,1)
    aux1 = aux1 - aux3
    If(Lprint)Write(6,'(" rr, f(rr), df(rr)..:",1p,3E15.6)')rr, aux1, aux2
    Step = -aux1/aux2
    If(Abs(Step).Gt.Step_max)Step=Step*Pas
    rr = rr +Step
    If(Lprint)Write(6,'("Iaxis, rr, f(rr), Step...:",I5,1p,3E15.6)')iaxis,rr, aux1, Step
    If(Abs(Step).Lt.Eps)Exit    
  Enddo
  DeAllocate(xSpls3); DeAllocate(ySpls3);DeAllocate(ASpls3);DeAllocate(BSpls3)
  If(it.Ge.Niter)Then
    Lcontrol(2)=.false.
    Return
  EndIf
Else
  Write(6,'(" Aquest valor de Iaxis no és vàlid...:",I10)')Iaxis      
  Stop 'From Bubble_radius (001)'
EndIf  
Return
End
