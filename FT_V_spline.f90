Module For_FT_V_spline
  Integer (Kind=4), Save :: N_p = 0, N_FT_V_spline=100001                
  Real    (Kind=8), Save :: r_max = 2000.d0, S1 = 0.d0, Sn = 0.d0
  Real    (Kind=8), Save :: Eps = 1.d-8, small_q=1.d-2, dr_FT_V_Spline=0.5d0
  Real    (Kind=8), Save, Allocatable ::  P_p(:), Vpp(:), A_Vpp(:), B_Vpp(:), C_Vpp(:), D_Vpp(:)
  Character (Len=80), Save  :: Selec_FT_V_spline
  Logical :: Lprint_FT_V_spline=.false., Lparabolic_fit_FT_V_spline=.false., Lredefine_small_q=.true.
End Module For_FT_V_spline
!
!  Funció que ens dona la T.F. d'una funció tipus L.J. (mirar el directori test/test_FFT_Real_to_Real pel seu funcionament practic)
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!

Double Precision Function FT_V_spline(Pp,Pmax,H_p,Selec,r_cutoff,f_cutoff) !

!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!                                                          
!       Pp: Moment on volem saber la T.F.
!     Pmax: Màxim moment on voldrem que es calculi la T.F.
!      H_P: Pas amb que calcularem la T.F. en el espai de moments
!    Seelc: Nom amb el que selecionarem la función en el espai de coordenades
! r_cutoff: Punt de tall per sota del qual soposarem que función es constant i val f_cutoff
! f_cutoff: Valor de la funcio per r=r_cutoff
!
Use For_FT_V_spline
Implicit Real*8(A-H,O-Z)
Real  (Kind=8) :: Pp, Pmax, H_p, r_cutoff, f_cutoff
Character (Len=80)  :: Selec 
Intent (In) Pp, Pmax, H_p, Selec, r_cutoff, f_cutoff

Real (Kind=8), Allocatable :: X(:), Y(:), A(:), B(:), C(:), D(:)

Data Pi/3.141592653589793238462643D0/
!
!  Nombre punts per fer els splines
!
N = N_FT_V_spline
If(Lprint_FT_V_spline)Then
  Open(Unit=10,File='FT_V_spline.out.p')
Endif

100 Continue

   N_p_Aux = (Pmax/H_p + 1.5)

If(N_p_Aux.Eq.N_p.And.selec.Eq.Selec_FT_V_spline)Then
!
!  Executem l'spline realitzat anteriorment (veure just per sota del Else)
!
  Call Findi(P_p,pp,N_p,k)
  FT_V_spline = A_Vpp(k) + pp*( B_Vpp(k) + pp*( C_Vpp(k) + pp*D_Vpp(k) ) )
  Return
Else
  Selec_FT_V_spline = selec
  N_P = N_p_Aux
!
!  Calculem els valors de Vpp, per realitzar l'spline posteriorment
!
  If(Allocated(P_p))Deallocate(P_p)
  If(Allocated(Vpp))Deallocate(Vpp)
  If(Allocated(A_Vpp))Deallocate(A_Vpp)
  If(Allocated(B_Vpp))Deallocate(B_Vpp)
  If(Allocated(C_Vpp))Deallocate(C_Vpp)
  If(Allocated(D_Vpp))Deallocate(D_Vpp)
  Allocate (P_p(N_p)); Allocate(Vpp(N_p)); Allocate(A_Vpp(N_p)); Allocate(B_Vpp(N_p)); Allocate(C_Vpp(N_p)); Allocate(D_Vpp(N_p));
  If(Allocated(X))Deallocate(X)
  If(Allocated(Y))Deallocate(Y)
  If(Allocated(A))Deallocate(A)
  If(Allocated(B))Deallocate(B)
  If(Allocated(C))Deallocate(C)
  If(Allocated(D))Deallocate(D)

Allocate(X(N))  ! Vectors pels splines
Allocate(Y(N))  !              "
Allocate(A(N))  !              "
Allocate(B(N))  !              "
Allocate(C(N))  !              "
Allocate(D(N))  !              "

!
!  Definicio dels diferents intervals dels splines
!
X(1) = R_cutoff + Eps
X(N) = r_max
!
!  Calculem els splines de la funció, per a posteriorment fer l'integració analítica
!
!  Calculem la derivada segona per millorar l'spline
!

r0 = X(1)
r1 = r0 + dr_FT_V_Spline
r2 = r1 + dr_FT_V_Spline
f0 = Select_Pot(selec,r0,r_cutoff,f_cutoff)
f1 = Select_Pot(selec,r1,r_cutoff,f_cutoff)
f2 = Select_Pot(selec,r2,r_cutoff,f_cutoff)

S1 = (f0 - 2.*f1 + f2)/dr_FT_V_Spline**2

!
!  Cas especial del LJ de la interacció de Orsay-Trento
!
If(Trim(selec).Eq.'LJ_OT')Then
  Lparabolic_fit_FT_V_spline=.false.
Endif

  r0  = X(N)
  rm1 = r0 - dr_FT_V_Spline
  rp1 = r0 + dr_FT_V_Spline
  fm1 = Select_Pot(selec,rm1,r_cutoff,f_cutoff)
  f0  = Select_Pot(selec,r0,r_cutoff,f_cutoff)
  fp1 = Select_Pot(selec,rp1,r_cutoff,f_cutoff)

  Sn = (fm1 - 2.*f0 + fp1)/dr_FT_V_Spline**2

  h  = (X(N)-X(1))/(N-1)
  Do i=1, N
    r = (i-1)*h + X(1)
    x(i)=r
    y(i)=Select_Pot(selec,r,r_cutoff,f_cutoff)
  EndDo
  Call Spline(X,Y,A,B,C,D,S1,Sn,N)

Qpi=4.0d0*Pi
Pio2= Pi*0.5d0
r0  = r_cutoff
f0  = f_cutoff
!
!  Per si volem ajustar a un paràbola invertida fins l'origen
!
If(Lparabolic_fit_FT_V_spline)Then
  r0 = r_cutoff
  f0 = f_cutoff
  r1 = r0 + dr_FT_V_Spline
  f1 = Select_Pot(selec,r1,r_cutoff,f_cutoff)
  aa = (f0-f1)/(r0**2-r1**2)
  cc = (r0**2*f1-r1**2*f0)/(r0**2-r1**2)
!
!  si la paràbola no és invertida hem de parar
!
  If(aa.Gt.0d0)Then
    Write(6,'("From FT_V_Spline: the parabola must be inverted..(a*r^2+b^r+c)..,a,b,c....:",/,1p,3E18.10)')aa,bb,cc
    Stop 'FromFT_V_Spline  004'
  Endif
Endif
rf  = x(N)
!
!  Aqui començema calcular la T.F. analítica amb els coeficients dels diferents splines
!
Do i=1,N_p
  p  = (i-1)*H_p
  P_p(i) = p
  q  = 2.0d0*Pi*P
  q2 = q*q
  q3 = q2*q
  c1oq=1.0d0/q
  c1oq2=c1oq*c1oq
  c1oq3=c1oq2*c1oq
  c1oq4=c1oq2*c1oq2
  c1oq5=c1oq4*c1oq
!
!  Cas especial per p=0
!
  If(P.eq.0.0d0)Then
    If(Lparabolic_fit_FT_V_spline)Then
      Vp =  aa*r0**5*0.2d0 + cc*r0**3/3.d0
    Else
      Vp = r0**3*f0/3.0d0
    Endif
    Fi =0.d0
!
! Integrem analiticament els splines
!
      Do j=1,N-1
        Fi = Fi + a(j)/3.0d0*(x((j+1))**3-x(j)**3)                &
                      + b(j)/4.0d0*(x((j+1))**4-x(j)**4)          &
                      + c(j)/5.0d0*(x((j+1))**5-x(j)**5)          &
                      + d(j)/6.0d0*(x((j+1))**6-x(j)**6)
      EndDo
    Vp = (Vp +  Fi)*Qpi
    If(Lprint_FT_V_spline)Then
!      Write(10,'(F10.3,1p,7E15.6)')p, Vp, Sto1, (Fi(j), j=1, N_intervals)
      Write(10,'(1p,E12.5,2x,7E15.6)')p, Vp
    Endif
    Vpp(i) = Vp
  Else
!
!  Casos per p#0
!
    If(Lparabolic_fit_FT_V_spline)Then
      Vp = c1oq*(Cos(q*r0)*(6.*aa*r0*c1oq2 - aa*r0**3 - cc*r0) +                                 &
                 Sin(q*r0)*((3.*aa*r0**2+cc)*c1oq -6.*aa*c1oq3))
    Else
      Vp = f0*c1oq*(Sin(q*r0)*c1oq - r0*Cos(q*r0))
    Endif
    Fi = 0.0d0
!
! Integrem analiticament els splines
!
      Do j=1, N-1
        rj = x(j)
        rj1 = x(j+1)
        sqrj = Sin(q*rj)
        cqrj = Cos(q*rj)
        sqrj1 = Sin(q*rj1)
        cqrj1 = Cos(q*rj1)
        Fi = Fi + (A(j) + (B(j) + (C(j) + D(j)*rj)*rj)*rj)*rj*c1oq*cqrj                   &
                      - (A(j) + (B(j) + (C(j) + D(j)*rj1)*rj1)*rj1)*rj1*c1oq*cqrj1        &
                      + (A(j) + (2.*B(j) + (3.*C(j) + 4.*D(j)*rj1)*rj1)*rj1)*c1oq2*sqrj1  &
                      - (A(j) + (2.*B(j) + (3.*C(j) + 4.*D(j)*rj)*rj)*rj)*c1oq2*sqrj      &
                      + (2.*B(j) + (6.*C(j) + 12.*D(j)*rj1)*rj1)*c1oq3*cqrj1              &
                      - (2.*B(j) + (6.*C(j) + 12.*D(j)*rj)*rj)*c1oq3*cqrj                 &
                      + (6.*C(j) + 24.*D(j)*rj)*C1oq4*sqrj                                &
                      - (6.*C(j) + 24.*D(j)*rj1)*C1oq4*sqrj1                              &
                      + 24.*D(j)*c1oq5*(cqrj - cqrj1)
      EndDo
    Vp    = (Vp +Fi)*C1oq*Qpi
    If(Lprint_FT_V_spline)Then
!      Write(10,'(F10.3,1p,7E15.6)')p, Vp, Sto, (Fi(j), j=1, N_intervals)
      Write(10,'(1p,E12.5,2x,7E15.6)')p, Vp
    Endif
    VPP(i) = Vp
  Endif
EndDo
If(Lprint_FT_V_spline)Close(10)
!
!  Abans de fer l'Spline, si ho creiem convenient redefinim el valors de Vpp per valors petits de p
!
If(Lredefine_small_q)Then
  iaux=1
  X(1) = 0.d0
  Y(1) =Vpp(1)
  Do i=2, N_p
    P = P_p(i)
    If(p.Gt.small_q)Then
      iaux = iaux +1
      x(iaux) = p
      y(iaux) = Vpp(i)
    Endif
    If(iaux.Eq.4)Exit
  EndDo
  If(Lprint_FT_V_spline)Write(6,'(" i, iaux...:",2I5)')i, iaux
  S1 = 0.d0
  Sn =0.d0
  Call Spline(X,Y,A,B,C,D,S1,Sn,4)
  Do i=2,N_p
    p= P_p(i)
    If(p.Le.small_q)Then
      If(Lprint_FT_V_spline)Write(6,'(" Valor abans de la seva refinicio.....:",1p,2E15.6)')p,Vpp(i)      
      Vpp(i) = A(1) + p*(B(1) + P*(C(1) + p*D(1)))
      If(Lprint_FT_V_spline)Write(6,'(" Valor despres de la seva refinicio...:",1p,2E15.6)')p,Vpp(i)      
    Else
      Exit
    EndIf
  EndDo
  If(Lprint_FT_V_spline)Then
    Open(Unit=10,File='FT_V_spline.out.p')
    Do i=1, N_p
      Write(10,'(1p,E12.5,2x,7E15.6)')p_p(i), Vpp(i)
    EndDo
    Close(10)
  Endif
EndIf  
!
!  Fem l'spline en l'espài de moments per poder-lo emprar posteriorment (veure les instruccions a partir del 100 Continue)
!
  Call Spline(P_p,Vpp,A_Vpp,B_Vpp,C_Vpp,D_Vpp,S1,Sn,N_p)

!
! Anem al començament per fer l'spline
!
  Go to 100
Endif
Return
End
