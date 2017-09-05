program DFT4HeImpd

!
!                         ******************
!                         *** DFT4HeImpd ***
!                         ******************
!
!

!
! This program computes time dependent Helium_4 drops with/without impurities using
! Density Functional Theory.
! The interation can be Orsay-Paris or Orsay-Trento
! The code is fully 3-dimensional
!
!-----------------------------------------------------------------------------
! Version 0  (alpha)   Barcelona April-19,2004   R. Mayol & M. Pi
! Version 99 (Epsilon) Barcelona April-19,2006   R. Mayol & M. Pi
!-----------------------------------------------------------------------------
!07oct2004/11:30AM  
!06oct2004/11:30AM
!23sep2004/11:30AM
!Feb2005
! (LE FALTA: Ficheros binarios)
!-----------------------------------------------------------------------------


Use For_FT_V_spline
use seleccio_de_potencial
use alphasterm
use deriva
use energies
use rho
use field
use fftmodule
use grid 
use gridk
use classicimp
use rkpc
use impur
use lenard4
use he4
use util1
use work1

implicit none

real (kind=4) :: t0,t1,t2,t3,t4,t5,t6   ! Variables used to printout run time execution.

logical              :: lfilepv         ! T-> print save file when change Paflo parameter.
logical              :: lpaflv=.false.  ! T-> allows change of Paflov coeffient
logical              :: lrkpc=.true.    ! T-> allows to use diferent evolution procedures
logical              :: lrk=.false.     ! T-> allows to only Runge-Kutta method
integer    (kind=4)  :: ndmax=2         ! maxima derivada a calcular
integer    (kind=4)  :: naux            ! Auxiliar variable
integer    (kind=4)  :: nstepp=1        ! Number of 'Paflov parameter'.
integer    (kind=4)  :: ix ,iy ,iz      ! Used as indexs of arrays
integer    (kind=4)  :: ipx,ipy,ipz     ! Used as indexs of arrays
integer    (kind=4)  :: ixl,ixr,iyl,iyr,izl,izr     ! Used to compute the bubble radius
integer    (kind=4)  :: iter,niter      ! Control the number of iterations
integer    (kind=4)  :: pener=50        ! Computes energy each 'pener' iterations
integer    (kind=4)  :: pdenpar=50      ! Writes partial densities each 'pdenpar' iterations
integer    (kind=4)  :: pchem=50        ! Writes partial chemical potential each 'pchem' iter.
integer    (kind=4)  :: ppot=1          ! Computes the potential each 'ppot' iterations
integer    (kind=4)  :: tstgrid
integer    (kind=4)  :: norder          ! Order of Taylor expansion for time evolution
integer    (kind=4)  :: nv_iter         ! Number of iterationsfor the first time evolution
real       (kind=8)  :: norma           ! Use for normalization
real       (kind=8)  :: mimpur          ! Impurity mass in uma
real       (kind=8)  :: rhomax          ! Impurity mass in uma
real       (kind=8)  :: v0x=0., v0y=0., v0z=0.   ! Initial velocity (Ang/ps)
real       (kind=8)  :: xbl, xbr, ybl, ybr, zbl, zbr
real       (kind=8)  :: denxmax         ! Maximum value for density of the impurity
real       (kind=8)  :: r_clust         ! Radius of the helium cluster
real       (kind=8)  :: rimpur=5.5      ! Radius of the impurity
real       (kind=8)  :: p,px2,py2,pz2   ! Temporary variables for momentum values
real       (kind=8)  :: mu4,mu4err      ! Value of Chemical potential and associated error
real       (kind=8)  :: epsx,epsxerr    ! Value of autovalue and associated error
real       (kind=8)  :: errmu4          ! Relative change betwen iteration for chemical potential
real       (kind=8)  :: deltat=0.d0     ! Time step
real       (kind=8)  :: deltatps=0.d0   ! Time step in picoseconds
real       (kind=8)  :: xlamda=0.0d0    ! Monopol costraint
real       (kind=8)  :: xlamdax=0.0d0    ! He Dipol costraint x direction
real       (kind=8)  :: xlamday=0.0d0    ! He Dipol costraint y direction
real       (kind=8)  :: xlamdaz=0.0d0    ! He Dipol costraint z direction
real       (kind=8)  :: xlamdax_x=0.0d0    ! Impurity Dipol costraint x direction
real       (kind=8)  :: xlamdax_y=0.0d0    ! Impurity Dipol costraint y direction
real       (kind=8)  :: xlamdax_z=0.0d0    ! Impurity Dipol costraint z direction
real       (kind=8)  :: eold            ! Auxiliar variables
real       (kind=8)  :: aux,aux1,aux2,aux3             ! Auxiliar variables
real       (kind=8)  :: aux4,aux5,aux6                 ! Auxiliar variables
integer    (kind=4),allocatable  :: nitera(:) ! In wich iteration change to the
real       (kind=8),allocatable  :: paflv(:)  !    corresponging Paflov coeffiecient
real       (kind=8),allocatable  :: xSpls3(:), ySpls3(:), ASpls3(:), BSpls3(:)  ! Vectors to perform a spline
real       (kind=8)  :: cnorm           !    Normalization constant
real       (kind=8)  :: pmod1        !    Work variable for controlling p-values
real       (kind=8)  :: pmod2        !    Work variable for controlling p-values
real       (kind=8)  :: xcm4,ycm4,zcm4,xcmx,ycmx,zcmx ! Center of mass Drop and Impurity
real       (kind=8)  :: xmsx,ymsx,zmsx                ! The impurity mean square radius 
real       (kind=8)  :: distx,disty,distz  ! Distance between center of masses
real       (kind=8)  :: rz,rx,ry,r,rmax,r2,qmax,Hq,FT_V_spline
real       (kind=8)  :: errHe, errimp,errvimp     ! Error evolution (only form Predictor-Corrector-Modificator)
real       (kind=8)  :: Zsurf = -25.d0, Morse_HeTiO2_1D, Morse_HeTiO2_3D, auxn4 ! Position of the surface
real       (kind=8)  :: Ktops=7.63822291d0, time0=0.0d0, time=0.0d0
real       (kind=8)  :: pstoK=1.d0/7.63822291d0, rb
Real       (Kind=8)  :: xi=1.d99, yi=1.d99, zi=1.d99, xf=1.d99, yf=1.d99, zf=1.d99
Real       (Kind=8), Allocatable  :: denr(:), rr(:)

integer    (kind=4)  :: n4=300         ! Number of helium_4 atoms
integer    (kind=4)  :: mode=0          ! Way to start the program (see readen subroutine)
integer    (kind=4)  :: instate=0          ! Way to start the program (see readen subroutine)
integer    (kind=4)  :: iter0=0         ! Starting point for iteration procedure
integer    (kind=4)  :: ncurr,pcurr
integer    (kind=4)  :: nr=501 ! For spherical density average construction
integer    (kind=4)  :: icurr = 0, ir

character  (len=40)  :: title         = 'Helium4 - 3dim.  '
character  (len=60)  :: fileout       = 'DFT4He3d.res'
character  (len=60)  :: filedenin     = 'he4.input.density'
character  (len=60)  :: filedenout    = 'he4.output.density'
character  (len=60)  :: fileimpin     = 'X.input.wf'
character  (len=60)  :: fileimpout    = 'X.output.wf'
character  (len=60)  :: filerimp      = 'rimp.out'
character  (len=60)  :: filevimp      = 'vimp.out'
character  (len=60)  :: fileaimp      = 'aimp.out'
character  (len=60)  :: filebubble    = 'bubble.out'
character  (len=60)  :: file_spherical_bubble    = 'spherical_bubble.out'
character  (len=60)  :: namefile,namefile1

character  (len=23)  :: curvfile
character  (len=3)   :: chariter


logical              :: Lbubble_radius=.false.        ! We will compute the bubble radius each time that we compute the energy
logical              :: Lspherical_bubble=.true.      ! We will compute the bubble radius from an speherical average density
logical              :: lsurf=.true.                  ! include TiO2 surface or not
logical              :: lsurf3D=.true.                ! include TiO2 surface or not
logical              :: Lcontrol(2)                   ! Auxiliar logicval variable
logical              :: LcontrolX,LcontrolY,LcontrolZ ! Auxiliar logicval variable

real       (kind=8)  :: Lambdah,tzmean,tzsurf,rt



!....................Variables read in a NAMELIST statement ..............................

namelist /input/title,fftwplan,nthread,nsfiles,                         &
                fileout,file_spherical_bubble,                          &
                filedenin, filedenout,                                  &
                fileimpin, fileimpout,                                  &
                n4,mode,filerimp,filevimp,fileaimp,filebubble,          &
                nx,ny,nz,xmax,ymax,zmax,Lbubble_radius,                 &
                xc,yc,zc,afermi,Lspherical_bubble,                      &
                eps4,sigma4,core4,l,nq,                                 &
                cp4,cpp4,den4c,alphas,h2o2m4,                           &
                denmin,psimin,npd,ndmax,icon,                           &
                niter,printpot,pchem,irespar,                           &
                pdenpar,pener,ppot,iron,ironx,                          &
                rimpur,mimpur,nr,                                       &
                norder,nv_iter,deltat,lrkpc,lrk,xlamda,lselection,      &
                xlamdax,xlamday,xlamdaz,deltatps,                       &
                xlamdax_x,xlamdax_y,xlamdax_z,iter0,Zsurf,time0,        &
                pcurr,icurr,lsurf,lsurf3D,Lambdah,tzmean,tzsurf,        &
                N_par,v0x,v0y,v0z,xi,xf,yi,yf,zi,zf,                    &
                selec_gs, r_cutoff_gs,  umax_gs

 
!................................ Start main Program ..............................
call timer(t0)

!.............................................
!... Inicializate some numerical constants ...
!.............................................

pi     = 4.0d0*datan(1.0d0) ! Initialization of pi
twopi  = 2.0d0*pi
fourpi = 4.0d0*pi
piq    = pi*pi


!...............................
!... Read  master-input file ...
!...............................
Lambdah = 0.d0
tzmean   = 0.d0
tzsurf   = 1.d0
read(5,nml=input,end=999)
If(xi.Gt.1.d10.Or.                             &
   xf.Gt.1.d10.Or.                              &
   yi.Gt.1.d10.Or.                              &
   yf.Gt.1.d10.Or.                              &
   zf.Gt.1.d10.Or.                              &
   zi.Gt.1.d10      ) Then
   Write(6,'("Has de definirme los limites de la impuerza...")')
   Stop ' Capullo'
Endif
open(10,file="DFT4He3d.namelist.read")
write(10,nml=input)
call flush(10)

!.......................................
!... Define the mass of the impurity ...
!.......................................
mAg_u = mimpur*mp_u


nn(1)  = nx ; nn(2)  = ny ; nn(3)  = nz;                ! Initialize things for PDERG
mmx(1) = nx ; mmx(2) = ny ; mmx(3) = nx ; mmx(4) = ny   ! (NO NOT MOVE THAT!!!!!!!)

!.............................................................
!.. Check if the size of the grid is among the valid values ..
!.............................................................

if(tstgrid(nx).ne.0 ) stop 'SEVERE ERROR: NX is not correct'
if(tstgrid(ny).ne.0 ) stop 'SEVERE ERROR: NY is not correct'
if(tstgrid(nz).ne.0 ) stop 'SEVERE ERROR: NZ is not correct'

!...................................................
!.. Controls Paflov parameters (read and storage ...
!...................................................



close(5)
close(10)

!.........................................
!.. Some consistency check on Delta t  ...
!.........................................
If(deltat.eq.0.d0)then
 if(deltatps.eq.0.d0)then    ! deltat ==  0, deltatps ==  0
  print*,'You must specify either Deltat (in kelvin) or Deltatps (in picosecond)'
  STOP
 else                        ! deltat ==  0, deltatps =/= 0
  deltat = deltatps/7.63822291d0
 endif
Else
 if(deltatps.eq.0.d0)then    ! deltat =/= 0, deltatps ==  0
  deltatps = deltat*7.63822291d0
 else                        ! deltat =/= 0, deltatps =/= 0
  if(deltatps.ne.(deltat*7.63822291d0))then
   print *,'Inconsistent deltat - deltatps'
   STOP
  endif
 endif
Endif
write(*,*)'Time step is ',deltat,' kelvins or ',deltatps,' picoseconds.'


!................................................
!.. Some consistency check on input variables ...
!................................................

nthread=abs(nthread)

Call Init_deriv_p(npd,ndmax,nthread)

! if(mode.eq.2) then
! !   if(.not.limp) then
!      write(6,*) ' '
!      write(6,*) ' Inconsistency input error.'
!      write(6,*) ' '
!      write(6,*) ' If mode=2. MUST be limp=.true.'
!      write(6,*) ' '
!      write(6,*) ' ACTION: ABORT EXECUTION'
!      write(6,*) ' '
!      STOP
! !   end if
! end if

hx    = 2.0d0*abs(xmax)/(nx)  ! Step in x-grid
hy    = 2.0d0*abs(ymax)/(ny)  ! Step in y-grid
hz    = 2.0d0*abs(zmax)/(nz)  ! Step in z-grid

dxyz  = hx*hy*hz              ! Element of volum in real space
nxyz  = nx*ny*nz              ! Total number of points
dVomAg = dxyz/mAg_u

hpx   = 1.0d0/(nx*hx)         ! Step in px-grid
hpy   = 1.0d0/(ny*hy)         ! Step in py-grid
hpz   = 1.0d0/(nz*hz)         ! Step in pz-grid

pmaxx = 1.0d0/(2.0d0*hx)      ! Maximum 'frequency' in X-grid
pmaxy = 1.0d0/(2.0d0*hy)      ! Maximum 'frequency' in Y-grid
pmaxz = 1.0d0/(2.0d0*hz)      ! Maximum 'frequency' in Z-grid


!...............................
!.. Dimensionate main ARRAYS ...
!...............................

call dimen()

!.....................................
!.....................................


!................................
!... Build grid in real space ...
!................................

do ix=1,nx  !.................... Grid X
 x(ix) = -xmax+hx*(ix-1)
end do
do iy=1,ny  !.................... Grid Y
 y(iy) = -ymax+hy*(iy-1)
end do
do iz=1,nz  !.................... Grid  Z
 z(iz) = -zmax+hz*(iz-1)
end do

!....................................
!... Build grid in momentum space ...
!....................................

!.... Build p-grid. In order to use FFTW the grid must
!     start from frequency zero to the maximum and then continue from
!     (-maximum) to zero (negative).

!............................................ grid Px
do ipx=1,nx/2+1
   px(ipx) =        hpx*(ipx-1)
end do
do ipx=nx/2+2,nx
   px(ipx) = -pmaxx+hpx*(ipx-(nx/2+1))
end do

!............................................ grid Py
do ipy=1,ny/2+1
   py(ipy) =        hpy*(ipy-1)
end do
do ipy=ny/2+2,ny
   py(ipy) = -pmaxy+hpy*(ipy-(ny/2+1))
end do

!............................................ grid Pz
do ipz=1,nz/2+1
   pz(ipz) =        hpz*(ipz-1)
end do
do ipz=nz/2+2,nz
   pz(ipz) = -pmaxz+hpz*(ipz-(nz/2+1))
end do

!............................................ Compule modulus of p
do ipz=1,nz
  pz2=pz(ipz)**2
  do ipy=1,ny
    py2=py(ipy)**2
    do ipx=1,nx/2+1
      px2               = px(ipx)**2
      pmod(ipx,ipy,ipz) = sqrt(px2+py2+pz2)
    end do
  end do
end do

pmod1=maxval(pmod)
pmod2=sqrt(pmaxx**2+pmaxy**2+pmaxz**2)


!write(6,*) '    Initialize Linear Interpolation for V_Pi and V_Sig'
call flush(8)
!call potenimpini() ! interpolation + first call to updatepoten

!................................
!... read density or build-it ...
!................................

call readenc(n4,densat4,filedenin,fileimpin,mode,rimpur,r_clust)

!...............................................................test
! do iz=1,nz ; do iy=1,ny ; do ix=1,nx 
!  if(.not.(den(ix,iy,iz).gt.0))print*,ix,iy,iz,den(ix,iy,iz)
! end do ; enddo ; enddo
!...............................................................test


! if(irespar.ne.0) then
!    if(limp) then
       call respar(x,y,z,nx,ny,nz,2,'hedenini','impurdenini',den,denx)
!    else
!      call respar(x,y,z,nx,ny,nz,1,'hedenini','hedenini',den,den)
!    end if
! end if

!....................................
!.. Print-out initial quantities  ...
!....................................

!open(6,file=fileout)

write(6,6010) title

select case(mode)
   case(0) !................................... Start a dynamic form static calcultaion
!       if(limp) then
!          write(6,6111) filedenin,filedenout,fileimpin,fileimpout
!       else
         write(6,6011) filedenin,filedenout
!       end if
   case(2) !................................... Continue a dynamic calculation from a prevous one
      write(6,6013) filedenout,fileimpout
   case default !.............................. Start still not programed.
      write(6,*) ' '
      write(6,*) ' The variable mode has no acceptable value.'
      write(6,*) ' '
      stop
end select

!...............................................................
!................................ Write the input parameters ...
!...............................................................

write(6,6018) nthread,niter
if(mode.ne.0) then
  write(6,6020) n4,r_clust
else
  write(6,6025) n4
end if
write(6,6030) nx,ny,nz,hx, hy, hz, x(1) ,y(1) ,z(1) ,x(nx),y(ny),z(nz)
write(6,6035)          hpx,hpy,hpz,px(1),py(1),pz(1),pmaxx,pmaxy,pmaxz,&
                       pmod1,pmod2
write(6,6037) cp4,cpp4,den4c,alphas,l,den0s,h2o2m4


!...............
!.. Impurity ...
!...............

    h2o2mx = h2o2m4*(4.002602d0/mimpur)

! if(limp) then
!    if(mimpur.lt.1.d-10) then
!      print *,'********************'
!      print *,'*** SEVERE ERROR ***'
!      print *,'********************'
!      print *,' '
!      print *,'If LIMP=.TRUE. Mimpur cannor be zero...'
!      print *,'  '
!      print *,'PROGRAM CANCELLED'
!    end if
!    h2o2mx = h2o2m4*(4.002602d0/mimpur)
!    write(6,6138) h2o2mx
!    cnorm = sum(denx)*dxyz
!    write(6,*) ' '
!    write(6,*) '    Wave function normalization:   ',cnorm
!    write(6,*) ' '
!    write(6,*) '    Calculation with the impurity: ',elem
!    write(6,fmt='(''         Atomic mass '',F10.3)') mimpur
!    write(6,*) ' '
!    write(6,*) '    Initial position of the impurity:'
!    write(6,fmt='(''            X_imp = '',F10.5,'' A'')') ximp
!    write(6,fmt='(''            Y_imp = '',F10.5,'' A'')') yimp
!    write(6,fmt='(''            Z_imp = '',F10.5,'' A'')') zimp
!    write(6,*) ' '
! else
!    write(6,*) '    Calculation without impurities.'
! end if

!...................................................................
!... Compute the FFT of Lennard-Jones                            ...
!... Prepara \alpha_s term in case of Orsay-Trento Interaction.  ...
!...................................................................

select case(core4)
   case('OP ')
     h4=h4op
     write(6,*) '    Use Orsay-Paris-Barcelona Interaction.'
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OT ')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction. (ONLY CORE)'
     write(6,*) '    Do not calculate Alpha_s in the field neither energy.'
     write(6,6040) core4,h4,eps4,sigma4,b4
   case('OTC')
     h4=h4ot
     write(6,*) '    Use Orsay-Trento Interaction.'
     write(6,*) '    Full Orsay-Trento calculation. (Field and Energy)'
     write(6,6040) core4,h4,eps4,sigma4,b4
!     allocate( denalf(nx  ,ny,nz))                                                            
!     allocate(  falfs(nx  ,ny,nz))
!     allocate(kalfs(nx/2+1,ny,nz))
!     allocate(intxalf(nx  ,ny,nz))
!     allocate(intyalf(nx  ,ny,nz))
!     allocate(intzalf(nx  ,ny,nz))
!     allocate(ualphas(nx  ,ny,nz))
   case default
     print *,' ***************** WARNING ************************'
     print *,' '
     print *,' I do not know how to work with the core: ',core4
     print *,' '
     print *,' **************************************************'
     print *,' '
     STOP ' ... The program stops due to a severe error.'
end select
     allocate( denalf(nx  ,ny,nz))                                                            
     allocate(  falfs(nx  ,ny,nz))
     allocate(kalfs(nx/2+1,ny,nz))
     allocate(intxalf(nx  ,ny,nz))
     allocate(intyalf(nx  ,ny,nz))
     allocate(intzalf(nx  ,ny,nz))
     allocate(ualphas(nx  ,ny,nz))



!...............................
!... Prepare plans for FFTWs ...
!...............................
write(6,*) '    Initialize Plans for FFT.'
call fftini(nx,ny,nz)

!...........................................................
!... Form  factor for Lennard-Jones and for the impurity ...
!...........................................................

write(6,*) '    Compute the FFT of the kernel of Lennard-Jones integrals.'

call fforma(core4,b4,eps4,sigma4,h4,nx,ny,nz,pmod,fvlj4)

!...........................................................
!... Form factor for Impurity potential and prepare      ...
!... arrays for the inclusion of the impurity            ...
!...                                                     ...
!.......................................................................
!... If (mode=2), contnue a dynamic calculation from a previuos one  ...
!.......................................................................

   qmax = Sqrt(pmaxx**2+pmaxy**2+pmaxz**2)
   Hq   = qmax/(Nq-1)
   Do iz=1, nz
     Do iy=1, ny
       Do ix=1, nx/2 + 1
         p=pmod(ix,iy,iz)
         Vq(ix,iy,iz) = cmplx(FT_V_spline(p,qmax,Hq,Selec_gs,r_cutoff_gs,umax_gs))
       EndDo
     EndDo
   EndDo

   open(32,file='Interpolated.Impurity.Potential.new.dat')

do iz=nz/2+2,nz
   write(32,32000) pz(iz),real(vq(1,1,iz))
end do
do iz=1,nz/2+1
   write(32,32000) pz(iz),real(vq(1,1,iz))
end do

32000 format(5x,0P,F10.5,3x,1P,E13.5)
close(32)
!   call fftfw(den ,fden)
!    call fftfw_denx()
!   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
!     wk1(ix,iy,iz) = fden(ix,iy,iz)*vq(ix,iy,iz)
!     wk1(ix,iy,iz) = wk1(ix,iy,iz)*vq(ix,iy,iz)
!   end forall
!   call fftbk_uimp()                  ! Integral{denx*V_x}
!   call fftbk(wk2,potx4)                  ! Integral{|Psi_x|**2 V_X}

!........................................
!.. Initialize coarse-graining kernel ...
!........................................

write(6,*) '    Initialize Coarse-graining kernel.'
call initcg(h4,wcgk)


if(core4.eq.'OTC') then
   forall(ix=1:nx/2+1,iy=1:ny,iz=1:nz)
      kalfs(ix,iy,iz) = exp(-(pi*l*pmod(ix,iy,iz))**2)
   end forall
end if

!
!    We compute the perturbate wave function
!
! Initial velocity : Altough xlamdax_i is in K units, xlamda
! is introduced in Angstrom/picosecond for the sake of lazyness.
auxn4 = sum(den)*dxyz

if(xlamda.ne.0.d0)then
  write(*,*)'xlambda not equal zero'
  xlamdaz = xlamda*7.63822291d0/(2.d0*h2o2m4)
!  xlamdaz = xlamda*0.630247988
endif

If(mode.eq.0)then
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         aux=x(ix)*xlamdax       &
            +y(iy)*xlamday       &
            +z(iz)*xlamdaz
         psi(ix,iy,iz) = sqrt(den(ix,iy,iz)) &
                       * cmplx(cos(aux),sin(aux))
       end do
     end do
   end do

  Call sorteo(denx,nx,ny,nz,x,y,z,rimp,N_par,xi,yi,zi,xf,yf,zf)

  v0x =  v0x*Ktops ! Because it is given in A/ps, we have to transform to A*K
  v0y =  v0y*Ktops ! Because it is given in A/ps, we have to transform to A*K
  v0z =  v0z*Ktops ! Because it is given in A/ps, we have to transform to A*K
  vimp(:,1) =  v0x
  vimp(:,2) =  v0y
  vimp(:,3) =  v0z

     aux1 = sum(rimp(:,1))/N_par
     aux2 = sum(rimp(:,2))/N_par
     aux3 = sum(rimp(:,3))/N_par
     Write(6,'("Impurity C.M.,<x>,<y>,<z>..:",1p,3E15.6)')aux1,aux2,aux3

      namefile='density.000.dat'
      namefile1='densityx.000.dat'

       call printoutc(0.d0,3,namefile,namefile1,psi,nx,ny,nz,hx,hy,hz, &
                    xmax,ymax,zmax,rimp,vimp,N_par,        &
                    deltat,iter)

Endif

open(81,file=filerimp)
open(82,file=filevimp)
open(83,file=fileaimp)


!.................................
!.. First call to total energy ...
!.................................

den=Conjg(psi)*psi

call potenimp()
call poten()              ! First Potential  (for Lagrange Equation)
call forceimp()

write(81,'("# Tiempo(ps), x(AA), y(AA), z(AA)")')
     aux1 = sum(rimp(:,1))/N_par
     aux2 = sum(rimp(:,2))/N_par
     aux3 = sum(rimp(:,3))/N_par
write(81,'(1x,1p,E15.6,3E18.10)')time0, aux1, aux2, aux3
     aux1 = sum(vimp(:,1))/N_par
     aux2 = sum(vimp(:,2))/N_par
     aux3 = sum(vimp(:,3))/N_par

write(82,'("# Tiempo(ps), Vx(AA/ps), Vy(AA/ps), Vz(AA/ps)")')
write(82,'(1x,1p,E15.6,3E18.10)')time0, aux1*pstoK, aux2*pstoK, aux3*pstoK
     aux1 = sum(aimp(:,1))/N_par
     aux2 = sum(aimp(:,2))/N_par
     aux3 = sum(aimp(:,3))/N_par
write(83,'("# Tiempo(ps), Ax(AA/ps**2), Ay(AA/ps**2), Az(AA/ps**2)")')
write(83,'(1x,1p,E15.6,3E18.10)')time0, aux1*pstoK*pstoK, aux2*pstoK*pstoK, aux3*pstoK*pstoK

call energy()             ! Calculate energies

! TEST: Treu el valor de uimp
call respar(x,y,z,nx,ny,nz,1,'uimp','den',uimp,den)

  write(6,'("Number of He4 atoms",1P,E15.6)')auxn4

  write(6,6050) etot4,etot4/n4,ekin4,elj4,ealphas,esolid,ecor4
  write(6,6060) eimpu,ekinx,ekinqx,eHeX,0.d0,etot

        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop
        write(6,'("C.M. He..:",1p,3E15.6)') xcm4,ycm4,zcm4
        aux1 = sum(rimp(:,1))/N_par
        aux2 = sum(rimp(:,2))/N_par
        aux3 = sum(rimp(:,3))/N_par
        write(6,'("C.M. X..:",1p,3E15.6)') aux1,aux2,aux3
   ipx = 1.5 + (aux1 +xmax)/hx
   ipy = 1.5 + (aux2 +ymax)/hy
   ipz = 1.5 + (aux3 +zmax)/hz
        aux1 = sum(vimp(:,1))/N_par
        aux2 = sum(vimp(:,2))/N_par
        aux3 = sum(vimp(:,3))/N_par
        write(6,'("C.M. V..:",1p,3E15.6)') aux1,aux2,aux3
        Write(6,'("Min val V..:",1p,3E15.6)')Minval(vimp(:,1)), Minval(vimp(:,2)),Minval(vimp(:,3))
        Write(6,'("Max val V..:",1p,3E15.6)')Maxval(vimp(:,1)), Maxval(vimp(:,2)),Maxval(vimp(:,3))
        Write(6,'("Min val a..:",1p,3E15.6)')Minval(aimp(:,1)), Minval(aimp(:,2)),Minval(aimp(:,3))
        Write(6,'("Max val a..:",1p,3E15.6)')Maxval(aimp(:,1)), Maxval(aimp(:,2)),Maxval(aimp(:,3))

        Open(Unit =1, file="denx-x.dat")
        do ix = 1, nx
          Write(1,'(1p,2E15.6)') x(ix),denx(ix,ipy,ipz)
        EndDo
        Close(Unit=1)

        Open(Unit =1, file="denx-y.dat")
        do iy = 1, ny
          Write(1,'(1p,2E15.6)') y(iy),denx(ipx,iy,ipz)
        EndDo
        Close(Unit=1)

        Open(Unit =1, file="denx-z.dat")
        do iz = 1, nz
          Write(1,'(1p,2E15.6)') z(iz),denx(ipx,ipy,iz)
        EndDo
        Close(Unit=1)
        If(Lbubble_radius)Then
          Open(Unit=61,File=filebubble)      
          Write(61,'("# Tiempo(ps), r(X)(AA), r(Y)(AA), r(Z)(AA)")')
          If(Lspherical_bubble)Then
            Open(Unit=62,File=file_spherical_bubble)      
            Write(62,'("# Tiempo(ps), r(X)(AA), r(Y)(AA), r(Z)(AA)")')
          EndIf        
!
!  We will compute the xmsx, ymsx & zmsx, to control the minimum value of bubble radius
!
          xmsx = 0.d0
          ymsx = 0.d0
          zmsx = 0.d0
          Do iz = 1, nz
            Do iy = 1, ny
              Do ix = 1, nx
                xmsx = xmsx + denx(ix,iy,iz)*x(ix)**2
                ymsx = ymsx + denx(ix,iy,iz)*y(iy)**2
                zmsx = zmsx + denx(ix,iy,iz)*z(iz)**2
              EndDo  
            EndDo  
          EndDo  
          xmsx = Sqrt(xmsx*dxyz)
          ymsx = Sqrt(ymsx*dxyz)
          zmsx = Sqrt(zmsx*dxyz)
          Write(6,'(" Impurity mean square radius:(X),(Y),(Z): ",1p,3E15.6)')xmsx, ymsx, zmsx
!
!  We compute the bubble radius around He*
!    
          Call Bubble_Radius(1, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,ixl,ixr,xbl,xbr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in X direction...:",1p,3E15.6)')(xbr-xbl)*0.5
             LcontrolX=.true.
          Else
             Write(6,'("The bubble radius in X direction cannot be found...:")')
             LcontrolX=.false.
          Endif        
          Call Bubble_Radius(2, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,iyl,iyr,ybl,ybr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in Y direction...:",1p,3E15.6)')(ybr-ybl)*0.5
             LcontrolY=.true.
          Else
             Write(6,'("The bubble radius in Y direction cannot be found...:")')
             LcontrolY=.false.
          Endif        
          Call Bubble_Radius(3, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,izl,izr,zbl,zbr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in Z direction...:",1p,3E15.6)')(zbr-zbl)*0.5
             LcontrolZ=.true.
          Else
             Write(6,'("The bubble radius in Z direction cannot be found...:")')
             LcontrolZ=.false.
          Endif        
          If(LcontrolX.And.LcontrolY.And.LcontrolZ)Then
            Write(61,'(1x,1p,E15.6,3E18.10)')time0, (xbr-xbl)*0.5, (ybr-ybl)*0.5, (zbr-zbl)*0.5
          Endif        
          If(Lspherical_bubble)Then
             Allocate(denr(nr)); Allocate(rr(nr))
             Call Bubble_Radius_esferic(den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,nr,rr,denr,ir,rb,Lcontrol)
             If(Lcontrol(1).And.Lcontrol(2))Then
               Write(6,'("Spherical bubble radius........:",1p,E15.6)')rb
               Write(62,'(1p,E15.6,E18.10)')time0, rb
             EndIf
          EndIf        
        Endif

  eold = etot
  call flush(6)

!-------------------------------------------------------------------------------
!---                            Iterative procedure                           --
!-------------------------------------------------------------------------------

! TIME CONSTANT !
! This time it's a cylinder.
do iz=1,nz
 do iy=1,ny
  do ix=1,nx
!    rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy)+z(iz)*z(iz))
   rt = dsqrt(x(ix)*x(ix)+y(iy)*y(iy))
!    timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((rt-tmean)/tsurf)),1.d0)
   timec(ix,iy,iz)=cmplx(Lambdah*(1.d0+tanh((z(iz)-tzmean)/tzsurf)),1.d0)
  enddo
 enddo
enddo


!plot it
! call respar(x,y,z,nx,ny,nz,1,'timec','den',timec,den)
open(unit=32,file='timec.dat')
do iz=1,nz
 do ix=1,nx
  write(32,*)x(ix),z(iz),real(timec(ix,ny/2+1,iz))
 enddo
 write(32,*)''
enddo
close(32)


lfilepv  = .false.

! write(6,7000)
call timer(t5)


!  Compute external potential
!   Surface TiO2
!if(lsurf)then
! write(*,*)'Compute surface potential'
! if(lsurf3D)then
!   do iz=1,nz
!    do iy=1,ny
!     do ix=1,nx
!          uext(ix,iy,iz) = min(3000.d0,Morse_HeTiO2_3D(x(ix),y(iy),z(iz)-Zsurf))
!!        uextimp(ix,iy,iz) = min(3000.d0,?????????(x(ix),y(iy),z(iz)-Zsurf))
!     enddo
!    enddo
!   enddo
! else
!   do iz=1,nz
!    do iy=1,ny
!     do ix=1,nx
!          uext(ix,iy,iz) = min(3000.d0,Morse_HeTiO2_1D(z(iz)-Zsurf))
!!        uextimp(ix,iy,iz) = min(3000.d0,?????????(z(iz)-Zsurf))
!     enddo
!    enddo
!   enddo
! endif
!else
! write(*,*)"Dont Compute surface potential"
!endif


!open(unit=42,file='uext.dat')
! do iz=1,nz
!     write(42,*)z(iz),uext(3,5,iz)
! enddo
!
!close(42)
!write(*,*)'Compute surface potential - done'

uext =0.d0

qr(:,:) = 0.d0
qv(:,:) = 0.d0
rimpold(:,:,:) = 0.d0
vimpold(:,:,:) = 0.d0
aimpold(:,:,:) = 0.d0


iter0=iter0+1
do iter=iter0,niter  ! <--------- Iterative procedure starts here.

    if(iter.le.(iter0+2).Or.lrk)then
      call steprk(deltat)

    else
      call steppc(deltat,errHe,errimp,errvimp)

      write(6,'("Err(He, rimp, vimp)..;",1p,3E15.6)')errHe,errimp,errvimp
    endif

    time = time0 + (iter-iter0+1)*deltatps
    call potenimp()
    call poten()
    call forceimp()
     aux1 = sum(rimp(:,1))/N_par
     aux2 = sum(rimp(:,2))/N_par
     aux3 = sum(rimp(:,3))/N_par

write(81,'(1x,1p,E15.6,3E18.10)')time, aux1, aux2, aux3
     aux1 = sum(vimp(:,1))/N_par
     aux2 = sum(vimp(:,2))/N_par
     aux3 = sum(vimp(:,3))/N_par

write(82,'(1x,1p,E15.6,3E18.10)')time, aux1*pstoK, aux2*pstoK, aux3*pstoK
     aux1 = sum(aimp(:,1))/N_par
     aux2 = sum(aimp(:,2))/N_par
     aux3 = sum(aimp(:,3))/N_par

write(83,'(1x,1p,E15.6,3E18.10)')time, aux1*pstoK*pstoK, aux2*pstoK*pstoK, aux3*pstoK*pstoK


   if(mod(iter,pener).eq.0) then          ! Compute New energy and max of density
      call energy()

      write(6,'("Iteration........:",I15)')iter
      auxn4 = sum(den)*dxyz
      write(6,'("Number of He4 atoms",1P,E15.6)')auxn4
      write(6,7010) etot4,(etot-eold),etot4/n4,ekin4,elj4,ealphas,esolid,ecor4

      write(6,7015) eimpu,ekinx,ekinqx,eHeX,0.d0,etot

      eold = etot

        call r_cm(den,n4,xcm4,ycm4,zcm4)    ! Center of mass of 4He Drop
        write(6,'("C.M. He..:",1p,3E15.6)') xcm4,ycm4,zcm4
     aux1 = sum(rimp(:,1))/N_par
     aux2 = sum(rimp(:,2))/N_par
     aux3 = sum(rimp(:,3))/N_par
   ipx = 1.5 + (aux1 +xmax)/hx
   ipy = 1.5 + (aux2 +ymax)/hy
   ipz = 1.5 + (aux3 +zmax)/hz
        write(6,'("C.M. X..:",1p,3E15.6)') aux1,aux2,aux3
        aux1 = sum(vimp(:,1))/N_par
        aux2 = sum(vimp(:,2))/N_par
        aux3 = sum(vimp(:,3))/N_par
        write(6,'("C.M. V..:",1p,3E15.6)') aux1,aux2,aux3
        Write(6,'("Min val V..:",1p,3E15.6)')Minval(vimp(:,1)), Minval(vimp(:,2)),Minval(vimp(:,3))
        Write(6,'("Max val V..:",1p,3E15.6)')Maxval(vimp(:,1)), Maxval(vimp(:,2)),Maxval(vimp(:,3))
        Write(6,'("Min val a..:",1p,3E15.6)')Minval(aimp(:,1)), Minval(aimp(:,2)),Minval(aimp(:,3))
        Write(6,'("Max val a..:",1p,3E15.6)')Maxval(aimp(:,1)), Maxval(aimp(:,2)),Maxval(aimp(:,3))

        Open(Unit =1, file="denx-x.dat")
        do ix = 1, nx
          Write(1,'(1p,2E15.6)') x(ix),denx(ix,ipy,ipz)
        EndDo
        Close(Unit=1)

        Open(Unit =1, file="denx-y.dat")
        do iy = 1, ny
          Write(1,'(1p,2E15.6)') y(iy),denx(ipx,iy,ipz)
        EndDo
        Close(Unit=1)

        Open(Unit =1, file="denx-z.dat")
        do iz = 1, nz
          Write(1,'(1p,2E15.6)') z(iz),denx(ipx,ipy,iz)
        EndDo
        Close(Unit=1)
        
        If(Lbubble_radius)Then
!
!  We will compute the xmsx, ymsx & zmsx, to control the minimum value of bubble radius
!
          xmsx = 0.d0
          ymsx = 0.d0
          zmsx = 0.d0
          Do iz = 1, nz
            Do iy = 1, ny
              Do ix = 1, nx
                xmsx = xmsx + denx(ix,iy,iz)*x(ix)**2
                ymsx = ymsx + denx(ix,iy,iz)*y(iy)**2
                zmsx = zmsx + denx(ix,iy,iz)*z(iz)**2
              EndDo  
            EndDo  
          EndDo  
          xmsx = Sqrt(xmsx*dxyz)
          ymsx = Sqrt(ymsx*dxyz)
          zmsx = Sqrt(zmsx*dxyz)
          Write(6,'(" Impurity mean square radius:(X),(Y),(Z): ",1p,3E15.6)')xmsx, ymsx, zmsx
!
!  We compute the bubble radius around He*
!   
          Call Bubble_Radius(1, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,ixl,ixr,xbl,xbr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in X direction...:",1p,3E15.6)')(xbr-xbl)*0.5
             LcontrolX=.true.
          Else
             Write(6,'("The bubble radius in X direction cannot be found...:")')
             LcontrolX=.false.
          Endif        
          Call Bubble_Radius(2, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,iyl,iyr,ybl,ybr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in Y direction...:",1p,3E15.6)')(ybr-ybl)*0.5
             LcontrolY=.true.
          Else
             Write(6,'("The bubble radius in Y direction cannot be found...:")')
             LcontrolY=.false.
          Endif        
          Call Bubble_Radius(3, den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,izl,izr,zbl,zbr,Lcontrol)
          If(Lcontrol(1).And.Lcontrol(2))Then
             Write(6,'("Bubble radius in Z direction...:",1p,3E15.6)')(zbr-zbl)*0.5
             LcontrolZ=.true.
          Else
             Write(6,'("The bubble radius in Z direction cannot be found...:")')
             LcontrolZ=.false.
          Endif        
          If(LcontrolX.And.LcontrolY.And.LcontrolZ)Then
            Write(61,'(1x,1p,E15.6,3E18.10)')time, (xbr-xbl)*0.5, (ybr-ybl)*0.5, (zbr-zbl)*0.5
          Endif        
          If(Lspherical_bubble)Then
             Call Bubble_Radius_esferic(den,x,y,z,nx,ny,nz,ipx,ipy,ipz,xmsx,ymsx,zmsx,nr,rr,denr,ir,rb,Lcontrol)
             If(Lcontrol(1).And.Lcontrol(2))Then
               Write(6,'("Spherical bubble radius........:",1p,E15.6)')rb
               Write(62,'(1p,E15.6,E18.10)')time, rb
             EndIf
          EndIf        
        Endif
 
   end if

!..............................................................................

   if(mod(iter,pcurr).eq.0) then        ! Save wavefunction for current

     ncurr = (iter-iter0+1)/pcurr + icurr
      select case (ncurr)
       case(1:9)
      write(chariter,8011)ncurr
       case(10:99)
      write(chariter,8012)ncurr
       case(100:999)
      write(chariter,8013)ncurr
      end select
      namefile='density.'//chariter//'.dat'
      namefile1='densityx.'//chariter//'.dat'
       call printoutc(time,3,namefile,namefile1,psi,nx,ny,nz,hx,hy,hz, &
                    xmax,ymax,zmax,rimp,vimp,N_par,        &
                    deltatps,iter)

   endif


   if(lfilepv) then
     select case(iter)
         case(1:9)
           write(namefile ,5010) iter
           write(namefile1,5015) iter
         case(10:99)
           write(namefile ,5020) iter
           write(namefile1,5025) iter
         case(100:999)
           write(namefile ,5030) iter
           write(namefile1,5035) iter
         case(1000:9999)
           write(namefile ,5040) iter
           write(namefile1,5045) iter
         case(10000:99999)
           write(namefile ,5050) iter
           write(namefile1,5055) iter
         case(100000:999999)
           write(namefile ,5060) iter
           write(namefile1,5065) iter
     end select


        call printoutc(time,2,namefile,namefile1,psi,nx,ny,nz,hx,hy,hz, &
                      xmax,ymax,zmax,rimp,vimp,N_par,                   &
                      deltatps,iter)


   end if
!..............................................................................

   call timer(t6)                         ! Compute use time
   t5=t6
end do
   call printoutc(time,1,filedenout,fileimpout,psi,nx,ny,nz,hx,hy,hz, &
                 xmax,ymax,zmax,rimp,vimp,N_par,         &
                 deltatps,iter)

call timer(t4)
print *,' Total  ',t4-t0

stop
999 stop 'DFT3He3d. Error in input master file. Too short'

!...............
!... Formats ...
!...............

3100 format(3x,0P,f9.4,2x,1P,E13.5)

3156 format(10E13.5)

6010 format(//,&
T10,'   ######  ####### ####### #       #     #          #####          ',/,  &
T10,'   #     # #          #    #    #  #     #  ###### #     #  #####  ',/,  &
T10,'   #     # #          #    #    #  #     #  #            #  #    # ',/,  &
T10,'   #     # #####      #    #    #  #######  #####   #####   #    # ',/,  &
T10,'   #     # #          #    ####### #     #  #            #  #    # ',/,  &
T10,'   #     # #          #         #  #     #  #      #     #  #    # ',/,  &
T10,'   ######  #          #         #  #     #  ######  #####   #####  ',//, &
T6,'Title of the run: ',A)

6011 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input  densitity file: ',A,/,&
               T6,'Output densitity file: ',A)

6111 format(//,T6,'CONTINUE a calculation:',//,&
               T6,'Input file with helium densitity    : ',A,/,&
               T6,'Input file with impurity wave func. : ',A,/,&
               T6,'Output file with helium densitity   : ',A,/,&
               T6,'Output file with impurity wave func.: ',A,/,' ')

6012 format(//,T6,'Start a new calculation:',//,&
               T6,'Output densitity file: ',A)
6013 format(//,T6,'Start a new calculation with an impurity:',//,&
               T6,'Output file for Helium density: ',A,/,        &
               T6,'Output file for the impurity wave function: ',A)
6018 format(//,T6,'Number of threads:    ',I6,/,&
               T6,'Number of iterations: ',i6)
6020 format(//,T6,'Number of particles:    ',0P,I10,/,&
               T6,'Radius of the cluster : ',F10.3,' A')
6025 format(//,T6,'Number of particles:    ',0P,I10,/,' ')
6030 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| REAL GRID         |     X-grid       Y-grid       Z-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Number of points  |',0P,T32,I4,T45,I4,T58,I4,T66,' |',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+')
6035 format(//,T6,'+-------------------+----------------------------------------+',/,&
               T6,'| MOMEMTUM GRID     |    Px-grid      Py-grid      Pz-grid   |',/,&
               T6,'+-------------------+----------------------------------------+',/,&
               T6,'| Step              |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Min value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'| Max value         |',1P,T27,E12.5,T40,E12.5,T53,E12.5,T66,' |',/,&
               T6,'+-------------------+----------------------------------------+' ,//,&
               T6,'Maximum modulus of p ',0P,F11.6,3x,F11.6,/,' ')
6037 format(//,T6,'Parameters for the funcional:',//,          &
               T8,'cp4 ...... ',1P,E13.5 ,' K \AA**6'      ,/, &
               T8,'cpp4 ..... ',1P,E13.5 ,' K \AA**9'      ,/, &
               T8,'den4c .... ',0P,F13.3 ,' \AA**{-3}'     ,/, &
               T8,'Alphas ... ',0P,F13.3 ,' K ^-1 \AA**3'  ,/, &
               T8,'L ........ ',0P,F13.2 ,' \AA'           ,/, &
               T8,'den0s .... ',0P,F13.2 ,' \AA**-3'       ,/, &
               T8,'h2o2m4 ... ',0P,F14.11,' hbar**2 / (2 m_4)' )
6138 format(' ',T8,'h2o2mx ... ',0P,F14.11,' hbar**2 / (2 m_x)' )
6038 format(//,T6,'Change of Paflov parameter allowed: ',//, &
     T18,'From     to      iter     Factor',/, &
     T19,'------  ------  ------  -----------')

6039 format(1x,0P,T17,i6,T25,I6,T33,i6,t42,f11.7)

6150 format(//,T6,'Pavlov parameter fixed for all the run to: ',F8.4)

6040 format( /,T6,'Lennard-Jones parameters:',//,&
               T10,'Core    ',A3,/,&
               T10,'h     ',F11.7,' A',/,&
               T10,'eps   ',F11.7,' K',/,&
               T10,'sigma ',F11.7,' A',/,&
               T10,'b     ',F11.3,' K A**3 '//,' ')
6050 format(//,T5,'FIRST ENERGY BALANCE: ',                    //     &
              ,T5,'TOTAL   energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Energy per particle (He) .....: ',F18.6,' K',/,    &
               T5,'Kinetic energy (He) ..........: ',F18.6,' K',/,    &
               T5,'Lennard-Jones energy (He) ....: ',F18.6,' K',/,    &
               T5,'Alpha_s term  energy (He) ....: ',F18.6,' K',/,    &
               T5,'Solid energy (He)  ...........: ',F18.6,' K',/,    &
               T5,'Correlation energy   (He) ....: ',F18.6,' K')
6060 format(1x,T5,'Impurity energy (X) ..........: ',F18.6,' K',/,    &
               T5,'Kinetic energy (X) ...........: ',F18.6,' K',/,    &
               T5,'Kinetic Q-energy (X) .........: ',F18.6,' K',/,    &
               T5,'Interaction energy (X-He) ....: ',F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ........: ',F18.6,' K',/,    &
               T5,'TOTAL ENERGY (He+X) ..........: ',F18.6,' K',/)
6065 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' A',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' A',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' A',/)
7000 format(//,1x,T2,'Iter     Mu(K)      Err(Mu)    Ttime  / Lap Time',/,&
 '--------------------------------------------------')
7001 format(//,1x,T2, &
     'Iter     Mu(K)      Err(Mu)    Autovalue(K)   err(K)   ETtime  / Lap Time',&
     /,74('-'))

7010 format(//,T5,'ITERATIVE PROCEDURE ',                                    //  &
              ,T5,'Total Energy (He).......... ',0P,F18.6,' K +- ',1P,e12.4,' K',&
             /,T5,'Energy per particle (He)... ',0P,F18.6,' K',/,                &
             /,T5,'Kinetic Energy (He)........ ',0P,F18.6,' K',                  &
             /,T5,'Lennard-Jones Energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Alpha_s term  energy (He).. ',0P,F18.6,' K',                  &
             /,T5,'Solid energy (He) ......... ',0P,F18.6,' K',                  &
             /,T5,'Correlation Energy  (He)... ',0P,F18.6,' K')
7015 format(   T5,'Impurity energy (X->He) ... ',0P,F18.6,' K',/,&
               T5,'Kinetic energy (X) ........ ',0P,F18.6,' K',/,&
               T5,'Kinetic Q-energy (X) ...... ',0P,F18.6,' K',/,&
               T5,'Interaction energy (X-He) . ',0P,F18.6,' K',/,    &
               T5,'Spin-Orbit energy (X) ..... ',0P,F18.6,' K',/,    &
               T5,'TOTAL energy (He+X) ....... ',0P,F18.6,' K',/,' ')

7016 format(1x,T5,'Impurity location  (x-axix) ..: ',F18.6,' K',/,    &
               T5,'                   (y-axis) ..: ',F18.6,' K',/,    &
               T5,'                   (z-axis) ..: ',F18.6,' K',/,' ')

7017 format(   T5,'Chemical Potential ........ ',0P,F18.6,' K +- ',1P,e12.4,'K')
7018 format(   T5,'Autovalue (impurity) ...... ',0P,F18.6,' K +- ',1P,e12.4,'K')

! 7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
!                      0P,F10.6,',',F10.6,',',F10.6,') A')
7100 format(/1x,T5,'Center of Mass of the Helium ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A')
7110 format(/1x,T5,'Center of Mass of the Helium ........(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Center of Mass of the Impurity ......(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,       &
             1x,T5,'Distances between centers of mass ...(', &
                    0P,F9.4,',',F9.4,',',F9.4,') A',/,' ')

7020 format(1x,1P,2E14.5)
7030 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F8.0,'/',F8.2)
7035 format(0P,I5,T7,F13.7,T21,1P,E9.2,T32,0P,F13.7,T46,1P,E9.2,T57,0P,F8.0,'/',F8.2)


8010 format('partial.density' ,i1)
8015 format('partial.densityx',i1)
8020 format('partial.density' ,i2)
8025 format('partial.densityx',i2)
8030 format('partial.density' ,i3)
8035 format('partial.densityx',i3)

8011 format('00',i1)
8012 format('0',i2)
8013 format(i3)

5010 format('density.',SS,i1,'.out')
5020 format('density.',SS,i2,'.out')
5030 format('density.',SS,i3,'.out')
5040 format('density.',SS,i4,'.out')
5050 format('density.',SS,i5,'.out')
5060 format('density.',SS,i6,'.out')

5015 format('densityx.',SS,i1,'.out')
5025 format('densityx.',SS,i2,'.out')
5035 format('densityx.',SS,i3,'.out')
5045 format('densityx.',SS,i4,'.out')
5055 format('densityx.',SS,i5,'.out')
5065 format('densityx.',SS,i6,'.out')
!         1         2         3         4         5         6         7         8
!|2345678901234567890123456789012345678901234567890123456789012345678901234567890

end program 
