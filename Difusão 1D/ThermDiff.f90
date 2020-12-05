!+---------------------------------------------------------------------------------------------------------------------------------+
!|                                               Solver para difusão térmica unidimensional                                        |
!+---------------------------------------------------------------------------------------------------------------------------------+
! Se atentar ao tamanho do cabeçalho ^, ele determina o tamanho máximo das linhas

!+---------------------------------------------------------------------------------------------------------------------------------+
!| Nescessita do módulo de diferenças finitas DiffSchemes                                                                          |
!| Compilar utilizando:                                                                                                            |
!| python -m numpy.f2py -c -llapack -lblas DiffImp.f90 Diff2Schemes.f90 DiffSchemes.f90 ADiffPoint.f90 ThermDiff.f90 -m ThermDiff  |
!+---------------------------------------------------------------------------------------------------------------------------------+

! ----------------------------------------------------------------------------
! Paredes adiabáticas:
! ----------------------------------------------------------------------------

Subroutine Adiabatica(T,T0,dx,dt,integrador,nt,nc)

use Diff2Schemes
use DiffSchemes
use ADiffPoint
use DiffImp

implicit none

!..............................................................................
! Variáveis de interface

real, intent(in) :: dx;     real, intent(in) :: dt;
integer, intent(in) :: nt;  integer, intent(in) :: nc
real, intent(in) :: T0(nc);
character, intent(in) :: integrador(5);

real, intent(out) :: T(nc,nt);

!..............................................................................
! Variáveis de trabalho

real, parameter :: alfa = 30;
real, parameter :: k = 30;
integer :: i;      integer :: j;      integer :: m;
real :: f(nc);     real :: g(nc);
real :: g2(nc);    real :: g3(nc);    real :: g4(nc);
real :: p1(nc);    real :: p2(nc);    real :: p3(nc);
real :: fonte0;    real :: fontef;

real :: matriz_C(nc,nc);  real :: matriz_B(nc,nc);  real :: Taux(nc);
real :: dTmem(nc);
character :: int;
character :: odf1;  character :: odf2;  character :: odf3;  character :: odf4;

int = integrador(1);
odf1 = integrador(2);
odf2 = integrador(3);
odf3 = integrador(4);
odf4 = integrador(5);

!..............................................................................
! Condição inicial

do i = 1,nc
T(i,1) = T0(i);
end do

do j = 2,nt

!..............................................................................
! Discretização

f = T(:,j-1)

select case (odf1)
case ('2')
call Center2_d2(g,f,dx,nc)
case ('4')
call Center4_d2(g,f,dx,nc)
case ('6')
call Thomas6_d2(g,f,dx,nc)
case ('8')
call Thomas8_d2(g,f,dx,nc)
case ('0')
call Thomas10_d2(g,f,dx,nc)
case default
  stop 'Discretização inválida'
end select

!..............................................................................
! Condições de contorno:

call F3(fonte0,f,1,dx,nc);
call B3(fontef,f,nc,dx,nc);

T(1,j) = T(1,j-1) + alfa*dt*g(1) - alfa*dt*fonte0/k;
T(nc,j) = T(nc,j-1) + alfa*dt*g(nc) + alfa*dt*fontef/k;

!..............................................................................
! Integração

select case (int)
case('A') ! Euler explícito **************************************************
do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*g(i);
end do
case('B') ! MacCormack *******************************************************

do i = 2,nc-1
p1(i) = T(i,j-1) + alfa*dt*g(i);
end do

p1(1) = T(1,j);
p1(nc) = T(nc,j);

select case (odf2)
case ('2')
call Center2_d2(g,p1,dx,nc)
case ('4')
call Center4_d2(g,p1,dx,nc)
case ('6')
call Thomas6_d2(g,p1,dx,nc)
case ('8')
call Thomas8_d2(g,p1,dx,nc)
case ('0')
call Thomas10_d2(g,p1,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = 0.5*(T(i,j-1) + p1(i) + alfa*dt*g(i));
end do

case('C') ! RK4 *************************************************************

do i = 2,nc-1
p1(i) = T(i,j-1) + 0.5*alfa*dt*g(i);
end do

p1(1) = T(1,j);
p1(nc) = T(nc,j);

select case (odf2)
case ('2')
call Center2_d2(g2,p1,dx,nc)
case ('4')
call Center4_d2(g2,p1,dx,nc)
case ('6')
call Thomas6_d2(g2,p1,dx,nc)
case ('8')
call Thomas8_d2(g2,p1,dx,nc)
case ('0')
call Thomas10_d2(g2,p1,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
p2(i) = T(i,j-1) + 0.5*alfa*dt*g2(i);
end do

p2(1) = T(1,j);
p2(nc) = T(nc,j);

select case (odf3)
case ('2')
call Center2_d2(g3,p2,dx,nc)
case ('4')
call Center4_d2(g3,p2,dx,nc)
case ('6')
call Thomas6_d2(g3,p2,dx,nc)
case ('8')
call Thomas8_d2(g3,p2,dx,nc)
case ('0')
call Thomas10_d2(g3,p2,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
p3(i) = T(i,j-1) + alfa*dt*g3(i);
end do

p3(1) = T(1,j);
p3(nc) = T(nc,j);

select case (odf4)
case ('2')
call Center2_d2(g4,p3,dx,nc)
case ('4')
call Center4_d2(g4,p3,dx,nc)
case ('6')
call Thomas6_d2(g4,p3,dx,nc)
case ('8')
call Thomas8_d2(g4,p3,dx,nc)
case ('0')
call Thomas10_d2(g4,p3,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*(g(i) + 2.*(g2(i) + g3(i)) + g4(i))/6.;
end do

case('D') ! Lax-Wendroff ****************************************************

select case (odf2)
case ('B')
call backward_d1(g2,g,dx,nc)
case ('F')
call forward_d1(g2,g,dx,nc)
case ('2')
call center2_d1(g2,g,dx,nc)
case ('4')
call center4_d1(g2,g,dx,nc)
case ('6')
call Thomas6_d1(g2,g,dx,nc)
case ('8')
call Thomas8_d1(g2,g,dx,nc)
case ('0')
call Thomas10_d1(g2,g,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*g(i) + ((alfa*dt)**2)*g2(i);
end do

case('E') ! Euler ***********************************************************

call ProxDiff(matriz_C,alfa*dt,dx,odf1,nc)

do i = 2,nc-1
T(i,j) =  0;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*T(m,j-1);
end do
end do

case('F') ! Trapezoidal *****************************************************

call ProxDiff(matriz_C,0.5*alfa*dt,dx,odf1,nc)
call AtualDiff(matriz_B,0.5*alfa*dt,dx,odf1,nc)

do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
end do

case('G') ! AM3 *************************************************************

if (j == 2) then
call ProxDiff(matriz_C,0.5*alfa*dt,dx,odf1,nc)
call AtualDiff(matriz_B,0.5*alfa*dt,dx,odf1,nc)
do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
end do
dTmem = g;
else
call ProxDiff(matriz_C,5.*alfa*dt/12.,dx,odf1,nc)
call AtualDiff(matriz_B,2.*alfa*dt/3.,dx,odf1,nc)
do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
T(i,j) = T(i,j) - alfa*dt*dTmem(i)/12.;
end do
dTmem = g;
end if

case default
stop "Integrador inexistente"
end select


end do

end Subroutine Adiabatica

! ----------------------------------------------------------------------------
! Paredes isotérmicas:
! ----------------------------------------------------------------------------

Subroutine Isotermica(T,T0,dx,dt,integrador,nt,nc)

use Diff2Schemes
use DiffSchemes
use ADiffPoint
use DiffImp

implicit none

!..............................................................................
! Variáveis de interface

real, intent(in) :: dx;     real, intent(in) :: dt;
integer, intent(in) :: nt;  integer, intent(in) :: nc
real, intent(in) :: T0(nc);
character, intent(in) :: integrador(5);

real, intent(out) :: T(nc,nt);

!..............................................................................
! Variáveis de trabalho

real, parameter :: alfa = 30;
real, parameter :: k = 30;
integer :: i;      integer :: j;      integer :: m;
real :: f(nc);     real :: g(nc);
real :: g2(nc);    real :: g3(nc);    real :: g4(nc);
real :: p1(nc);    real :: p2(nc);    real :: p3(nc);

real :: matriz_C(nc,nc);  real :: matriz_B(nc,nc);  real :: Taux(nc);
real :: dTmem(nc);
character :: int;
character :: odf1;  character :: odf2;  character :: odf3;  character :: odf4;

int = integrador(1);
odf1 = integrador(2);
odf2 = integrador(3);
odf3 = integrador(4);
odf4 = integrador(5);

!..............................................................................
! Condição inicial

do i = 1,nc
T(i,1) = T0(i);
end do

do j = 2,nt

!..............................................................................
! Discretização

g = 0.;
f = T(:,j-1)

select case (odf1)
case ('2')
call Center2_d2(g,f,dx,nc)
case ('4')
call Center4_d2(g,f,dx,nc)
case ('6')
call Thomas6_d2(g,f,dx,nc)
case ('8')
call Thomas8_d2(g,f,dx,nc)
case ('0')
call Thomas10_d2(g,f,dx,nc)
case default
  stop 'Discretização inválida'
end select

!..............................................................................
! Condições de contorno:

T(1,j) = T(1,j-1);
T(nc,j) = T(nc,j-1);

!..............................................................................
! Integração

select case (int)
case('A') ! Euler explícito **************************************************
do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*g(i);
end do
case('B') ! MacCormack *******************************************************

do i = 2,nc-1
p1(i) = T(i,j-1) + alfa*dt*g(i);
end do

p1(1) = T(1,j);
p1(nc) = T(nc,j);

select case (odf2)
case ('2')
call Center2_d2(g,p1,dx,nc)
case ('4')
call Center4_d2(g,p1,dx,nc)
case ('6')
call Thomas6_d2(g,p1,dx,nc)
case ('8')
call Thomas8_d2(g,p1,dx,nc)
case ('0')
call Thomas10_d2(g,p1,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = 0.5*(T(i,j-1) + p1(i) + alfa*dt*g(i));
end do

case('C') ! RK4 *************************************************************

do i = 2,nc-1
p1(i) = T(i,j-1) + 0.5*alfa*dt*g(i);
end do

p1(1) = T(1,j);
p1(nc) = T(nc,j);

select case (odf2)
case ('2')
call Center2_d2(g2,p1,dx,nc)
case ('4')
call Center4_d2(g2,p1,dx,nc)
case ('6')
call Thomas6_d2(g2,p1,dx,nc)
case ('8')
call Thomas8_d2(g2,p1,dx,nc)
case ('0')
call Thomas10_d2(g2,p1,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
p2(i) = T(i,j-1) + 0.5*alfa*dt*g2(i);
end do

p2(1) = T(1,j);
p2(nc) = T(nc,j);

select case (odf3)
case ('2')
call Center2_d2(g3,p2,dx,nc)
case ('4')
call Center4_d2(g3,p2,dx,nc)
case ('6')
call Thomas6_d2(g3,p2,dx,nc)
case ('8')
call Thomas8_d2(g3,p2,dx,nc)
case ('0')
call Thomas10_d2(g3,p2,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
p3(i) = T(i,j-1) + alfa*dt*g3(i);
end do

p3(1) = T(1,j);
p3(nc) = T(nc,j);

select case (odf4)
case ('2')
call Center2_d2(g4,p3,dx,nc)
case ('4')
call Center4_d2(g4,p3,dx,nc)
case ('6')
call Thomas6_d2(g4,p3,dx,nc)
case ('8')
call Thomas8_d2(g4,p3,dx,nc)
case ('0')
call Thomas10_d2(g4,p3,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*(g(i) + 2.*(g2(i) + g3(i)) + g4(i))/6.;
end do

case('D') ! Lax-Wendroff ****************************************************

select case (odf2)
case ('B')
call backward_d1(g2,g,dx,nc)
case ('F')
call forward_d1(g2,g,dx,nc)
case ('2')
call center2_d1(g2,g,dx,nc)
case ('4')
call center4_d1(g2,g,dx,nc)
case ('6')
call Thomas6_d1(g2,g,dx,nc)
case ('8')
call Thomas8_d1(g2,g,dx,nc)
case ('0')
call Thomas10_d1(g2,g,dx,nc)
case default
  stop 'Discretização inválida'
end select

do i = 2,nc-1
T(i,j) = T(i,j-1) + alfa*dt*g(i) + ((alfa*dt)**2)*g2(i);
end do

case('E') ! Euler ***********************************************************

call ProxDiff(matriz_C,alfa*dt,dx,odf1,nc)

do i = 2,nc-1
T(i,j) =  0;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*T(m,j-1);
end do
end do

case('F') ! Trapezoidal *****************************************************

call ProxDiff(matriz_C,0.5*alfa*dt,dx,odf1,nc)
call AtualDiff(matriz_B,0.5*alfa*dt,dx,odf1,nc)

do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
end do

case('G') ! AM3 *************************************************************

if (j == 2) then
call ProxDiff(matriz_C,0.5*alfa*dt,dx,odf1,nc)
call AtualDiff(matriz_B,0.5*alfa*dt,dx,odf1,nc)
do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
end do
dTmem = g;
else
call ProxDiff(matriz_C,5.*alfa*dt/12.,dx,odf1,nc)
call AtualDiff(matriz_B,2.*alfa*dt/3.,dx,odf1,nc)
do i = 1,nc
Taux(i) =  0;
do m = 1,nc
Taux(i) = Taux(i) + matriz_B(i,m)*T(m,j-1);
end do
end do

do i = 2,nc-1
T(i,j) = 0.;
do m = 1,nc
T(i,j) = T(i,j) + matriz_C(i,m)*Taux(m);
end do
T(i,j) = T(i,j) - alfa*dt*dTmem(i)/12.;
end do
dTmem = g;
end if

case default
stop "Integrador inexistente"
end select

end do

end Subroutine Isotermica
