Module DiffSchemes

implicit none

contains

Subroutine Backward_d1(g,f,h,n)

  integer, intent(in) :: n
  real, intent(in) :: f(n)
  real, intent(in) :: h
  real, intent(inout) :: g(n)

  integer :: i;

  g(1) = (f(2) - f(1))/h;

  do i = 2,n
  g(i) = (f(i) - f(i-1))/h;
  end do

end Subroutine Backward_d1

Subroutine Forward_d1(g,f,h,n)

  integer, intent(in) :: n
  real, intent(in) :: f(n)
  real, intent(in) :: h
  real, intent(inout) :: g(n)

  integer :: i;

  do i = 1,n-1
  g(i) = (f(i+1) - f(i))/h;
  end do

  g(n) = (f(n) - f(n-1))/h;

end Subroutine Forward_d1

Subroutine Center2_d1(g,f,h,n)

integer, intent(in) :: n
real, intent(in) :: f(n)
real, intent(in) :: h
real, intent(inout) :: g(n)

integer :: i;

g(1) = (f(2) - f(1))/h;

do i = 2,n-1
g(i) = 0.5*(f(i+1) - f(i-1))/h;
end do

g(n) = (f(n) - f(n-1))/h;

end Subroutine Center2_d1

Subroutine Center4_d1(g,f,h,n)

integer, intent(in) :: n
real, intent(in) :: f(n)
real, intent(in) :: h
real, intent(inout) :: g(n)

integer :: i;

g(1) = (f(2) - f(1))/h;
g(2) = 0.5*(f(3) - f(1))/h;

do i = 3,n-2
g(i) = (-f(i+2)/12. + 2.*f(i+1)/3. - 2.*f(i-1)/3. + f(i-2)/12.)/h;
end do

g(n-1) = 0.5*(f(n) - f(n-2))/h;
g(n) = (f(n) - f(n-1))/h;

end Subroutine Center4_d1

Subroutine Thomas6_d1(g,f,h,n)

! Variáveis de entrada:

integer, intent(in) :: n
real, intent(in) :: f(n)
real, intent(in) :: h
real, intent(inout) :: g(n)

! Variáveis de trabalho

integer :: i;

real :: aux1(n);
real :: aux2(n);
real :: x(n);
real :: alfa;
real :: a;
real :: b;

!******************************************************************************

alfa = 1./3.;
a = 14./9.;
b = 1./9.;

!******************************************************************************

! Discretização:

g = 0.;
x = 0.;
aux1 = 0.;
aux2 = 0.;

x(1) = (f(2) - f(1))/h;
x(2) = 0.5*(f(3) - f(1))/h;
x(3) = (-f(5)/12. + 2.*f(4)/3. - 2.*f(2)/3. + f(1)/12.)/h;

x(n) = (f(n) - f(n-1))/h;
x(n-1) = 0.5*(f(n) - f(n-2))/h;
x(n-2) = (-f(n)/12. + 2.*f(n-1)/3. - 2.*f(n-3)/3. + f(n-4)/12.)/h;

do i = 4,n-3
x(i) = (b*(f(i+2) - f(i-2))/4. + a*(f(i+1) - f(i-1))/2.)/h
end do

do i = 1,n
if (i<4 .or. i>n-3) then
aux1(i) = 0.;
aux2(i) = x(i);
else
aux1(i) = alfa/(1  - aux1(i-1)*alfa);
aux2(i) = (x(i) - aux2(i-1)*alfa)/(1  - aux1(i-1)*alfa);
end if
end do

g(n) = aux2(n);

do i = n-1,1,-1
g(i) = aux2(i) - aux1(i)*g(i+1);
end do

end Subroutine Thomas6_d1

Subroutine Thomas8_d1(g,f,h,n)

! Variáveis de entrada:

integer, intent(in) :: n
real, intent(in) :: f(n)
real, intent(in) :: h
real, intent(inout) :: g(n)

! Variáveis de trabalho

integer :: i;

real :: aux1(n);
real :: aux2(n);
real :: x(n);
real :: alfa;
real :: a;
real :: b;
real :: c;

!******************************************************************************

alfa = 3./8.;
a = (1./6.)*(alfa + 9.);
b = (1./15.)*(32.*alfa - 9.);
c = (1./10.)*(-3.*alfa + 1.);

!******************************************************************************

! Discretização:

g = 0.;
x = 0.;
aux1 = 0.;
aux2 = 0.;

x(1) = (f(2) - f(1))/h;
x(2) = 0.5*(f(3) - f(1))/h;
x(3) = (-f(5)/12. + 2.*f(4)/3. - 2.*f(2)/3. + f(1)/12.)/h;

x(n) = (f(n) - f(n-1))/h;
x(n-1) = 0.5*(f(n) - f(n-2))/h;
x(n-2) = (-f(n)/12. + 2.*f(n-1)/3. - 2.*f(n-3)/3. + f(n-4)/12.)/h;

do i = 4,n-3
x(i) = (c*(f(i+3)-f(i-3))/6. + b*(f(i+2) - f(i-2))/4. + a*(f(i+1) - f(i-1))/2.)/h
end do

do i = 1,n
if (i<4 .or. i>n-3) then
aux1(i) = 0.;
aux2(i) = x(i);
else
aux1(i) = alfa/(1  - aux1(i-1)*alfa);
aux2(i) = (x(i) - aux2(i-1)*alfa)/(1  - aux1(i-1)*alfa);
end if
end do

g(n) = aux2(n);

do i = n-1,1,-1
g(i) = aux2(i) - aux1(i)*g(i+1);
end do

end Subroutine Thomas8_d1

Subroutine Thomas10_d1(g,f,h,n)

! *****************************************************************************
! Variáveis de entrada/saída

integer, intent(in) :: n;
real, intent(in) :: f(n);
real, intent(in) :: h;
real, intent(out) :: g(n);

! Variáveis do problema

real, parameter :: beta = 1./36.;
real, parameter :: alfa = 4./9.;
real, parameter :: a = 40./27.;
real, parameter :: b = 25./54.;
real, parameter :: c = 0;

! Variáveis de trabalho:

integer :: i;
integer :: j;
integer :: k;
real :: y(n);
real :: z(n);
real :: alfaV(n);
real :: betaV(n);
real :: gama;
real :: mu;

! *****************************************************************************
! Algoritmo de solução
! *****************************************************************************

y = 0.;

! *****************************************************************************
! Lado direito

y(1) = (f(2) - f(1))/h;
y(2) = 0.5*(f(3) - f(1))/h;
y(3) = (-f(5)/12. + 2.*f(4)/3. - 2.*f(2)/3. + f(1)/12.)/h;

y(n-2) = (-f(n)/12. + 2.*f(n-1)/3. - 2.*f(n-3)/3. + f(n-4)/12.)/h;
y(n-1) = 0.5*(f(n) - f(n-2))/h;
y(n) = (f(n) - f(n-1))/h;

do i = 4,n-3
y(i) = (c*(f(i+3) - f(i-3))/6. + b*(f(i+2) - f(i-2))/4. + a*(f(i+1) - f(i-1))/2.)/h;
end do

! *****************************************************************************
! Lado esquerdo

g(1) = y(1);
g(2) = y(2);
g(3) = y(3);
g(n-2) = y(n-2);
g(n-1) = y(n-1);
g(n) = y(n);

mu = 0.;
gama = 0.;
alfaV(2) = 0.;
betaV(2) = 0.;
alfaV(3) = 0.;
betaV(3) = 0.;
z = y;

do i = 4,n-3
gama = alfa - alfaV(i-2)*beta;
mu = 1 - betaV(i-2)*beta - alfaV(i-1)*gama;
alfaV(i) = (alfa - betaV(i-1)*gama)/mu;
betaV(i) = beta/mu;
z(i) = (y(i) - z(i-2)*beta - z(i-1)*gama)/mu;
end do

do i = n-3,4,-1
g(i) = z(i) - alfaV(i)*g(i+1) - betaV(i)*g(i+2)
end do

end Subroutine Thomas10_d1

end module DiffSchemes
