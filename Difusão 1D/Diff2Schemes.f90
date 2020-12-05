Module Diff2Schemes

implicit none

contains

Subroutine Thomas10_d2(g,f,h,n)

! *****************************************************************************
! Variáveis de entrada/saída

integer, intent(in) :: n;
real, intent(in) :: f(n);
real, intent(in) :: h;
real, intent(out) :: g(n);

! Variáveis do problema

real, parameter :: alfa = 334./899.;
real, parameter :: beta = 43./1798.;
real, parameter :: a = 1065./1798.;
real, parameter :: b = 1038./899.;
real, parameter :: c = 79./1798.;

! Variáveis de trabalho:

integer :: i;
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

y(2) = (f(3) - 2.*f(2) + f(1))/(h**2);
y(3) = (f(4) - 2.*f(3) + f(2))/(h**2);
y(1) = 2.*y(2) - y(3);

y(n-1) = (f(n) - 2.*f(n-1) + f(n-2))/(h**2);
y(n-2) = (f(n-1) - 2.*f(n-2) + f(n-3))/(h**2);
y(n) = 2.*y(n-1) - y(n-2);

do i = 4,n-3
y(i) = (a*(f(i+1) - 2.*f(i) + f(i-1)) + b*(f(i+2) - 2.*f(i) + f(i-2))/4. + c*(f(i+3) - 2.*f(i) + f(i-3))/9.)/(h**2.);
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

end  Subroutine Thomas10_d2

Subroutine Thomas8_d2(g,f,h,n)

! *****************************************************************************
! Variáveis de entrada/saída

integer, intent(in) :: n;
real, intent(in) :: f(n);
real, intent(in) :: h;
real, intent(out) :: g(n);

! Variáveis do problema

real, parameter :: alfa = 344./1179.;
real :: beta;
real :: a;
real :: b;

! Variáveis de trabalho:

integer :: i;
real :: y(n);
real :: z(n);
real :: alfaV(n);
real :: betaV(n);
real :: gama;
real :: mu;

! *****************************************************************************
! Algoritmo de solução
! *****************************************************************************

beta = (38.*alfa - 9.)/214.;
a = (696. - 1191.*alfa)/428.;
b = (2454.*alfa - 294.)/535.;

y = 0.;

! *****************************************************************************
! Lado direito

y(2) = (f(3) - 2.*f(2) + f(1))/(h**2);
y(3) = (f(4) - 2.*f(3) + f(2))/(h**2);
y(1) = 2.*y(2) - y(3);

y(n-1) = (f(n) - 2.*f(n-1) + f(n-2))/(h**2);
y(n-2) = (f(n-1) - 2.*f(n-2) + f(n-3))/(h**2);
y(n) = 2.*y(n-1) - y(n-2);

do i = 4,n-3
y(i) = (a*(f(i+1) - 2.*f(i) + f(i-1)) + b*(f(i+2) - 2.*f(i) + f(i-2))/4.)/(h**2.);
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

end  Subroutine Thomas8_d2

Subroutine Thomas4_d2(g,f,h,n)

  integer, intent(in) :: n;
  real, intent(in) :: f(n);
  real, intent(in) :: h;
  real, intent(inout) :: g(n);
  real :: alfa;
  real :: a;  real :: b;
  real :: d(n);
  real :: aux1(n);
  real :: aux2(n);
  integer :: i;

  alfa = 1./10.;
  a = (4./3.)*(1 - alfa);
  b = (1./3.)*(-1. + 10.*alfa);

  d = 0.;

  d(2) = (f(3) - 2.*f(2) + f(1))/(h**2);
  d(3) = (f(4) - 2.*f(3) + f(2))/(h**2);
  d(1) = 2.*d(2) - d(3);

  d(n-1) = (f(n) - 2.*f(n-1) + f(n-2))/(h**2);
  d(n-2) = (f(n-1) - 2.*f(n-2) + f(n-3))/(h**2);
  d(n) = 2.*d(n-1) - d(n-2);

  do i = 4,n-3
  d(i) = (a*(f(i+3) - 2.*f(i) + f(i-3))/9. + b*(f(i+2) - 2.*f(i) + f(i-2))/4.)/(h**2.)
  end do

  do i = 1,n
  if (i<4 .or. i>n-3) then
  aux1(i) = 0.;
  aux2(i) = d(i);
  else
  aux1(i) = alfa/(1  - aux1(i-1)*alfa);
  aux2(i) = (d(i) - aux2(i-1)*alfa)/(1  - aux1(i-1)*alfa);
  end if
  end do

  g(n) = aux2(n);

  do i = n-1,1,-1
  g(i) = aux2(i) - aux1(i)*g(i+1);
  end do

end Subroutine Thomas4_d2

Subroutine Center4_d2(g,f,h,n)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: f(n)
  real, intent(in) :: h
  real, intent(inout) :: g(n)

  integer :: i;

  real, parameter :: a = -5./2.;
  real, parameter :: b = 4./3.;
  real, parameter :: c = -1./12.

  do i = 3,n-2
  g(i) = (c*(f(i+2) + f(i-2)) + b*(f(i+1) + f(i-1)) + a*f(i))/(h**2);
  end do

  g(2) = (f(3) - 2.*f(2) + f(1))/(h**2);
  g(3) = (f(4) - 2.*f(3) + f(2))/(h**2);
  g(1) = 2.*g(2) - g(3);

  g(n-1) = (f(n) - 2.*f(n-1) + f(n-2))/(h**2);
  g(n-2) = (f(n-1) - 2.*f(n-2) + f(n-3))/(h**2);
  g(n) = 2*g(n-1) - g(n-2);

end Subroutine Center4_d2

Subroutine Center2_d2(g,f,h,n)

  integer, intent(in) :: n;
  real, intent(in) :: f(n);
  real, intent(in) :: h;
  real, intent(inout) :: g(n);

  integer :: i;

g = 0.;

do i = 2,n-1
g(i) = (f(i+1) - 2.*f(i) + f(i-1))/(h**2);
end do

g(1) = 2.*g(2) - g(3);
g(n) = 2.*g(n-1) - g(n-2);

end Subroutine Center2_d2

Subroutine Thomas6_d2(g,f,h,n)

  integer, intent(in) :: n;
  real, intent(in) :: f(n);
  real, intent(in) :: h;
  real, intent(inout) :: g(n);
  real :: alfa;
  real :: a;  real :: b;
  real :: d(n);
  real :: aux1(n);
  real :: aux2(n);
  integer :: i;

  alfa = 2./11.;
  a = 12./11.;
  b = 3./11.;

  d = 0.;

  d(2) = (f(3) - 2.*f(2) + f(1))/(h**2);
  d(3) = (f(4) - 2.*f(3) + f(2))/(h**2);
  d(1) = 2.*d(2) - d(3);

  d(n-1) = (f(n) - 2.*f(n-1) + f(n-2))/(h**2);
  d(n-2) = (f(n-1) - 2.*f(n-2) + f(n-3))/(h**2);
  d(n) = 2.*d(n-1) - d(n-2);

  do i = 4,n-3
  d(i) = (a*(f(i+1) - 2.*f(i) + f(i-1))/2. + b*(f(i+2) - 2.*f(i) + f(i-2))/4.)/(h**2.)
  end do

  do i = 1,n
  if (i<4 .or. i>n-3) then
  aux1(i) = 0.;
  aux2(i) = d(i);
  else
  aux1(i) = alfa/(1  - aux1(i-1)*alfa);
  aux2(i) = (d(i) - aux2(i-1)*alfa)/(1  - aux1(i-1)*alfa);
  end if
  end do

  g(n) = aux2(n);

  do i = n-1,1,-1
  g(i) = aux2(i) - aux1(i)*g(i+1);
  end do

end Subroutine Thomas6_d2

end Module Diff2Schemes
