!+---------------------------------------------------------------------------------------------------------------------------------+
!|                    Solver para propagação de onda unidimensional - métodos implícitos de marcha no tempo                        |
!+---------------------------------------------------------------------------------------------------------------------------------+
! Se atentar ao tamanho do cabeçalho ^, ele determina o tamanho máximo das linhas

!+---------------------------------------------------------------------------------------------------------------------------------+
!| Nescessita do módulo de diferenças finitas DiffSchemes                                                                          |
!| Compilar utilizando:                                                                                                            |
!| python -m numpy.f2py -c Diff2Schemes.f90 DiffSchemes.f90 SolversImp.f90 -m SolverImp                                            |
!+---------------------------------------------------------------------------------------------------------------------------------+

Subroutine AM3(P,L,v,Ploc,dt,nt,nc,odfi,odfe)

use DiffSchemes

implicit none

! Variáveis de entrada:

real, intent(in) :: L;
real, intent(in) :: Ploc;
real, intent(in) :: v;
real, intent(in) :: dt;
integer, intent(in) :: nt;
integer, intent(in) :: nc;
integer, intent(in) :: odfi;
integer, intent(in) :: odfe;

! Tensores de saída

real, intent(out) :: P(nc,nt);

! Variáveis de trabalho

integer :: i;
integer :: j;
real :: h;

real :: g(nc);
real :: f(nc);
real :: gmem(nc);
real :: ptil(nc)

real :: a;  real :: b;  real :: c;

real :: aux1(nc)
real :: aux2(nc)

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

select case (odfi)
case (-1)
a = - 5.*v*dt/(12.*h);
b = 1. + 5.*v*dt/(12.*h);
c = 0.;
case (1)
a = 0.;
b = 1. - 5.*v*dt/(12.*h);
c = 5.*v*dt/(12.*h);
case (2)
a = -5.*v*dt/(24.*h);
b = 1.;
c = 5.*v*dt/(24.*h);
case default
  stop "Discretização inválida"
end select

! Primeiro passo
g = 0.
f = P(:,1);

select case (odfe)
case (-1)
call backward_d1(g,f,h,nc)
case (1)
call forward_d1(g,f,h,nc)
case (2)
call center2_d1(g,f,h,nc)
case (4)
call center4_d1(g,f,h,nc)
case (6)
call Thomas6_d1(g,f,h,nc)
case (8)
call Thomas8_d1(g,f,h,nc)
case (10)
call Thomas10_d1(g,f,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
P(i,2) = P(i,1) - v*dt*g(i);
end do

do j = 3,nt

gmem = g;
g = 0.;
f = P(:,j-1);

select case (odfe)
case (-1)
call backward_d1(g,f,h,nc)
case (1)
call forward_d1(g,f,h,nc)
case (2)
call center2_d1(g,f,h,nc)
case (4)
call center4_d1(g,f,h,nc)
case (6)
call Thomas6_d1(g,f,h,nc)
case (8)
call Thomas8_d1(g,f,h,nc)
case (10)
call Thomas10_d1(g,f,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
ptil(i) = P(i,j-1) - 2.*dt*v*g(i)/3. + v*dt*gmem(i)/12.;
end do

aux1 = 0.;
aux2 = 0.;

aux1(1) = (5.*v*dt/(12.*h))/(1. - 5.*v*dt/(12.*h));
aux2(1) = ptil(1)/(1. - 5.*v*dt/(12.*h));

do i = 2,nc-1
aux1(i) = c/(b  - a*aux1(i-1));
aux2(i) = (ptil(i) - a*aux2(i-1))/(b  - a*aux1(i-1));
end do

aux2(nc) = (ptil(nc) + 5.*v*dt*aux2(i-1)/(12.*h))/(1. + 5.*v*dt/(12.*h) + 5.*v*dt*aux1(i-1)/(12.*h));

P(nc,j) = aux2(nc);

do i = nc-1,1,-1
P(i,j) = aux2(i) - aux1(i)*ptil(i+1);
end do

end do

end Subroutine AM3

Subroutine Trapezoidal(P,L,v,Ploc,dt,nt,nc,odf1,odf2)

implicit none

! Variáveis de entrada:

real, intent(in) :: L;
real, intent(in) :: Ploc;
real, intent(in) :: v;
real, intent(in) :: dt;
integer, intent(in) :: nt;
integer, intent(in) :: nc;
integer, intent(in) :: odf1;
integer, intent(in) :: odf2;

! Tensores de saída

real, intent(out) :: P(nc,nt);

! Variáveis de trabalho

integer :: i;
integer :: j;
real :: h;

real :: a1;  real :: b1;  real :: c1;
real :: a2;  real :: b2;  real :: c2;

real :: aux1(nc)
real :: aux2(nc)
real :: ptil(nc)

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

select case (odf1)
case (-1)
a1 = 0.5*v*dt/h;
b1 = 1. - 0.5*v*dt/h;
c1 = 0.;
case (1)
a1 = 0.;
b1 = 1. + 0.5*v*dt/h;
c1 = - 0.5*v*dt/h;
case (2)
a1 = v*dt/(4.*h);
b1 = 1.;
c1 = - v*dt/(4.*h);
case default
  stop "Discretização inválida"
end select

select case (odf2)
case (-1)
a2 = - 0.5*v*dt/h;
b2 = 1. + 0.5*v*dt/h;
c2 = 0.;
case (1)
a2 = 0.;
b2 = 1. - 0.5*v*dt/h;
c2 = 0.5*v*dt/h;
case (2)
a2 = -v*dt/(4.*h);
b2 = 1.;
c2 = v*dt/(4.*h);
case default
  stop "Discretização inválida"
end select

do j = 2,nt

ptil = 0.;

ptil(1) = (1 + 0.5*v*dt/h)*P(1,j-1) - 0.5*v*dt*P(2,j-1)/h;

do i = 2,nc-1
ptil(i) = a1*P(i-1,j-1) + b1*P(i,j-1) + c1*P(i+1,j-1);
end do

ptil(nc) = (1 - 0.5*v*dt/h)*P(nc,j-1) + 0.5*v*dt*P(nc-1,j-1)/h;

aux1 = 0.;
aux2 = 0.;

aux1(1) = (0.5*v*dt/h)/(1. - 0.5*v*dt/h);
aux2(1) = ptil(1)/(1. - 0.5*v*dt/h);

do i = 2,nc-1
aux1(i) = c2/(b2  - a2*aux1(i-1));
aux2(i) = (ptil(i) - a2*aux2(i-1))/(b2  - a2*aux1(i-1));
end do

aux2(nc) = (ptil(nc) + 0.5*v*dt*aux2(i-1)/h)/(1. + 0.5*v*dt/h + 0.5*v*dt*aux1(i-1)/h);
P(nc,j) = aux2(nc);

do i = nc-1,1,-1
P(i,j) = aux2(i) - aux1(i)*P(i+1,j);
end do

end do

end Subroutine Trapezoidal

Subroutine Eulerimp(P,L,v,Ploc,dt,nt,nc,odf)

implicit none

! Variáveis de entrada:

real, intent(in) :: L;
real, intent(in) :: Ploc;
real, intent(in) :: v;
real, intent(in) :: dt;
integer, intent(in) :: nt;
integer, intent(in) :: nc;
integer, intent(in) :: odf;

! Tensores de saída

real, intent(out) :: P(nc,nt);

! Variáveis de trabalho

integer :: i;
integer :: j;
real :: h;

real :: a;  real :: b;  real :: c;

real :: aux1(nc)
real :: aux2(nc)

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

select case (odf)
case (-1)
a = -v*dt/h;
b = 1. + v*dt/h;
c = 0.;
case (1)
a = 0.;
b = 1. - v*dt/h;
c = v*dt/h;
case (2)
a = -v*dt/(2.*h);
b = 1.;
c = v*dt/(2.*h);
case default
  stop "Discretização inválida"
end select

do j = 2,nt

aux1 = 0.;
aux2 = 0.;

aux1(1) = (v*dt/h)/(1. - v*dt/h);
aux2(1) = P(1,j-1)/(1. - v*dt/h);

do i = 2,nc-1
aux1(i) = c/(b  - a*aux1(i-1));
aux2(i) = (P(i,j-1) - a*aux2(i-1))/(b  - a*aux1(i-1));
end do

aux2(nc) = (P(nc,j-1) + v*dt*aux2(i-1)/h)/(1. + v*dt/h + v*dt*aux1(i-1)/h);

P(nc,j) = aux2(nc);

do i = nc-1,1,-1
P(i,j) = aux2(i) - aux1(i)*P(i+1,j-1);
end do

end do

end Subroutine Eulerimp
