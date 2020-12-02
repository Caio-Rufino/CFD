!+---------------------------------------------------------------------------------------------------------------------------------+
!|                                               Solver para propagação de onda unidimensional                                     |
!+---------------------------------------------------------------------------------------------------------------------------------+
! Se atentar ao tamanho do cabeçalho ^, ele determina o tamanho máximo das linhas

!+---------------------------------------------------------------------------------------------------------------------------------+
!| Nescessita do módulo de diferenças finitas DiffSchemes                                                                          |
!| Compilar utilizando:                                                                                                            |
!| python -m numpy.f2py -c Diff2Schemes.f90 DiffSchemes.f90 Solvers.f90 -m Solver                                                  |
!+---------------------------------------------------------------------------------------------------------------------------------+

Subroutine RK4(P,L,v,Ploc,dt,nt,nc,odf1,odf2,odf3,odf4)

  use DiffSchemes

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
  integer, intent(in) :: odf3;
  integer, intent(in) :: odf4;

  ! Tensores de saída

  real, intent(out) :: P(nc,nt);

  ! Variáveis de trabalho

  integer :: i;
  integer :: j;
  real :: h;

  real :: f1(nc);
  real :: f2(nc);
  real :: f3(nc);
  real :: f4(nc);
  real :: g1(nc);
  real :: g2(nc);
  real :: g3(nc);
  real :: g4(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

do j = 2,nt

f1 = P(:,j-1);
g1 = 0.;  g2 = 0.;  g3 = 0.;  g4 = 0.;

select case (odf1)
case (-1)
call backward_d1(g1,f1,h,nc)
case (1)
call forward_d1(g1,f1,h,nc)
case (2)
call center2_d1(g1,f1,h,nc)
case (4)
call center4_d1(g1,f1,h,nc)
case (6)
call Thomas6_d1(g1,f1,h,nc)
case (8)
call Thomas8_d1(g1,f1,h,nc)
case (10)
call Thomas10_d1(g1,f1,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
f2(i) = f1(i) - 0.5*v*dt*g1(i);
end do

select case (odf2)
case (-1)
call backward_d1(g2,f2,h,nc)
case (1)
call forward_d1(g2,f2,h,nc)
case (2)
call center2_d1(g2,f2,h,nc)
case (4)
call center4_d1(g2,f2,h,nc)
case (6)
call Thomas6_d1(g2,f2,h,nc)
case (8)
call Thomas8_d1(g2,f2,h,nc)
case (10)
call Thomas10_d1(g2,f2,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
f3(i) = f1(i) - 0.5*v*dt*g2(i);
end do

select case (odf3)
case (-1)
call backward_d1(g3,f3,h,nc)
case (1)
call forward_d1(g3,f3,h,nc)
case (2)
call center2_d1(g3,f3,h,nc)
case (4)
call center4_d1(g3,f3,h,nc)
case (6)
call Thomas6_d1(g3,f3,h,nc)
case (8)
call Thomas8_d1(g3,f3,h,nc)
case (10)
call Thomas10_d1(g3,f3,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
f4(i) = f1(i) - v*dt*g3(i);
end do

select case (odf4)
case (-1)
call backward_d1(g4,f4,h,nc)
case (1)
call forward_d1(g4,f4,h,nc)
case (2)
call center2_d1(g4,f4,h,nc)
case (4)
call center4_d1(g4,f4,h,nc)
case (6)
call Thomas6_d1(g4,f4,h,nc)
case (8)
call Thomas8_d1(g4,f4,h,nc)
case (10)
call Thomas10_d1(g4,f4,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
P(i,j) = f1(i) - v*dt*(g1(i) + 2.*(g2(i) + g3(i)) + g4(i))/6.;
end do

end do

end Subroutine RK4

Subroutine Gazdag(P,L,v,Ploc,dt,nt,nc,odf)

  use DiffSchemes

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

  real :: f(nc);
  real :: g(nc);
  real :: gmem(nc);
  real :: gp(nc);
  real :: pred(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Primeiro passo

g = 0.
f(:) = P(:,1);

select case (odf)
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

gmem = g;
g = 0.;
f(:) = P(:,j-1);

select case (odf)
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

! Marcha no tempo:

do j = 3,nt

! Preditor

do i = 1,nc
pred(i) = P(i,j-1) - 0.5*v*dt*(3.*g(i) - gmem(i))
end do

! Corretor

gmem = g;
g = 0.;

select case (odf)
case (-1)
call backward_d1(g,pred,h,nc)
case (1)
call forward_d1(g,pred,h,nc)
case (2)
call center2_d1(g,pred,h,nc)
case (4)
call center4_d1(g,pred,h,nc)
case (6)
call Thomas6_d1(g,pred,h,nc)
case (8)
call Thomas8_d1(g,pred,h,nc)
case (10)
call Thomas10_d1(g,pred,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
P(i,j) = P(i,j-1) - 0.5*v*dt*(g(i) + gmem(i));
end do

end do

end Subroutine Gazdag

Subroutine Burstein(P,L,v,Ploc,dt,nt,nc,odf)

  use DiffSchemes

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

  real :: f(nc);
  real :: g(nc);
  real :: pred(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Marcha no tempo:

do j = 2,nt

g = 0.;
f(:) = P(:,j);

  select case (odf)
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

! Preditor

do i = 1,nc
pred(i) = P(i,j-1) - 0.5*v*dt*g(i);
end do

! Corretor

g = 0.;

  select case (odf)
  case (-1)
  call backward_d1(g,pred,h,nc)
  case (1)
  call forward_d1(g,pred,h,nc)
  case (2)
  call center2_d1(g,pred,h,nc)
  case (4)
  call center4_d1(g,pred,h,nc)
  case (6)
  call Thomas6_d1(g,pred,h,nc)
  case (8)
  call Thomas8_d1(g,pred,h,nc)
  case (10)
  call Thomas10_d1(g,pred,h,nc)
  case default
    stop 'Discretização inválida'
  end select

do i = 1,nc
P(i,j) = P(i,j-1) - v*dt*g(i);
end do

end do

end Subroutine Burstein

Subroutine ABM3(P,L,v,Ploc,dt,nt,nc,odfp,odfc)

  use DiffSchemes

  implicit none

  ! Variáveis de entrada:

  real, intent(in) :: L;
  real, intent(in) :: Ploc;
  real, intent(in) :: v;
  real, intent(in) :: dt;
  integer, intent(in) :: nt;
  integer, intent(in) :: nc;
  integer, intent(in) :: odfp;
  integer, intent(in) :: odfc;

  ! Tensores de saída

  real, intent(out) :: P(nc,nt);

  ! Variáveis de trabalho

  integer :: i;
  integer :: j;
  real :: h;

  real :: f(nc);
  real :: g(nc);
  real :: gmem(nc);
  real :: gp(nc);
  real :: pred(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Primeiro passo

g = 0.
f(:) = P(:,1);

select case (odfp)
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

! Marcha no tempo:

do j = 3,nt

! Preditor

gmem = g;
g = 0.;
f(:) = P(:,j-1);

select case (odfp)
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
pred(i) = P(i,j-1) - 0.5*v*dt*(3.*g(i) - gmem(i))
end do

! Corretor

gp = 0.;

select case (odfc)
case (-1)
call backward_d1(gp,pred,h,nc)
case (1)
call forward_d1(gp,pred,h,nc)
case (2)
call center2_d1(gp,pred,h,nc)
case (4)
call center4_d1(gp,pred,h,nc)
case (6)
call Thomas6_d1(gp,pred,h,nc)
case (8)
call Thomas8_d1(gp,pred,h,nc)
case (10)
call Thomas10_d1(gp,pred,h,nc)
case default
  stop 'Discretização inválida'
end select

do i = 1,nc
P(i,j) = P(i,j-1) - v*dt*(5.*gp(i) + 8.*g(i) - gmem(i))/12.;
end do

end do

end Subroutine ABM3

Subroutine AB3(P,L,v,Ploc,dt,nt,nc,odf)

  use DiffSchemes

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

  real :: f(nc);
  real :: g(nc);
  real :: gmem1(nc);
  real :: gmem2(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Primeiro passo

g = 0.
f(:) = P(:,1);

select case (odf)
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

! Segundo passo

gmem1 = g;
g = 0.;

f(:) = P(:,2);

select case (odf)
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
P(i,3) = P(i,2) - 0.5*v*dt*(3.*g(i) - gmem1(i));
end do

! Demais passos

do j = 4,nt

gmem2 = gmem1;
gmem1 = g;
g = 0.;
f(:) = P(:,j-1);

select case (odf)
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
P(i,j) = P(i,j-1) - v*dt*(23.*g(i) - 16.*gmem1(i) + 5.*gmem2(i))/12.;
end do

end do

end Subroutine AB3

Subroutine AB2(P,L,v,Ploc,dt,nt,nc,odf)

  use DiffSchemes

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

  real :: f(nc);
  real :: g(nc);
  real :: gmem(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Primeiro passo

g = 0.
f(:) = P(:,1);

select case (odf)
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
f(:) = P(:,j-1);

select case (odf)
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
P(i,j) = P(i,j-1) - 0.5*v*dt*(3.*g(i) - gmem(i));
end do

end do

end Subroutine AB2

Subroutine LeapFrog(P,L,v,Ploc,dt,nt,nc,odf)

  use DiffSchemes
  use Diff2Schemes

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

  real :: f(nc);
  real :: g(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Primeiro passo:

g = 0.
f(:) = P(:,1);

select case (odf)
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

! Segundo passo (utiliza o mesmo diferencial)

do i = 1,nc
P(i,3) = P(i,1) - 2*v*dt*g(i);
end do

! Demais passos:

do j = 4,nt

  g = 0.
  f(:) = P(:,j-2);

  select case (odf)
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
P(i,j) = P(i,j-2) - 2*v*dt*g(i);
end do

end do

end Subroutine LeapFrog

Subroutine LaxWendroff(P,L,v,Ploc,dt,nt,nc,odf1,odf2)

  use DiffSchemes
  use Diff2Schemes

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

  real :: f(nc);
  real :: g1(nc);
  real :: g2(nc);

  !******************************************************************************
  ! Condições iniciais e de contorno:

  h = L/nc;

  do i = 1,nc
  P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
  end do

  !******************************************************************************
  ! Loop de cálculo:

  ! Marcha no tempo

  do j = 1,nt-1

  ! Discretização de primeira

  f(:) = P(:,j);
  g1 = 0.;

  select case (odf1)
  case (-1)
  call backward_d1(g1,f,h,nc)
  case (1)
  call forward_d1(g1,f,h,nc)
  case (2)
  call center2_d1(g1,f,h,nc)
  case (4)
  call center4_d1(g1,f,h,nc)
  case (6)
  call Thomas6_d1(g1,f,h,nc)
  case (8)
  call Thomas8_d1(g1,f,h,nc)
  case (10)
  call Thomas10_d1(g1,f,h,nc)
  case default
    stop 'Discretização inválida'
  end select

! Discretização de segunda

g2 = 0.;

select case (odf2)
case (2)
call Center2_d2(g2,f,h,nc)
case (4)
call Center4_d2(g2,f,h,nc)
case (6)
call Thomas6_d2(g2,f,h,nc)
case (8)
call Thomas8_d2(g2,f,h,nc)
case (10)
call Thomas10_d2(g2,f,h,nc)
case default
  stop 'Discretização inválida'
end select

! Integrador

do i = 1,nc
P(i,j+1) = P(i,j) - v*dt*g1(i) - ((v*dt)**2)*g2(i);
end do

end do

end Subroutine LaxWendroff

Subroutine MacCormack(P,L,v,Ploc,dt,nt,nc,odfp,odfc)

  use DiffSchemes

  implicit none

  ! Variáveis de entrada:

  real, intent(in) :: L;
  real, intent(in) :: Ploc;
  real, intent(in) :: v;
  real, intent(in) :: dt;
  integer, intent(in) :: nt;
  integer, intent(in) :: nc;
  integer, intent(in) :: odfp;
  integer, intent(in) :: odfc;

  ! Tensores de saída

  real, intent(out) :: P(nc,nt);

  ! Variáveis de trabalho

  integer :: i;
  integer :: j;
  real :: h;

  real :: f(nc);
  real :: g(nc);
  real :: pred(nc);

  !******************************************************************************
  ! Condições iniciais e de contorno:

  h = L/nc;

  do i = 1,nc
  P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
  end do

  !******************************************************************************
  ! Loop de cálculo:

  ! Marcha no tempo

  do j = 1,nt-1

  ! Discretização do preditor

  f(:) = P(:,j);
  g = 0.;

  select case (odfp)
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

! Integrador preditor

do i = 1,nc;
  pred(i) = P(i,j) - v*dt*g(i);
end do

! Discretização do corretor

g = 0.;

select case (odfc)
case (-1)
call backward_d1(g,pred,h,nc)
case (1)
call forward_d1(g,pred,h,nc)
case (2)
call center2_d1(g,pred,h,nc)
case (4)
call center4_d1(g,pred,h,nc)
case (6)
call Thomas6_d1(g,pred,h,nc)
case (8)
call Thomas8_d1(g,pred,h,nc)
case (10)
call Thomas10_d1(g,pred,h,nc)
case default
  stop 'Discretização inválida'
end select

! Integrador corretor

do i = 1,nc

P(i,j+1) = 0.5*(P(i,j) + pred(i) - v*g(i)*dt);

end do

end do

end Subroutine MacCormack

Subroutine EulerExp(P,L,v,Ploc,dt,nt,nc,odf)

use DiffSchemes

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

real :: f(nc);
real :: g(nc);

!******************************************************************************
! Condições iniciais e de contorno:

h = L/nc;

do i = 1,nc
P(i,1) = 1./(1. + 25.*((i*h - L*Ploc)**2));
end do

!******************************************************************************
! Loop de cálculo:

! Marcha no tempo

do j = 1,nt-1

! Discretização:

f(:) = P(:,j);
g = 0.;

select case (odf)
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
P(i,j+1) = P(i,j) - v*g(i)*dt;
end do

end do

end subroutine EulerExp
