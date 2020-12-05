module DiffImp

implicit none

contains

Subroutine ProxDiff(matriz_C,c,dx,odf,nc)

external :: sgetri;
external :: sgetrf;

! Variáveis de troca:

integer, intent(in) :: nc;
character, intent(in) :: odf;
real, intent(in) :: c;
real, intent(in) :: dx;
real, intent(inout) :: matriz_C(nc,nc);

! Variáveis de trabalho:

integer :: i;
integer :: m;
integer :: n;

real :: matriz_A(nc,nc);  real :: matriz_B(nc,nc);

real :: work(nc);         integer :: pivo(nc);      integer :: info;

real :: A_coef_a;   real :: A_coef_b;   real :: A_coef_c;   real :: A_coef_d;
real :: B_coef_alfa;  real :: B_coef_beta;

! Algoritmo ****************************************************************

select case (odf)
case ('2')
A_coef_a = 1./(dx**2);
A_coef_b = 0;
A_coef_c = 0;

B_coef_alfa = 0;
B_coef_beta = 0;
case ('4')
A_coef_a = -5./(2.*(dx**2));
A_coef_b = 4./(3.*(dx**2));
A_coef_c = -1./(12.*(dx**2));

B_coef_alfa = 0;
B_coef_beta = 0;
case ('6')
A_coef_a = 12./(11.*(dx**2));
A_coef_b = 3./(11.*(dx**2));
A_coef_c = 0;

B_coef_alfa = 2./11.;
B_coef_beta = 0;
case ('8')
B_coef_alfa = 344./1179.;

A_coef_a = (696. - 1191.*B_coef_alfa)/(428.*(dx**2));
A_coef_b = (2454.*B_coef_alfa - 294.)/(535.*(dx**2));
A_coef_c = 0;

B_coef_beta = (38.*B_coef_alfa - 9.)/214.;
case ('0')
A_coef_a = 1065./(1789.*(dx**2));
A_coef_b = 1038./(899.*(dx**2));
A_coef_c = 79./(1798.*(dx**2));

B_coef_alfa = 334./899.;
B_coef_beta = 43./1798.;
case default
  stop 'Discretização inválida'
end select

A_coef_d = - 2.*A_coef_a - 2.*A_coef_b - 2.*A_coef_c;

matriz_A = 0.;
matriz_B = 0.;

do i = 4,nc-3
matriz_A(i,i-3) = A_coef_c;   matriz_A(i,i+3) = A_coef_c;
matriz_A(i,i-2) = A_coef_b;   matriz_A(i,i+2) = A_coef_b;
matriz_A(i,i-1) = A_coef_a;   matriz_A(i,i+1) = A_coef_a;
matriz_A(i,i) = A_coef_d;
end do

matriz_A(2,1) = 1./(dx**2);
matriz_A(2,2) = -2./(dx**2);
matriz_A(2,3) = 1./(dx**2);

matriz_A(3,2) = 1./(dx**2);
matriz_A(3,3) = -2./(dx**2);
matriz_A(3,4) = 1./(dx**2);

matriz_A(1,1) = 2./(dx**2);
matriz_A(1,2) = -5./(dx**2);
matriz_A(1,3) = 4./(dx**2);
matriz_A(1,4) = -1./(dx**2);

matriz_A(nc-1,nc-2) = 1./(dx**2);
matriz_A(nc-1,nc-1) = -2./(dx**2);
matriz_A(nc-1,nc) = 1./(dx**2);

matriz_A(nc-2,nc-3) = 1./(dx**2);
matriz_A(nc-2,nc-2) = -2./(dx**2);
matriz_A(nc-2,nc-1) = 1./(dx**2);

matriz_A(nc,nc) = 2./(dx**2);
matriz_A(nc,nc-1) = -5./(dx**2);
matriz_A(nc,nc-2) = 4./(dx**2);
matriz_A(nc,nc-3) = -1./(dx**2);

do i = 4,nc-3
matriz_B(i,i-2) = B_coef_beta;	matriz_B(i,i+2) = B_coef_beta;
matriz_B(i,i-1) = B_coef_alfa;	matriz_B(i,i+1) = B_coef_alfa;
matriz_B(i,i) = 1;
end do

matriz_B(nc,nc) = 1;
matriz_B(nc-1,nc-1) = 1;
matriz_B(nc-2,nc-2) = 1;

matriz_B(3,3) = 1;
matriz_B(2,2) = 1;
matriz_B(1,1) = 1;

if (odf /= '2' .and. odf /= '4') then
call sgetrf(nc,nc,matriz_B,nc,pivo,info);
call sgetri(nc,matriz_B,nc,pivo,work,nc,info);
end if

do i = 1,nc
do m = 1,nc
matriz_C(i,m) = 0;
do n = 1,nc
matriz_C(i,m) = matriz_C(i,m) + matriz_B(i,n)*matriz_A(n,m);
end do
if (i == m) then
matriz_C(i,m) = 1 - c*matriz_C(i,m);
else
matriz_C(i,m) = - c*matriz_C(i,m);
end if
end do
end do

call sgetrf(nc,nc,matriz_C,nc,pivo,info);
call sgetri(nc,matriz_C,nc,pivo,work,nc,info);

end Subroutine ProxDiff

Subroutine AtualDiff(matriz_C,c,dx,odf,nc)

external :: sgetri;
external :: sgetrf;

! Variáveis de troca:

integer, intent(in) :: nc;
character, intent(in) :: odf;
real, intent(in) :: c;
real, intent(in) :: dx;
real, intent(inout) :: matriz_C(nc,nc);

! Variáveis de trabalho:

integer :: i;
integer :: m;
integer :: n;

real :: matriz_A(nc,nc);  real :: matriz_B(nc,nc);

real :: work(nc);         integer :: pivo(nc);      integer :: info;

real :: A_coef_a;   real :: A_coef_b;   real :: A_coef_c;   real :: A_coef_d;
real :: B_coef_alfa;  real :: B_coef_beta;

! Algoritmo ****************************************************************

select case (odf)
case ('2')
A_coef_a = 1./(dx**2);
A_coef_b = 0;
A_coef_c = 0;

B_coef_alfa = 0;
B_coef_beta = 0;
case ('4')
A_coef_a = -5./(2.*(dx**2));
A_coef_b = 4./(3.*(dx**2));
A_coef_c = -1./(12.*(dx**2));

B_coef_alfa = 0;
B_coef_beta = 0;
case ('6')
A_coef_a = 12./(11.*(dx**2));
A_coef_b = 3./(11.*(dx**2));
A_coef_c = 0;

B_coef_alfa = 2./11.;
B_coef_beta = 0;
case ('8')
B_coef_alfa = 344./1179.;

A_coef_a = (696. - 1191.*B_coef_alfa)/(428.*(dx**2));
A_coef_b = (2454.*B_coef_alfa - 294.)/(535.*(dx**2));
A_coef_c = 0;

B_coef_beta = (38.*B_coef_alfa - 9.)/214.;
case ('0')
A_coef_a = 1065./(1789.*(dx**2));
A_coef_b = 1038./(899.*(dx**2));
A_coef_c = 79./(1798.*(dx**2));

B_coef_alfa = 334./899.;
B_coef_beta = 43./1798.;
case default
  stop 'Discretização inválida'
end select

A_coef_d = - 2.*A_coef_a - 2.*A_coef_b - 2.*A_coef_c;

matriz_A = 0.;
matriz_B = 0.;

do i = 4,nc-3
matriz_A(i,i-3) = A_coef_c;   matriz_A(i,i+3) = A_coef_c;
matriz_A(i,i-2) = A_coef_b;   matriz_A(i,i+2) = A_coef_b;
matriz_A(i,i-1) = A_coef_a;   matriz_A(i,i+1) = A_coef_a;
matriz_A(i,i) = A_coef_d;
end do

matriz_A(2,1) = 1./(dx**2);
matriz_A(2,2) = -2./(dx**2);
matriz_A(2,3) = 1./(dx**2);

matriz_A(3,2) = 1./(dx**2);
matriz_A(3,3) = -2./(dx**2);
matriz_A(3,4) = 1./(dx**2);

matriz_A(1,1) = 2./(dx**2);
matriz_A(1,2) = -5./(dx**2);
matriz_A(1,3) = 4./(dx**2);
matriz_A(1,4) = -1./(dx**2);

matriz_A(nc-1,nc-2) = 1./(dx**2);
matriz_A(nc-1,nc-1) = -2./(dx**2);
matriz_A(nc-1,nc) = 1./(dx**2);

matriz_A(nc-2,nc-3) = 1./(dx**2);
matriz_A(nc-2,nc-2) = -2./(dx**2);
matriz_A(nc-2,nc-1) = 1./(dx**2);

matriz_A(nc,nc) = 2./(dx**2);
matriz_A(nc,nc-1) = -5./(dx**2);
matriz_A(nc,nc-2) = 4./(dx**2);
matriz_A(nc,nc-3) = -1./(dx**2);

do i = 4,nc-3
matriz_B(i,i-2) = B_coef_beta;	matriz_B(i,i+2) = B_coef_beta;
matriz_B(i,i-1) = B_coef_alfa;	matriz_B(i,i+1) = B_coef_alfa;
matriz_B(i,i) = 1;
end do

matriz_B(nc,nc) = 1;
matriz_B(nc-1,nc-1) = 1;
matriz_B(nc-2,nc-2) = 1;

matriz_B(3,3) = 1;
matriz_B(2,2) = 1;
matriz_B(1,1) = 1;

if (odf /= '2' .and. odf /= '4') then
call sgetrf(nc,nc,matriz_B,nc,pivo,info);
call sgetri(nc,matriz_B,nc,pivo,work,nc,info);
end if

do i = 1,nc
do m = 1,nc
matriz_C(i,m) = 0;
do n = 1,nc
matriz_C(i,m) = matriz_C(i,m) + matriz_B(i,n)*matriz_A(n,m);
end do
if (i == m) then
matriz_C(i,m) = 1 + c*matriz_C(i,m);
else
matriz_C(i,m) = c*matriz_C(i,m);
end if
end do
end do

end Subroutine AtualDiff

end module DiffImp
