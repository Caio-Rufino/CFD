module ADiffPoint

implicit none

contains

Subroutine F1B2(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if (i < 3 .or. (n-i) < 1) then
  stop "Ponto inválido para receber F1B2"
end if

dfi = ((1./3.)*f(i+1) + (1./2.)*f(i) - f(i-1) + (1./6.)*f(i-2))/h;

end Subroutine F1B2

Subroutine F2B1(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if (i < 2 .or. (n-i) < 2) then
  stop "Ponto inválido para receber F2B1"
end if

dfi = ((-1./3.)*f(i-1) - (1./2.)*f(i) + f(i+1) - (1./6.)*f(i+2))/h;

end Subroutine F2B1

Subroutine B4(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if (i < 4) then
  stop "Ponto inválido para receber Backward assimétrico de 4 pontos"
end if

dfi = ((11./6.)*f(i) - 3.*f(i-1) + (3./2.)*f(i-2) - (1./3.)*f(i-3))/h;

end Subroutine B4

Subroutine F4(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if ((n - i) < 3) then
  stop "Ponto inválido para receber Forward assimétrico de 4 pontos"
end if

dfi = ((-11./6.)*f(i) + 3.*f(i+1) - (3./2.)*f(i+2) + (1./3.)*f(i+3))/h;

end Subroutine F4

Subroutine F3(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if ((n - i) < 2) then
  stop "Ponto inválido para receber Forward assimétrico de 3 pontos"
end if

dfi = (3.*f(i) - 4.*f(i+1) + f(i+2))/(2.*h);

end Subroutine F3

Subroutine B3(dfi,f,i,h,n)

integer, intent(in) :: n;
real, intent(in) :: h;
integer, intent(in) :: i;
real, intent(in) :: f(n);
real, intent(inout) :: dfi;

if ( i < 3 ) then
  stop "Ponto inválido para receber Backward assimétrico de 3 pontos"
end if

dfi = (-3.*f(i) + 4.*f(i-1) - f(i-2))/(2.*h);

end Subroutine B3

end module ADiffPoint
