module func_base_lyap
    implicit none
    public

contains

subroutine F(a,y,res)
  implicit none
  real(8), intent(in) :: a(3), y(3)
  real(8), intent(out) :: res(3)

  res(1) = a(1)*(y(2)-y(1))
  res(2) = ( (a(3)-a(1))*y(1) - y(1)*y(3) + a(3)*y(2) )
  res(3) = ( y(1)*y(2) - a(2)*y(3) )

end subroutine F
!=========================================================

subroutine Jac_F(a,y,W,res)
  implicit none
  real(8), intent(in) :: a(3),y(3),W(3,3)
  real(8), intent(out) :: res(3,3)

  real(8) :: J(3,3)

  J(1,1) = -a(1) ; J(1,2) = a(1); J(1,3) = 0.0d0 
  J(2,1) = (a(3)-a(1)) - y(3) ; J(2,2) = a(3); J(2,3) = y(1) 
  J(3,1) = y(2) ; J(3,2) = y(1); J(3,3) = a(2)

  res = matmul(J,W)

end subroutine Jac_F
!=========================================================
subroutine mapa_x(a,x,h,x_new)
  implicit none
  real(8), intent(in) :: a(3), x(3), h
  real(8), intent(out) :: x_new(3)

  real(8) :: k1(3), k2(3), k3(3), k4(3), xtemp(3)


  call F(a,x,k1)

  xtemp = x + k1*(h/2.0d0)

  call F(a,xtemp,k2)
  xtemp = x + k2*(h/2.0d0)

  call F(a,xtemp,k3)

  xtemp = x + h*k3

  call F(a,xtemp,k4)


  x_new = x + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0

end subroutine mapa_x
!=========================================================
subroutine mapa(a,x,W,h,x_new,W_new)
  implicit none
  real(8), intent(in) :: a(3), x(3), W(3,3), h
  real(8), intent(out) :: x_new(3), W_new(3,3)

  real(8) :: k1(3), k2(3), k3(3), k4(3), xtemp(3)
  real(8) :: k1_w(3,3), k2_w(3,3), k3_w(3,3), k4_w(3,3), wtemp(3,3)


  call F(a,x,k1)
  call Jac_F(a,x,W,k1_w)

  xtemp = x + k1*(h/2.0d0)
  wtemp = W + k1_w*(h/2.0d0)


  call F(a,xtemp,k2)
  call Jac_F(a,xtemp,wtemp,k2_w)

  xtemp = x + k2*(h/2.0d0)
  wtemp = W + k2_w*(h/2.0d0)

  call F(a,xtemp,k3)
  call Jac_F(a,xtemp,wtemp,k3_w)

  xtemp = x + h*k3
  wtemp = W + h*k3_w

  call F(a,xtemp,k4)
  call Jac_F(a,xtemp,wtemp,k4_w)

  x_new = x + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
  W_new = W + h*(k1_w + 2d0*k2_w + 2d0*k3_w + k4_w)/6d0

end subroutine mapa
!=========================================================
subroutine transi_x(a,x,h,N,sol)
  implicit none
  real(8), intent(in) :: a(3), x(3), h
  integer, intent(in) :: N
  real(8), intent(out) :: sol(3)

  integer :: i
  real(8) :: tran(3), temp(3)

  tran = x

  ! transiente
  do i=1,N
    call mapa_x(a,tran,h,temp)
    tran = temp
  end do

  sol = tran

end subroutine transi_x
!=========================================================
subroutine fecho_orbita_lyap(a,x,h,sol,M, exp_lyap)
  implicit none

  real(8), intent(in) :: a(3),x(3)
  real(8), intent(in) :: h
  integer, intent(out) :: M
  real(8), intent(out) :: exp_lyap(3)
  real(8), allocatable, intent(out) :: sol(:,:)


  integer :: i, N
  real(8) :: temp_x(3), tran_x(3),sol_x(3) !temporary x for the runing code // transient x for the transient of our problem // final solution of transient
  real(8) :: sum_log(3)
  real(8) :: W(3,3), tran_W(3,3), temp_W(3,3)

  W = 0.0d0
  W(1,1) = 1.0d0;  W(2,2) = 1.0d0;  W(3,3) = 1.0d0;
  sum_log = 0.0d0


  N = int(4.0d1/h)
  M = int(6.0d1/h)

  allocate(sol(M,3))

  tran_x = x !initial condition of x

  call transi_x(a,tran_x,h,N,sol_x)

  sol(1,:) = sol_x !! x after transient

  tran_W = W !initial condition of W
  
  do i = 2,M
    call mapa(a,sol(i-1,:),tran_W,h,temp_x,temp_W)
    sol(i,:) = temp_x
    call Gram_S(temp_W,sum_log)
    tran_W = temp_W
  end do
    

  exp_lyap = sum_log / (M*h)
end subroutine fecho_orbita_lyap


!=========================================================
subroutine Gram_S(W, sum_log)
  implicit none 

  real(8), intent(inout) :: W(3,3)
  real(8), intent(inout) :: sum_log(3)

  real(8) :: u(3,3), norm_1, norm_2, norm_3 !auxiliar variable to construct the orthogonal basis


  !!first vector normalized 
  u(:,1) = W(:,1)
  norm_1 = sqrt(dot_product(u(:,1),u(:,1)))
  sum_log(1) = sum_log(1) + log(norm_1)
  W(:,1) = u(:,1)/norm_1

  !!second vector normalized
  u(:,2) = W(:,2) - u(:,1)*dot_product(w(:,2), u(:,1))
  norm_2 = sqrt(dot_product(u(:,2),u(:,2)))
  sum_log(2) = sum_log(2) + log(norm_2)
  W(:,2) = u(:,2)/norm_2

  !!third vector normalized
  u(:,3) = W(:,3) - u(:,1)*dot_product(w(:,3), u(:,1)) - u(:,2)*dot_product(w(:,3), u(:,2))
  norm_3 = sqrt(dot_product(u(:,3),u(:,3)))
  sum_log(3) = sum_log(3) + log(norm_3)
  W(:,3) = u(:,3)/norm_3

end subroutine Gram_S

end module func_base_lyap