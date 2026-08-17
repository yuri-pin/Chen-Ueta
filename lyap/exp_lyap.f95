program lyap_exponents
    use func_base_lyap
    implicit none
    real(8) :: b1,b2
    b1 = 2.50d0
    b2 = 2.502d0

    call calcula_exp(b1,b2)

contains
subroutine calcula_exp(bi,bf)
    implicit none
    real(8), intent(in) :: bi, bf
    real(8) :: db, h
    real(8) :: a(3),x(3)
    real(8), allocatable :: sol(:,:)
    real(8) :: exp_lyap(3)



    integer :: i, n_it, M
    real(8) :: f_val


    a = (/40d0,1.8d0,33d0/)
    h  = 1.0d-3  !time step of integrartion
    db = 1.0d-3  !step size of parameter a(2)
    n_it = int((bf-bi)/db)
    


    x = (/-0.7d0,0.3d0,-0.6d0/)
    do i = 0,n_it
        f_val = bi + i*db
        a(2) = f_val

        call fecho_orbita_lyap(a,x,h,sol,M,exp_lyap)
        x = sol(M,:)
        print*, exp_lyap
    end do

end subroutine calcula_exp



end program lyap_exponents