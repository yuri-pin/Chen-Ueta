program lyap_exponents
    use func_base_lyap
    implicit none
    real(8) :: b1,b2
    character(len=100) :: nome_arquivo

    b1 = 2.5d0
    b2 = 4.5d0


    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "bifurcacao[", b1, ",", b2, "].bin"


    call calcula_exp(b1,b2, nome_arquivo)

contains
subroutine calcula_exp(bi,bf,arq_saida)
    implicit none
    real(8), intent(in) :: bi, bf
    character(len=*), intent(in) :: arq_saida
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
    
    open(unit=7,file="exponents/"//trim(arq_saida),status="replace", form = "unformatted", access = "stream" )


    x = (/-0.7d0,0.3d0,-0.6d0/)
    do i = 0,n_it
        f_val = bi + i*db
        a(2) = f_val

        call fecho_orbita_lyap(a,x,h,sol,M,exp_lyap)
        x = sol(M,:)
        print*, a,exp_lyap
        
        write(7) a, exp_lyap
    end do
    close(7)
end subroutine calcula_exp



end program lyap_exponents