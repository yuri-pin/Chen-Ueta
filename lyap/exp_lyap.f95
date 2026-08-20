program lyap_exponents
    use func_base_lyap
    implicit none
    real(8) :: b1,b2
    character(len=100) :: nome_arquivo

    b1 = 2.5d0
    b2 = 4.5d0


    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "[", b1, ",", b2, "]"


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

    real(8) :: val_ini(5,3)
    real(8) :: point_init(5)
    character(len=200) :: final_name

    integer :: i , j , n_it, M
    real(8) :: f_val
    
    !=========================================
    val_ini(1,:) = (/-0.7d0,0.3d0,-0.6d0/)
    point_init(1) = bi

    val_ini(2,:) = (/-0.1d0,0.6d0,1d0/)
    point_init(2) = 3.5d0

    val_ini(3,:) = (/-12.41583528d0,-11.86578167d0,39.49211319d0/)
    point_init(3) = 3.62d0

    val_ini(4,:) = (/0.7136508d0, 4.1874892d-2, -0.3248538d0/)
    point_init(4) = 3.48d0
    
    val_ini(5,:) = (/0.71597913987670680d0 ,     -0.77463467974981670d0 ,       0.75211526849548882d0/)
    point_init(5) = 3.48d0

    !=========================================

 


    a = (/40d0,1.8d0,33d0/)
    h  = 1.0d-3  !time step of integrartion
    db = 1.0d-3  !step size of parameter a(2)
    n_it = int((bf-bi)/db)
    

    do j = 1, size(point_init)
        write(final_name, '(a, a, i1, a)') trim(arq_saida), "_exp", j, ".bin"

        open(unit=7+j,file="exponents/"//trim(final_name),status="replace", form = "unformatted", access = "stream" )


        x = val_ini(j,:)
        do i = 0,n_it
            f_val = point_init(j) + i*db
            if (f_val > bf) exit

            a(2) = f_val

            call fecho_orbita_lyap(a,x,h,sol,M,exp_lyap)
            x = sol(M,:)
            
            write(7+j) a, exp_lyap
        end do

        x = val_ini(j,:)
        do i = 0,n_it
            f_val = point_init(j) - i*db
            if (f_val < bi) exit

            a(2) = f_val

            call fecho_orbita_lyap(a,x,h,sol,M,exp_lyap)
            x = sol(M,:)
            
            write(7+j) a, exp_lyap
        end do
        close(7 +j)
        print*, "um arquivo finalizou"
    end do
end subroutine calcula_exp



end program lyap_exponents