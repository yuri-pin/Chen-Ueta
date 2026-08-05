program constant_bifurcation
  !Esse código buscará gerar os pontos máximos em x para gerar o mapa de bifurcação.
  !Além disso, esse código toma a mesma condição inicial para todo o conjunto de valores de a(2)
    use func_base
    implicit none

    real(8) :: b1, b2
    character(len=100) :: nome_arquivo!, arg1, arg2

    ! 1. Lê os dois argumentos digitados no terminal
    !call get_command_argument(1, arg1)
    !call get_command_argument(2, arg2)

    ! 2. Converte as strings lidas para valores do tipo real(8)
    !read(arg1, *) b1
    !'read(arg2, *) b2
     
    !Fiz isso para evitar ter que escrever as condições inicias desejadas toda vez, mas dá para alterar
    b1 = 0.0d0
    b2 = 1.0d1

    ! 3. Monta o nome do arquivo dinamicamente usando b1 e b2
    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "bifurcacao[", b1, ",", b2, "].csv"

    ! 4. Chama a rotina principal passando o nome montado
    call rotina_principal(b1, b2, nome_arquivo)


    b1 = 2.5d0
    b2 = 4.5d0

    ! 3. Monta o nome do arquivo dinamicamente usando b1 e b2
    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "bifurcacao[", b1, ",", b2, "].csv"

    ! 4. Chama a rotina principal passando o nome montado
    call rotina_principal(b1, b2, nome_arquivo)


contains     
subroutine rotina_principal(bi, bf, arq_saida)
  implicit none

  real(8), intent(in) :: bi, bf  ! esses são respectivamente o parametro a(2) inicial e final da nossa 
  character(len=*), intent(in) :: arq_saida ! Recebe o nome do arquivo dinamicamente
  real(8) :: a(3), h, db, epsilon, f_val
  real(8) :: x0(3)
  real(8), allocatable :: sol(:,:)
  logical :: fecho
  integer :: i, nsteps, M

  open(unit=15,file="csv_const_bif/"//trim(arq_saida),status="replace")

  a = (/40d0,1.8d0,33d0/)
  h = 0.0005d0
  db = 0.005d0


  epsilon = 1d-5
  !x0 = (/-0.1d0,0.5d0,-0.6d0/)
  x0 = (/0.23989898989898994d0, -0.19444444444444445d0 ,-0.6d0/)
  nsteps = int((bf-bi)/db)

  do i=0,nsteps
    f_val = 0.0d0 + i*db
    a(2) = f_val

    call fecho_orbita(a,x0,h,epsilon,fecho,sol,M)
    call max_orbita_x(a,fecho,sol,M,15)

    deallocate(sol)
  end do

  close(15)
end subroutine rotina_principal

end program constant_bifurcation