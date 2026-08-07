program biperiodicidade
    use func_base
    implicit none

    real(8) :: b1, b2
    character(len=100) :: nome_ida,nome_volta !,arg1, arg2

    ! 1. Lê os dois argumentos digitados no terminal
    !call get_command_argument(1, arg1)
    !call get_command_argument(2, arg2)

    ! 2. Converte as strings lidas para valores do tipo real(8)
    !read(arg1, *) b1
    !read(arg2, *) b2


    !Fiz isso para evitar ter que escrever as condições inicias desejadas toda vez, mas dá para alterar
    b1 = 0.0d0
    b2 = 1.0d1

    ! 3. Monta o nome do arquivo dinamicamente usando b1 e b2
    write(nome_ida, '(a, f3.1, a, f4.1, a)') &
        "ida[", b1, ",", b2, "].bin"

    write(nome_volta, '(a, f3.1, a, f4.1, a)') &
        "volta[", b1, ",", b2, "].bin"

    ! 4. Chama a rotina principal passando o nome montado
    call rotina_biperiodicidade(b1, b2, nome_ida, nome_volta)
   
    !Fiz isso para evitar ter que escrever as condições inicias desejadas toda vez, mas dá para alterar
    b1 = 2.5d0
    b2 = 4.5d0

    ! 3. Monta o nome do arquivo dinamicamente usando b1 e b2
    write(nome_ida, '(a, f3.1, a, f4.1, a)') &
        "ida[", b1, ",", b2, "].bin"

    write(nome_volta, '(a, f3.1, a, f4.1, a)') &
        "volta[", b1, ",", b2, "].bin"

    ! 4. Chama a rotina principal passando o nome montado
    call rotina_biperiodicidade(b1, b2, nome_ida, nome_volta)

contains     
subroutine rotina_biperiodicidade(bi, bf, arq_ida,arq_volta)
  implicit none

  real(8), intent(in) :: bi, bf  ! esses são respectivamente o parametro a(2) inicial e final da nossa 
  character(len=*), intent(in) :: arq_ida, arq_volta ! Recebe o nome do arquivo dinamicamente
  real(8) :: a(3), h, db, epsilon, f_val
  real(8) :: x(3)
  real(8), allocatable :: sol(:,:)
  logical :: fecho
  integer :: M
  integer :: i, nsteps

  a = (/40d0,1.8d0,33d0/)
  h = 0.001d0
  db = 0.001d0
  epsilon = 1d-5

  nsteps = int((bf - bi)/db)

  !==================== IDA ====================
  open(unit=11,file="csv_ida_volta/"//trim(arq_ida),status="replace", form = "unformatted", access = "stream" )

  x = (/-0.7d0,0.3d0,-0.6d0/)

  do i=0,nsteps
    f_val = bi + i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,11)

    deallocate(sol)
  end do

  close(11)

  !==================== VOLTA ====================
  open(unit=12,file="csv_ida_volta/"//trim(arq_volta),status="replace", form = "unformatted", access = "stream" )

  ! usa condição final da IDA corretamente
  ! (garante continuidade real da bifurcação)
  ! x já contém o último valor da IDA aqui

  do i=0,nsteps
    f_val = bf - i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)

    ! ATUALIZA depois de calcular
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,12)

    deallocate(sol)
  end do

  close(12)

end subroutine rotina_biperiodicidade

end program biperiodicidade