program attaction_basin
    use func_base
    implicit none

    real(8) :: z,b
    character(len=100) :: nome_arquivo

    z = -0.6d0

    ! first attraction basin
    b = 3.5d0

    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "basin[", b, "].bin"

    call rotina_grid_periodo(z,b, nome_arquivo)


    ! second attraction basin
    b = 3.62d0

    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "basin[", b, "].bin"

    call rotina_grid_periodo(z,b,nome_arquivo)

    ! third attraction basin
    b = 3.195d0

    write(nome_arquivo, '(a, f3.1, a, f4.1, a)') &
        "basin[", b, "].bin"

    call rotina_grid_periodo(z,b,nome_arquivo)


contains
subroutine rotina_grid_periodo(z,b,arq_saida)
  use omp_lib 
  implicit none



  real(8), intent(in) :: z,b
  character(len=*), intent(in) :: arq_saida ! Recebe o nome do arquivo dinamicamente

  real(8) :: a(3), h, epsilon, f_xval, f_yval, dx, dy
  real(8) :: x_initial(3)
  real(8), allocatable :: sol(:,:)
  logical :: fecho
  integer :: M, i, j, n_x, n_y


  ! Parâmetros teste 1
  a = (/40.0d0, b, 33.0d0/)
  h = 0.001d0
  epsilon = 1d-5
  
  n_x = 1000
  n_y = 1000
  dx = 0.5d0/dble(n_x - 1)
  dy = 0.5d0/dble(n_y - 1)

  open(unit=61, file="csv_attraction_basin/"//trim(arq_saida), status="replace", form = "unformatted", access = "stream" )

  ! ATENÇÃO: O 'j' DEVE estar no private. 
  ! O 'sol' sendo private garante uma matriz por core.
  !$omp parallel do private(f_xval, f_yval, x_initial, sol, M, fecho, j) &
  !$omp shared(a, h, epsilon, dx, dy, n_x, n_y) &
  !$omp NUM_THREADS(5)
  do i = 1, n_x
     f_xval = -0.25d0 + dble(i-1)*dx
     
     do j = 1, n_y
        f_yval = -0.25d0 + dble(j-1)*dy
        x_initial = (/ f_xval, f_yval, z /)

        ! 1. Calcula a órbita
        ! Certifique-se que fecho_orbita use 'allocate' internamente
        call fecho_orbita(a, x_initial, h, epsilon, fecho, sol, M)

        if (allocated(sol)) then
           ! 2. Bloco crítico para escrita
           !$omp critical (write_file)
           call cont_perio(61,x_initial, sol, M, 1.0d-1, 5)

           !$omp end critical (write_file)
           
           deallocate(sol)
        end if
     end do
     
     
  end do
  !$omp end parallel do

  close(61)
  print *, "Varredura paralela concluída. Dados salvos no unit 61."





end subroutine rotina_grid_periodo


end program