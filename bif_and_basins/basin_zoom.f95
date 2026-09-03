program basin_zoom
    use func_base
    implicit none

    real(8) :: z,b,zoom
    real(8) :: x1,x2,y1,y2 
    integer :: n,i
    character(len=100) :: nome_arquivo

    ! z-coordenate inicial condition
    z = -0.6d0

    ! Attraction basin parameter 
    
    b = 3.195d0

    !Initial set of initial conditions
    x1 = -0.25d0
    x2 =  0.25d0
    y1 = -0.25d0
    y2 =  0.25d0

    n = 2 !times that we zoom in, at the attraction basin
    do i = 0,n
        write(nome_arquivo, '(i1, a, f5.3, a)') i, "basin[", b, "]_zoom.bin"
        zoom = 0.1d0**i ! zoom of ten times
        call rotina_grid_periodo(z,b,x1*zoom,x2*zoom,y1*zoom,y2*zoom,nome_arquivo)
        
    end do

contains
subroutine rotina_grid_periodo(z,b,x1,x2,y1,y2,arq_saida)
  use omp_lib 
  implicit none



  real(8), intent(in) :: z,b
  real(8), intent(in) :: x1,x2,y1,y2   !interval that we are working
  !x1<x2 and y1<y2 it must happen
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
  dx = (x2 - x1)/dble(n_x - 1)
  dy = (y2 - y1)/dble(n_y - 1)

  open(unit=64, file="csv_zoom_basin/"//trim(arq_saida), status="replace", form = "unformatted", access = "stream" )

  ! ATENÇÃO: O 'j' DEVE estar no private. 
  ! O 'sol' sendo private garante uma matriz por core.
  !$omp parallel do private(f_xval, f_yval, x_initial, sol, M, fecho, j) &
  !$omp shared(a, h, epsilon, dx, dy, n_x, n_y) &
  !$omp NUM_THREADS(5)
  do i = 1, n_x
     f_xval = x1 + dble(i-1)*dx
     
     do j = 1, n_y
        f_yval = y1 + dble(j-1)*dy
        x_initial = (/ f_xval, f_yval, z /)

        ! 1. Calcula a órbita
        ! Certifique-se que fecho_orbita use 'allocate' internamente
        call fecho_orbita(a, x_initial, h, epsilon, fecho, sol, M)

        if (allocated(sol)) then
           ! 2. Bloco crítico para escrita
           !$omp critical (write_file)
           call cont_perio(64,x_initial, sol, M, 1.0d-1, 5)

           !$omp end critical (write_file)
           
           deallocate(sol)
        end if
     end do
     
     
  end do
  !$omp end parallel do

  close(64)
  print *, "Varredura paralela concluída. Dados salvos no unit 61."





end subroutine rotina_grid_periodo


end program