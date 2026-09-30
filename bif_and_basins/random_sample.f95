program random
    use func_base
    implicit none

    real(8) :: b1,b2

    b1 = 3.0d0
    b2 = 4.0d0

    call rotina_condicoes_iniciais(b1,b2)





contains 
subroutine rotina_condicoes_iniciais(bi,bf)
  implicit none


  real(8), intent(in) :: bi,bf ! These are the extreme of my interval
  integer, parameter :: n_elementos = 120 ! number of times that we choose a random number
  real(8) :: a(3), h, db, epsilon, f_val, bs(3)
  real(8) :: x0(3)
  real(8), allocatable :: sol(:,:)
  real(8) :: x_init(3)

  logical :: fecho
  integer :: i, j, M, nsteps, n_b, k
  character (60) :: unit_id
  integer :: unit_txt

  a = (/40d0, 1.8d0, 33d0/) ! initial condition of parameters
  h = 0.001d0 ! time step of numerical integrations of our problem
  db = 0.001d0 ! Increment of parameter a(2)
  epsilon = 1d-5 

  

  bs =(/3.195d0, 3.48d0,3.62d0/) ! values of a(2) which we will inicializate our code
  n_b = size(bs)

  nsteps = int((bf - bi)/db)
  unit_txt = 99
  open(unit = unit_txt, file = "csv_random/all_init_cond.txt", status = "replace", action = "write")

  do k = 1, n_b
    !========================
    ! loop sobre condições
    !========================
    do i=1+(k-1)*int(n_elementos/n_b),(k)*int(n_elementos/n_b)

      !========================
      ! gera condições iniciais
      !========================
      call random_seed()
      call random_number(x_init)
      x_init = 0.5d0*x_init - 0.25d0   ! transforma em [-0.25,0.25]

      print *, "condicoes iniciais:"
      print *, x_init ! random initial  condition

      write(unit_id,'(A,I0,A)') "cond_", i, ".bin"

      write(unit_txt, '(ES16.8E2, 2X, I10, 2X, *(ES16.8E2,1X))') bs(k), i, x_init(:)

      open(unit=100+i,file="csv_random/"//trim(adjustl(unit_id)),status="replace", form = "unformatted", access = "stream" )

      x0 = x_init(:)

      !====================
      ! variação de f
      !====================
      do j=0,nsteps

        f_val = bs(k) + j*db ! forward evolution of parameter a(2)
        !f_val = 3.62+ j*db 
        if (f_val > bf) exit

        a(2) = f_val
        call fecho_orbita(a,x0,h,epsilon,fecho,sol,M) ! calculate the orbit
        x0 = sol(M,:)

        call max_orbita_x(a,fecho,sol,M,100 + i) ! take all the local maximum of the orbit

        deallocate(sol)
      end do
      
      x0 = x_init(:)
      ! volta
      do j=0,nsteps
        f_val = bs(k) - j*db ! backward evolution of parameter a(2)
        !f_val = 3.62- j*db 
        if (f_val < bi) exit

        a(2) = f_val
        call fecho_orbita(a,x0,h,epsilon,fecho,sol,M)
        x0 = sol(M,:)

        call max_orbita_x(a,fecho,sol,M,100+i)

        deallocate(sol)

      end do

      close(100+i)

    end do

  end do
  close(unit_txt)
end subroutine rotina_condicoes_iniciais
end program