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


  real(8), intent(in) :: bi,bf
  integer, parameter :: n_elementos = 40
  real(8) :: a(3), h, db, epsilon, f_val
  real(8) :: x0(3)
  real(8), allocatable :: sol(:,:)
  real(8) :: x_init(3)

  logical :: fecho
  integer :: i, j, M, nsteps
  character (60) :: unit_id

  a = (/40d0, 1.8d0, 33d0/)
  h = 0.001d0
  db = 0.001d0
  epsilon = 1d-5


  nsteps = int((bf - bi)/db)

  

  !========================
  ! loop sobre condições
  !========================
  do i=1,n_elementos

    !========================
    ! gera condições iniciais
    !========================
    call random_seed()
    call random_number(x_init)
    x_init = 0.5d0*x_init - 0.25d0   ! transforma em [-1,1]

    print *, "condicoes iniciais:"
    print *, x_init

    write(unit_id,'(A,I0,A)') "cond_", i, ".bin"

    open(unit=100+i,file="csv_random/"//trim(adjustl(unit_id)),status="replace", form = "unformatted", access = "stream" )

    x0 = x_init(:)

    !====================
    ! variação de f
    !====================
    do j=0,nsteps

      f_val = 3.48d0 + j*db
      !f_val = 3.62+ j*db 
      if (f_val > bf) exit

      a(2) = f_val
      call fecho_orbita(a,x0,h,epsilon,fecho,sol,M)
      x0 = sol(M,:)

      call max_orbita_x(a,fecho,sol,M,100 + i)

      deallocate(sol)
    end do
    
    x0 = x_init(:)
    ! volta
    do j=0,nsteps
      f_val = 3.48d0 - j*db
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

end subroutine rotina_condicoes_iniciais
end program