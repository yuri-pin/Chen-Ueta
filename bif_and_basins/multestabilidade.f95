program multiestabilidade
    use func_base
    implicit none

    real(8) :: b1,b2
    b2 = 4.5d0 !final parameter a(2)
    b1 = 2.50d0 !initial parameter a(2)

    call rotina_multiperiodicidade(b1,b2)
   

contains     
subroutine rotina_multiperiodicidade(bi,bf)
  implicit none

  real(8), intent(in) :: bi,bf
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
  open(unit=1,file="csv_multiestabilidade/ida.bin",status="replace", form = "unformatted", access = "stream" )

  x = x = (/-0.7d0,0.3d0,-0.6d0/)

  do i=0,nsteps
    f_val = bi + i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,1)

    deallocate(sol)
  end do

  close(1)

  !==================== VOLTA ====================
  open(unit=2,file="csv_multiestabilidade/volta.bin",status="replace", form = "unformatted", access = "stream" )

  ! usa condição final da IDA corretamente
  ! (garante continuidade real da bifurcação)
  ! x já contém o último valor da IDA aqui

  do i=0,nsteps
    f_val = bf - i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)

    ! ATUALIZA depois de calcular
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,2)

    deallocate(sol)
  end do

  close(2)

  !==================== PERÍODO 3 ====================
  open(unit=3,file="csv_multiestabilidade/periodo3.bin",status="replace", form = "unformatted", access = "stream" )

  x = (/-0.1d0,0.6d0,1d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.5d0 + i*db
    if (f_val > bf) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,3)

    deallocate(sol)
  end do
  
  x = (/-0.1d0,0.6d0,1d0/)
  ! volta
  do i=0,nsteps
    f_val = 3.5d0 - i*db
    if (f_val < bi) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,3)

    deallocate(sol)
  end do

  close(3)
    !==================== PERÍODO 5 ====================
  open(unit=4,file="csv_multiestabilidade/periodo5.bin",status="replace", form = "unformatted", access = "stream" )

  x = (/-12.41583528d0,-11.86578167d0,39.49211319d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.62d0 + i*db
    if (f_val > bf) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,4)

    deallocate(sol)
  end do
  
  x = (/-12.41583528d0,-11.86578167d0,39.49211319d0/)
  ! volta
  do i=0,nsteps
    f_val = 3.62d0 - i*db
    if (f_val < bi) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,4)

    deallocate(sol)
  end do

  close(4)
!==================== PERÍODO 2 ====================
  open(unit=5,file="csv_multiestabilidade/periodo2.bin",status="replace", form = "unformatted", access = "stream" )

  x = (/0.7136508d0, 4.1874892d-2, -0.3248538d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.48d0 + i*db
    if (f_val > bf) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,5)

    deallocate(sol)
  end do
  
  x = (/0.7136508d0, 4.1874892d-2, -0.3248538d0/)

  ! volta
  do i=0,nsteps
    f_val = 3.48d0 - i*db
    if (f_val < bi) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,5)

    deallocate(sol)
  end do

  close(5)

  !==================== PERÍODO 2' ====================
  open(unit=6,file="csv_multiestabilidade/periodo2'.bin",status="replace", form = "unformatted", access = "stream" )

  x = (/0.71597913987670680d0 ,     -0.77463467974981670d0 ,       0.75211526849548882d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.48d0 + i*db
    if (f_val > bf) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,6)

    deallocate(sol)
  end do
  
  x = (/0.71597913987670680d0 ,     -0.77463467974981670d0 ,       0.75211526849548882d0/)

  ! volta
  do i=0,nsteps
    f_val = 3.48d0 - i*db
    if (f_val < bi) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,6)

    deallocate(sol)
  end do

  close(6)

end subroutine rotina_multiperiodicidade

end program multiestabilidade