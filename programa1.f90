program main
  implicit none

  !call rotina_biperiodicidade()
  call rotina_principal()
  !call rotina_condicoes_iniciais()
  !call rotina_grid_periodo()

contains
!=========================================================
subroutine F(a,y,res)
  implicit none
  real(8), intent(in) :: a(3), y(3)
  real(8), intent(out) :: res(3)

  res(1) = a(1)*(y(2)-y(1))
  res(2) = ( (a(3)-a(1))*y(1) - y(1)*y(3) + a(3)*y(2) )
  res(3) = ( y(1)*y(2) - a(2)*y(3) )
end subroutine F
!=========================================================
subroutine mapa(a,x,h,x_new)
  implicit none
  real(8), intent(in) :: a(3), x(3), h
  real(8), intent(out) :: x_new(3)

  real(8) :: k1(3), k2(3), k3(3), k4(3), xtemp(3)

  call F(a,x,k1)

  xtemp = x + k1*(h/2.0d0)
  call F(a,xtemp,k2)

  xtemp = x + k2*(h/2.0d0)
  call F(a,xtemp,k3)

  xtemp = x + h*k3
  call F(a,xtemp,k4)

  x_new = x + h*(k1 + 2d0*k2 + 2d0*k3 + k4)/6d0
end subroutine mapa
!=========================================================
subroutine dist(x,y,d)
  implicit none
  real(8), intent(in) :: x(3), y(3)
  real(8), intent(out) :: d

  d = maxval(abs(x-y))
end subroutine dist
!=========================================================
subroutine norma(x,d)
  implicit none
  real(8), intent(in) :: x(3)
  real(8), intent(out) :: d

  d = maxval(abs(x))
end subroutine norma
!=========================================================
subroutine fecho_orbita(a,x,h,epsilon,fecho,sol,M)
  implicit none
  real(8), intent(in) :: a(3), x(3), h, epsilon
  logical, intent(out) :: fecho
  integer, intent(out) :: M
  real(8), allocatable, intent(out) :: sol(:,:)

  integer :: i, N
  real(8) :: tran(3), temp(3)

  ! NÃO vamos mais usar fecho
  fecho = .false.

  N = int(40d0/h)
  M = int(60d0/h)

  allocate(sol(M,3))

  tran = x

  ! transiente
  do i=1,N
    call mapa(a,tran,h,temp)
    tran = temp
  end do

  sol(1,:) = tran

  ! integração completa (sem testes, sem interrupção)
  do i=2,M
    call mapa(a,sol(i-1,:),h,temp)
    sol(i,:) = temp
  end do

end subroutine fecho_orbita
!=========================================================
subroutine max_orbita_x(a,fecho,sol,M,unit)
  implicit none
  real(8), intent(in) :: a(3), sol(M,3)
  logical, intent(in) :: fecho
  integer, intent(in) :: M, unit

  integer :: i
  real(8), allocatable :: an_sol(:,:)

  allocate(an_sol(M+2,3))

  an_sol(1,:) = sol(M,:)
  an_sol(2:M+1,:) = sol
  an_sol(M+2,:) = sol(1,:)

  if (fecho) then
    write(unit,*) a(1),",",a(2),",",a(3),",",0.0d0
  else
    do i=21, M-20
      if (an_sol(i-1,1) < an_sol(i,1) .and. an_sol(i,1) > an_sol(i+1,1)) then
        write(unit,*) a(1),",",a(2),",",a(3),",",an_sol(i,1)
      end if
    end do
  end if

  deallocate(an_sol)
end subroutine max_orbita_x

!=========================================================
!=========================================================
subroutine cont_perio(unit,cond_ini, sol, M, tol, periodo_max)
  implicit none
  real(8), intent(in) :: cond_ini(3)
  integer, intent(in) :: M, unit
  real(8), intent(in) :: sol(M,3)
  real(8), intent(in) :: tol
  integer, intent(in) :: periodo_max ! <-- Novo parâmetro: ex: 10, 20...
  
  real(8) :: raw_max(M), pt_fix(M), temp
  ! O array de saída terá tamanho: 2 (cond. ini) + 1 (contador) + periodo_max
  real(8) :: conta_per(3 + periodo_max)
  integer :: n_raw, n_unique, i, j
  logical :: eh_novo

  n_raw = 0
  n_unique = 0
  conta_per = 0.0d0

  ! 1. Identifica os máximos locais
  do i = 21, M - 20
    if (sol(i-1,1) < sol(i,1) .and. sol(i,1) > sol(i+1,1)) then
      n_raw = n_raw + 1
      raw_max(n_raw) = sol(i,1)
    end if
  end do

  ! 2. Filtra os valores únicos baseado na tolerância
  if (n_raw > 0) then
    n_unique = 1
    pt_fix(1) = raw_max(1)

    do i = 2, n_raw
       eh_novo = .true.
       do j = 1, n_unique
          if (abs(raw_max(i) - pt_fix(j)) <= tol) then
             eh_novo = .false.
             exit
          end if
       end do
       if (eh_novo) then
          n_unique = n_unique + 1
          pt_fix(n_unique) = raw_max(i)
       end if
    end do

    ! 3. Ordena os pontos (Bubble Sort)
    if (n_unique > 1) then
      do i = 1, n_unique - 1
         do j = i + 1, n_unique
            if (pt_fix(i) > pt_fix(j)) then
               temp = pt_fix(i); pt_fix(i) = pt_fix(j); pt_fix(j) = temp
            end if
         end do
      end do
    end if
  end if

  ! 4. Preenche o array de saída
  conta_per(1) = cond_ini(1)     ! x0
  conta_per(2) = cond_ini(2)     ! y0
  conta_per(3) = dble(n_unique)  ! Período real encontrado
  
  ! Preenche apenas até o limite de n_unique ou periodo_max
  do i = 1, min(n_unique, periodo_max)
     conta_per(3 + i) = pt_fix(i)
  end do

  ! 5. Escreve no arquivo unit 61 com estrutura CSV
  ! Escreve x0, y0 e n_unique primeiro
  write(unit, '(F14.8, A, F14.8, A, F8.0)', advance='no') conta_per(1), ',', conta_per(2), ',', conta_per(3)
  
  ! Escreve os pontos (ou zeros se não houver ponto para aquela coluna)
  do i = 1, periodo_max
     write(unit, '(A, F14.8)', advance='no') ',', conta_per(3 + i)
  end do
  
  ! Pula para a próxima linha após terminar o registro
  write(unit, *) 

end subroutine cont_perio
!=========================================================
subroutine rotina_principal()
  implicit none

  real(8) :: a(3), h, db, epsilon, f_val
  real(8) :: x0(3)
  real(8), allocatable :: sol(:,:)
  logical :: fecho
  integer :: i, nsteps, M

  open(unit=10,file="bifurcacao.csv",status="replace")

  a = (/40d0,1.8d0,33d0/)
  h = 0.0005d0
  db = 0.05d0


  epsilon = 1d-5
  !x0 = (/-0.1d0,0.5d0,-0.6d0/)
  x0 = (/0.23989898989898994d0, -0.19444444444444445d0 ,-0.6d0/)
  nsteps = int((1.0d1-0.0d0)/db)

  do i=0,nsteps
    f_val = 0.0d0 + i*db
    a(2) = f_val

    call fecho_orbita(a,x0,h,epsilon,fecho,sol,M)
    call max_orbita_x(a,fecho,sol,M,10)

    deallocate(sol)
  end do

  close(10)
end subroutine rotina_principal
!=========================================================
subroutine rotina_biperiodicidade()
  implicit none

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

  nsteps = int((10d0 - 0d0)/db)

  !==================== IDA ====================
  open(unit=21,file="ida.csv",status="replace")

  x = (/-0.7d0,0.3d0,-0.6d0/)

  do i=0,nsteps
    f_val = 0d0 + i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,21)

    deallocate(sol)
  end do

  close(21)

  !==================== VOLTA ====================
  open(unit=22,file="volta.csv",status="replace")

  ! usa condição final da IDA corretamente
  ! (garante continuidade real da bifurcação)
  ! x já contém o último valor da IDA aqui

  do i=0,nsteps
    f_val = 10d0 - i*db
    a(2) = f_val

    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)

    ! ATUALIZA depois de calcular
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,22)

    deallocate(sol)
  end do

  close(22)

  !==================== PERÍODO 3 ====================
  open(unit=23,file="periodo3.csv",status="replace")

  x = (/-0.1d0,0.6d0,1d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.5d0 + i*db
    if (f_val > 4.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,23)

    deallocate(sol)
  end do
  
  x = (/-0.1d0,0.6d0,1d0/)
  ! volta
  do i=0,nsteps
    f_val = 3.5d0 - i*db
    if (f_val < 2.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,23)

    deallocate(sol)
  end do

  close(23)
    !==================== PERÍODO 5 ====================
  open(unit=24,file="periodo5.csv",status="replace")

  x = (/-12.41583528d0,-11.86578167d0,39.49211319d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.62d0 + i*db
    if (f_val > 4.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,24)

    deallocate(sol)
  end do
  
  x = (/-12.41583528d0,-11.86578167d0,39.49211319d0/)
  ! volta
  do i=0,nsteps
    f_val = 3.62d0 - i*db
    if (f_val < 2.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,24)

    deallocate(sol)
  end do

  close(24)
!==================== PERÍODO 2 ====================
  open(unit=25,file="periodo2.csv",status="replace")

  x = (/0.7136508d0, 4.1874892d-2, -0.3248538d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.48d0 + i*db
    if (f_val > 4.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,25)

    deallocate(sol)
  end do
  
  x = (/0.7136508d0, 4.1874892d-2, -0.3248538d0/)

  ! volta
  do i=0,nsteps
    f_val = 3.48d0 - i*db
    if (f_val < 2.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,25)

    deallocate(sol)
  end do

  close(25)

  !==================== PERÍODO 2' ====================
  open(unit=26,file="periodo2'.csv",status="replace")

  x = (/0.71597913987670680d0 ,     -0.77463467974981670d0 ,       0.75211526849548882d0/)

  ! ida
  do i=0,nsteps
    f_val = 3.48d0 + i*db
    if (f_val > 4.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,26)

    deallocate(sol)
  end do
  
  x = (/0.71597913987670680d0 ,     -0.77463467974981670d0 ,       0.75211526849548882d0/)

  ! volta
  do i=0,nsteps
    f_val = 3.48d0 - i*db
    if (f_val < 2.5d0) exit

    a(2) = f_val
    call fecho_orbita(a,x,h,epsilon,fecho,sol,M)
    x = sol(M,:)

    call max_orbita_x(a,fecho,sol,M,26)

    deallocate(sol)
  end do

  close(26)

end subroutine rotina_biperiodicidade
!=========================================================
subroutine rotina_condicoes_iniciais()
  implicit none

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


  nsteps = int((4d0 - 3d0)/db)

  

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

    write(unit_id,'(A,I0,A)') "cond_", i, ".csv"

    open(unit=100+i,file="sorteio_cond/"//trim(adjustl(unit_id)),status="replace")

    x0 = x_init(:)

    !====================
    ! variação de f
    !====================
    do j=0,nsteps

      f_val = 3.48d0 + j*db
      !f_val = 3.62+ j*db 
      if (f_val > 4.0d0) exit

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
      if (f_val < 3.0d0) exit

      a(2) = f_val
      call fecho_orbita(a,x0,h,epsilon,fecho,sol,M)
      x0 = sol(M,:)

      call max_orbita_x(a,fecho,sol,M,100+i)

      deallocate(sol)

    end do

    close(100+i)

  end do

end subroutine rotina_condicoes_iniciais

!=========================================================
subroutine rotina_grid_periodo()
  use omp_lib 
  implicit none
  real(8) :: a(3), h, epsilon, f_xval, f_yval, dx, dy
  real(8) :: x_initial(3)
  real(8), allocatable :: sol(:,:)
  logical :: fecho
  integer :: M, i, j, n_x, n_y

  ! Parâmetros teste 1
  a = (/40.0d0, 3.5d0, 33.0d0/)
  h = 0.001d0
  epsilon = 1d-5
  
  n_x = 1000
  n_y = 1000
  dx = 0.5d0/dble(n_x - 1)
  dy = 0.5d0/dble(n_y - 1)

  open(unit=61, file="per1.csv", status="replace")

  ! ATENÇÃO: O 'j' DEVE estar no private. 
  ! O 'sol' sendo private garante uma matriz por core.
  !$omp parallel do private(f_xval, f_yval, x_initial, sol, M, fecho, j) &
  !$omp shared(a, h, epsilon, dx, dy, n_x, n_y)
  do i = 1, n_x
     f_xval = -0.25d0 + dble(i-1)*dx
     
     do j = 1, n_y
        f_yval = -0.25d0 + dble(j-1)*dy
        x_initial = (/ f_xval, f_yval, -0.6d0 /)

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


!_________________________________________________________________________________
!_________________________________________________________________________________
! Parâmetros teste 2
  a = (/40.0d0, 3.62d0, 33.0d0/)
  h = 0.001d0
  epsilon = 1d-5
  
  n_x = 1000
  n_y = 1000
  dx = 0.5d0/dble(n_x - 1)
  dy = 0.5d0/dble(n_y - 1)

  open(unit=62, file="per2.csv", status="replace")

  ! ATENÇÃO: O 'j' DEVE estar no private. 
  ! O 'sol' sendo private garante uma matriz por core.
  !$omp parallel do private(f_xval, f_yval, x_initial, sol, M, fecho, j) &
  !$omp shared(a, h, epsilon, dx, dy, n_x, n_y)
  do i = 1, n_x
     f_xval = -0.25d0 + dble(i-1)*dx
     
     do j = 1, n_y
        f_yval = -0.25d0 + dble(j-1)*dy
        x_initial = (/ f_xval, f_yval, -0.6d0 /)

        ! 1. Calcula a órbita
        ! Certifique-se que fecho_orbita use 'allocate' internamente
        call fecho_orbita(a, x_initial, h, epsilon, fecho, sol, M)

        if (allocated(sol)) then
           ! 2. Bloco crítico para escrita
           !$omp critical (write_file)
           call cont_perio(62,x_initial, sol, M, 1.0d-1, 5)

           !$omp end critical (write_file)
           
           deallocate(sol)
        end if
     end do
     
     
  end do
  !$omp end parallel do

  close(62)
  print *, "Varredura paralela concluída. Dados salvos no unit 62."
!_________________________________________________________________________________
!_________________________________________________________________________________
! Parâmetros teste 3
  a = (/40.0d0, 3.195d0, 33.0d0/)
  h = 0.001d0
  epsilon = 1d-5
  
  n_x = 1000
  n_y = 1000
  dx = 0.5d0/dble(n_x - 1)
  dy = 0.5d0/dble(n_y - 1)

  open(unit=63, file="per3.csv", status="replace")

  ! ATENÇÃO: O 'j' DEVE estar no private. 
  ! O 'sol' sendo private garante uma matriz por core.
  !$omp parallel do private(f_xval, f_yval, x_initial, sol, M, fecho, j) &
  !$omp shared(a, h, epsilon, dx, dy, n_x, n_y)
  do i = 1, n_x
     f_xval = -0.25d0 + dble(i-1)*dx
     
     do j = 1, n_y
        f_yval = -0.25d0 + dble(j-1)*dy
        x_initial = (/ f_xval, f_yval, -0.6d0 /)

        ! 1. Calcula a órbita
        ! Certifique-se que fecho_orbita use 'allocate' internamente
        call fecho_orbita(a, x_initial, h, epsilon, fecho, sol, M)

        if (allocated(sol)) then
           ! 2. Bloco crítico para escrita
           !$omp critical (write_file)
           call cont_perio(63,x_initial, sol, M, 1.0d-1, 5)

           !$omp end critical (write_file)
           
           deallocate(sol)
        end if
     end do
     
     
  end do
  !$omp end parallel do

  close(63)
  print *, "Varredura paralela concluída. Dados salvos no unit 63."
end subroutine rotina_grid_periodo


end program main