module func_base
    implicit none
    public

contains

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

end module func_base