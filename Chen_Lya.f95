!
!Determina o espectro de Lyapunov para sistemas a tempo contínuo, 
!conforme metodo de Wolf et al:
!A. Wolf et al. "Determining Lyapunov exponents from a time series." Physica D, 16.3 (1985).
!Aqui consta uma versão atualizada do algoritmo de Wolf.  
!Este codigo varre um eixo de parametros num intervalo definido.
!
PROGRAM LyapunovEspectro_EixoParametro_2023MAR_V02 !Versao didatica
	USE, INTRINSIC :: omp_lib
	IMPLICIT NONE
	CHARACTER(LEN=30) :: pastaSaida = 'dadosLyapunov'
	CHARACTER(LEN=30) :: arquivoSaida = 'espectroSigma'
	INTEGER, PARAMETER :: nucleos = 16 !Para paralelizacao
	INTEGER, PARAMETER :: numVariaveis = 3, numParametros = 3 !Quantidade de variaveis e parametros no sistema
	INTEGER, PARAMETER :: transiente = 2.d5, pontos = 1.d6  !Passos de integraao a descartar; Comprimento da serie temporal
	INTEGER, PARAMETER :: grade = 16001 !Discretizacao do eixo referente ao parametro varrido
	DOUBLE PRECISION, PARAMETER :: passoRK4 = 1.d-3 !Passo de Integracao
	DOUBLE PRECISION, PARAMETER :: rhoFixo = 3.3d1, betaFixo = 4.0d1!Parametros fixos
	DOUBLE PRECISION, PARAMETER :: sigmaMinimo = 2.5d0, sigmaMaximo = 4.5d0 !Parametro varrido
	DOUBLE PRECISION, PARAMETER, DIMENSION(numVariaveis) :: condicaoInicial = (/1d-1,2d-1,3d-1/)!(/-12.41583528d0,-11.86578167d0,39.49211319d0/)
	DOUBLE PRECISION, DIMENSION(grade) :: sigmaGrade !Especificar parametro varrido
	DOUBLE PRECISION, DIMENSION(numVariaveis,grade) :: espectroLyapunov !Expoentes de Lyapunov para cada variavel
	DOUBLE PRECISION, DIMENSION(numVariaveis,numVariaveis) :: identidade !Matriz identidade: condicao inicial do espaco tangente.
	INTEGER :: ii !Contador global enprega letra duplicada
!	
	CALL SA_DiscretizaParametro(grade, sigmaMinimo, sigmaMaximo, sigmaGrade) !Especificar parametro varrido
!	
	identidade = 0.d0 !Define matriz identidade
	FORALL (ii = 1:numVariaveis) identidade(ii,ii) = 1.d0
!
	CALL OMP_SET_NUM_THREADS(nucleos)
!
!	Bloco de loop paralelo com OpenMP
	!$OMP PARALLEL DO
	DO ii = 1, grade
		!Informar parametros na ordem correta
		CALL SE_EvoluiSistema((/sigmaGrade(ii), rhoFixo, betaFixo/),&
			condicaoInicial, espectroLyapunov(:,ii))
!		
!---FLAG MISTURADO (OMP PARALELO)
		WRITE (*, *) ii, grade
	END DO
	!$OMP END PARALLEL DO
!
	CALL SYSTEM ('mkdir -p '//TRIM(pastaSaida)) !Abre pasta
!	
	CALL SD_SalvaEspectro(sigmaGrade) !Especificar o parametro varrido
!	
CONTAINS
!
!-SOMENTE A SUBROTINA << SE_Sistema >> NECESSITA SER MODIFICADA (NELA EH DECLARADO O SISTEMA DE EDOs)-------------------------------
!
!-Subrotinas de Evolucao: SE_..-----------------------------------------------------------------------------------------------------
	SUBROUTINE SE_Sistema(pLocal, x, dx, xx, dxx)!Declarar aqui o sistema de EDOs autonomo e o linearizado----------------------------
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numParametros) :: pLocal !Parametros
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numVariaveis) :: x !Estado atual do sistema
		DOUBLE PRECISION, INTENT(OUT), DIMENSION(numVariaveis) :: dx !Deriva no tempo (dx/dt) de cada variavel
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numVariaveis,numVariaveis) :: xx !Matriz sistema linearizado
		DOUBLE PRECISION, INTENT(OUT), DIMENSION(numVariaveis,numVariaveis) :: dxx !Derivada sistema linearizado
		DOUBLE PRECISION :: sigma, rho, beta !MANTER NOMES DIFERENTES DOS USADOS NA DECLARACAO GLOBAL (CABECALHO)
!
!---Atribuicao dos parametros conforme aparecem nas equacoes: somente para conforto no entendimento.
!---Esta etapa nao eh necessaria e diminui a ineficiencia. 
!---Mais eficiente trabalhar com pLocal(i) nas EDOs (OU SOMENTE p(i)).
!
		sigma = pLocal(1) !Deve estar na mesma ordem da entrada em << SG_DeterminaMaximoLocal >>
		rho = pLocal(2)
		beta = pLocal(3)
!		
!---Sistema de EDOs autonomo:
!
		dx(1) = beta*(x(2) - x(1))
		dx(2) = (rho - beta)*x(1) + rho*x(2) - x(1)*x(3)
		dx(3) = x(1)*x(2) - sigma*x(3)
!		
!---Sistema linearizado:
!
		dxx(1,:) = beta*(xx(2,:) - xx(1,:))
		dxx(2,:) = (rho - beta)*xx(1,:) + rho*xx(2,:) - xx(1,:)*x(3) - x(1)*xx(3,:)
		dxx(3,:) = xx(1,:)*x(2) + x(1)*xx(2,:) - sigma*xx(3,:)
	END SUBROUTINE SE_Sistema
!	
!-AS SUBROTINAS ABAIXO NAO PRECISAM SER MODIFICADOS (PARA O USO GERAL DESTE CODIGO)-------------------------------------------------
!
	SUBROUTINE SE_EvoluiSistema(pLocal, x0, espectroLocal)!Evolui o sistema com descarte do transiente--------------------------------
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numParametros) :: pLocal !Parametros
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numVariaveis) :: x0 !Condicao Inicial
		DOUBLE PRECISION, INTENT(OUT), DIMENSION(numVariaveis) :: espectroLocal !Expoentes de Lyapunov
		DOUBLE PRECISION, DIMENSION(numVariaveis) :: x, norma !Estado do sistema; Norma vetor coluna Gram-Schmidt
		DOUBLE PRECISION, DIMENSION(numVariaveis,numVariaveis) :: xx !Matriz sistema linearizado
		INTEGER :: i !Contador
!
		x = x0 !Atribui condicoes Iniciais
		xx = identidade !Condicao inicial do sistema linearizado
!
		espectroLocal = 0.d0 !Zera valor inicial a fim acumular a soma sem erros
!
		DO i = 1, pontos + transiente!Evolui sistema com descarte de transiente
			CALL SF_RK4(pLocal, x, xx) !Chama subrotina de integraçao
			CALL SF_OrtonormalizaGramSchmidt(norma, xx) !Chama ortonormalizacao de Gram-Schmidt
!			
			!Apos descarte do transiente, acumula a expansao dos vetores coluna de xx
			IF (i > transiente) WHERE (norma > 0.d0) espectroLocal = espectroLocal + DLOG(norma)
		END DO
!
		espectroLocal = espectroLocal/(passoRK4*DFLOAT(pontos)) !Divide a soma dos logs pelo tempo
	END SUBROUTINE SE_EvoluiSistema
!-Subrotinas de Ferramentas: SF_..--------------------------------------------------------------------------------------------------
	SUBROUTINE SF_RK4(pLocal, x, xx)!Integrador Runge-Kutta de quarta ordem-----------------------------------------------------------
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN), DIMENSION(numParametros) :: pLocal !Parametros
		DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numVariaveis) :: x !Estado do sistema
		DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numVariaveis,numVariaveis) :: xx !Matriz sistema linearizado
		DOUBLE PRECISION, DIMENSION(numVariaveis) :: k1, k2, k3, k4 !Coeficientes do metodo de integracao
		DOUBLE PRECISION, DIMENSION(numVariaveis,numVariaveis) :: kk1, kk2, kk3, kk4
!
		CALL SE_Sistema(pLocal, x, k1, xx, kk1) !Quatro chamadas conforme as etapas intermediarias do metodo RK4
		CALL SE_Sistema(pLocal, x + k1*passoRK4*0.5d0, k2, xx + kk1*passoRK4*0.5d0, kk2)
		CALL SE_Sistema(pLocal, x + k2*passoRK4*0.5d0, k3, xx + kk2*passoRK4*0.5d0, kk3)
		CALL SE_Sistema(pLocal, x + k3*passoRK4, k4, xx + kk3*passoRK4, kk4)
!
		x = x + passoRK4*(k1 + 2.d0*(k2 + k3) + k4)/6.d0
		xx = xx + passoRK4*(kk1 + 2.d0*(kk2 + kk3) + kk4)/6.d0
	END SUBROUTINE SF_RK4
	SUBROUTINE SF_OrtonormalizaGramSchmidt(norma, xx)!--------------------------------------------------------------------------------
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(OUT), DIMENSION(numVariaveis) :: norma
		DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numVariaveis,numVariaveis) :: xx
		INTEGER :: i, j
!		
		DO i = 1, numVariaveis !Loop nos vetores coluna da matriz xx
			DO j = 1, i - 1 !Subtrai projecões do vetor coluna i nos vetores coluna anteriores (j < i)
				xx(:,i) = xx(:,i) - xx(:,j)*SUM(xx(:,i)*xx(:,j)) 
			END DO
!			
			norma(i) = NORM2(xx(:,i)) !Determina norma euclidiana do vetor coluna
			xx(:,i) = xx(:,i)/norma(i) !Normaliza vetor
		END DO	
	END SUBROUTINE SF_OrtonormalizaGramSchmidt
!-Subrotinas Auxiliares: SA_..------------------------------------------------------------------------------------------------------
	SUBROUTINE SA_DiscretizaParametro(elementos, minimo, maximo, vetor)!--------------------------------------------------------------
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: elementos
		DOUBLE PRECISION, INTENT(IN) :: minimo, maximo
		DOUBLE PRECISION, INTENT(OUT), DIMENSION(elementos) :: vetor
		DOUBLE PRECISION :: intervalo
		INTEGER :: i
!
		intervalo = (maximo - minimo)/DFLOAT(elementos - 1)
		FORALL(i = 1:elementos) vetor(i) = minimo + intervalo*DFLOAT(i - 1) !Pode ficar mais eficiente: DO CONCURRENT
	END SUBROUTINE SA_DiscretizaParametro
!-Subrotinas de dados: SD_..--------------------------------------------------------------------------------------------------------	
	SUBROUTINE SD_SalvaEspectro(parametroA)!------------------------------------------------------------------------------------------
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN), DIMENSION(grade) :: parametroA
		INTEGER :: i, j
!
		OPEN (101, FILE = './'//TRIM(pastaSaida)//'/'//TRIM(arquivoSaida)//'.dat', STATUS = 'unknown')
!
		!Cabecalho
		WRITE(101, '(A14)', ADVANCE = 'no') '#parametro'
		DO i = 1, numVariaveis 
			WRITE(101, '(A12,I2.2)', ADVANCE = 'no') 'lambda_', i
		END DO
		WRITE(101, *)
!		 
		DO i = 1, grade
			WRITE(101, '(F14.6)', ADVANCE = 'no') parametroA(i)
!
			DO j = 1, numVariaveis
				WRITE(101, '(F14.6)', ADVANCE = 'no') espectroLyapunov(j,i)
			END DO
!
			WRITE(101, *)
		END DO
!			
		CLOSE(101)
	END SUBROUTINE SD_SalvaEspectro
END PROGRAM LyapunovEspectro_EixoParametro_2023MAR_V02
