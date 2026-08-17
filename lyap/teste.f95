program teste
    implicit none 
    real(8) :: W(2,2)

    W(1,1) = 1; W(1,2) = 2
    W(2,1) = 3; W(2,2) = 4

    print*, W(:,1)

    print*, W(1, :)



end program teste