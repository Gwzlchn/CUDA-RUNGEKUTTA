

program read
implicit none
DOUBLE PRECISION  :: mat(10000,12)
INTEGER:: i
open(12 , File = 'Initialization.dat')
do i = 1,10000
    !TODO_statement
    read(12,*) mat(i,:)
enddo


write(*,*) mat(1000,:)
write(*,*) mat(2000,:)


end program read