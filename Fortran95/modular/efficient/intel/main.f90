module eigen

    implicit none

    contains

        subroutine eig(inpath, outpath, n)

            implicit none
            
            real, dimension(:, :), allocatable :: matrix
            character(len = 20), intent(in) :: inpath, outpath
            integer :: i, j
            character(len = 1) :: jobz = 'V'
            character(len = 1) :: uplo = 'L'
            integer, intent(in) :: n
            integer :: lda
            real, dimension(:), allocatable :: w
            real, dimension(:), allocatable :: work
            integer :: lwork
            integer :: info
            
            lda = n
            lwork = 10 * n
            
            allocate(matrix(n, n), w(n), work(lwork))
            
            open(1, file = inpath)

            do i = 1, n
            
                read (1,*) matrix(n, :)
                
            end do

            close(1)           
            
            ! eigenvalues for symmetric matrices
            
            call dsyev(jobz, uplo, n, matrix, lda, w, work, lwork, info)
            
            open(1, file = outpath)
                
                write (1, *) "Eigenvalues: "
                write (1, *)
                
                do i = 1, n
                
                    write (1, *) i, w(i)
                    
                end do

                write (1, *)
                write (1, *) "Eigenvectors: " 
                write (1, *)               

                do i = 1, n
                    do j = 1, n
                    
                        write (1, *) i, matrix(j, i)
                        
                    end do
                    write (1, *)
                    
                end do
            
            close(1)
            
            deallocate(matrix, w, work)

        end subroutine
        
        subroutine v(x, y)

            implicit none
    
            real, intent(in) :: x
            real, intent(out) :: y
    
            y = -1 / x

        end subroutine
        
        subroutine linspace(startpoint, endpoint, num, outpath, step)
        
            implicit none
            
            real, intent(in) :: startpoint, endpoint
            integer, intent(in) :: num
            character(len = 20), intent(in) :: outpath
            real, intent(out) :: step
            
            real, dimension(:), allocatable :: answer
            integer :: i
            
            allocate(answer(num))
          
            step = (endpoint - startpoint) / (num - 1)
            
            open(1, file = outpath)
            
            do i = 1, num
            
                answer(i) = startpoint + (i - 1) * step
                write (1, *) answer(i)
                
            end do
            
            close(1)
            
            deallocate(answer)
            
        end subroutine
        
        subroutine matrix_generator(num, stepsize, m, h, inpath, outpath)
        
            implicit none
            
            integer, intent(in) :: num
            real, intent(in) :: stepsize, m, h
            character(len = 20), intent(in) :: inpath, outpath
            
            real, dimension(:), allocatable :: potentiarray
            real, dimension(:, :), allocatable :: arr
            
            integer :: i, j
            
            allocate(potentiarray(num), arr(num-2, num-2))
            
                open(1, file = inpath)
    
                    do i = 1, num
                    
                        read (1, *) potentiarray(i)
                    
                    end do
                
                close(1)
            
                do i = 1, num - 2
    
                    arr(i, i) = 2 + 2 * m * stepsize * stepsize * potentiarray(i) / (h * h) 
                    
                    if (i < (num - 2)) then
                        
                        arr(i + 1, i) = -1
                        arr(i, i + 1) = -1
                    
                    end if
    
                end do
                
                open(1, file = outpath)
                
                    do i = 1, num - 2
                
                        write (1, *) (arr(i, j), j = 1, num - 2)
                
                    end do                
                
                close(1)
                
            deallocate(potentiarray, arr)
        
        
        end subroutine

end module eigen

program eigenvalue

    use eigen    

    implicit none
    
    real :: startpoint, endpoint, stepsize, m, h
    real, dimension(:), allocatable :: array, potentiarray
    character(len = 20) :: arraypath, potentiarraypath
    character(len = 20) :: matinpath, matoutpath
    character(len = 20) :: eigeninpath, eigenoutpath
    integer, parameter :: num = 52
    integer :: i, j
    
    arraypath = "array.out"
    potentiarraypath = "potentiarray.out"
    
    matinpath = "potentiarray.out"
    matoutpath = "matrix-gen.out"
    
    eigeninpath = "matrix-gen.out"
    eigenoutpath = "eigenstuff.out"
    
    m = 1
    h = 1
    
    startpoint = 1
    endpoint = 10
    
    allocate(array(num), potentiarray(num))
    
    call linspace(startpoint, endpoint, num, arraypath, stepsize)
    
    open(1, file = arraypath)
    
        do i = 1, num
        
            read (1, *) array(i)
        
        end do
    
    close(1)
    
    do i = 1, num
    
        call v(array(i), potentiarray(i))
    
    end do
    
    open(1, file = potentiarraypath)
    
        do i = 1, num
        
            write (1, *) potentiarray(i)
        
        end do
    
    close(1)
    
    call matrix_generator(num, stepsize, m, h, matinpath, matoutpath)
    
    call eig(eigeninpath, eigenoutpath, num-2)
    
    deallocate(array, potentiarray)

end program eigenvalue
