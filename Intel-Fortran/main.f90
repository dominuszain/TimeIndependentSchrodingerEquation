module useful_subs

    implicit none
    
    contains
        subroutine linspace(startpoint, endpoint, num, array, step)
        
            implicit none
            
            real, intent(in) :: startpoint, endpoint
            integer, intent(in) :: num
            real, dimension(:), intent(out) :: array
            real, intent(out) :: step
            
            integer :: i
            
            step = (endpoint - startpoint) / (num - 1)
            
            do i = 1, num
            
                array(i) = startpoint + (i - 1) * step
                
            end do
            
        end subroutine
        
        subroutine v(x, y)

            implicit none
    
            real, intent(in) :: x
            real, intent(out) :: y
    
            y = -1 / x

        end subroutine
        
        subroutine matrix_generator(num, stepsize, m, h, potentiarray, arr)
        
            implicit none
            
            integer, intent(in) :: num
            real, intent(in) :: stepsize, m, h
            
            real, dimension(:), intent(in) :: potentiarray
            real, dimension(:, :), intent(out) :: arr
            
            integer :: i, j
            
            do i = 1, num - 2

                arr(i, i) = 2 + 2 * m * stepsize * stepsize * potentiarray(i) / (h * h)
                
                if (i < (num - 2)) then
                    
                    arr(i + 1, i) = -1
                    arr(i, i + 1) = -1
                
                end if

            end do
            
            !do i = 1, num - 2
            
            !    print *, (arr(i, j), j = 1, num - 2)
            
            !end do
        
        end subroutine
        
        subroutine eig(matrix, n, w)

            implicit none
            
            real, dimension(:, :), intent(inout) :: matrix
            character(len = 1) :: jobz = 'V'
            character(len = 1) :: uplo = 'L'
            integer, intent(in) :: n
            integer :: lda
            real, dimension(:), intent(out) :: w
            real, dimension(:), allocatable :: work
            integer :: lwork
            integer :: info
            
            integer :: i, j
            
            lda = n
            lwork = 10 * n
            
            allocate(work(lwork))       
            
            ! eigenvalues for symmetric matrices
            
            call dsyev(jobz, uplo, n, matrix, lda, w, work, lwork, info)
            
            deallocate(work)

        end subroutine

end module

program tise_on_steroids
    
    use useful_subs
    
    implicit none
    
    real :: startpoint, endpoint, step, m, h
    integer :: num, i, j
    real, dimension(:), allocatable :: array, potentiarray, eigenvalues
    real, dimension(:, :), allocatable :: matrix
    
    m = 1
    h = 1
    
    startpoint = 1
    endpoint = 100
    num = 1000
    
    allocate(array(num), potentiarray(num), matrix(num - 2, num - 2), eigenvalues(num - 2))
    
    call linspace(startpoint, endpoint, num, array, step)
    
    do i = 1, num
    
        call v(array(i), potentiarray(i))
        !print *, array(i), potentiarray(i)
    
    end do
    
    call matrix_generator(num, step, m, h, potentiarray, matrix)
    
    !print *
    !do i = 1, num - 2
    
    !    print *, (matrix(i, j), j = 1, num - 2)
    
    !end do
    
    call eig(matrix, num - 2, eigenvalues)
    
    open(1, file = "eigenout.txt")
    !print *
    !print *, eigenvalues
    write (1, *) eigenvalues
    write (1, *)
    
    !print *
    do i = 1, num - 2
        
    !    print *, (matrix(i, j), j = 1, num - 2)
        write (1, *) (matrix(i, j), j = 1, num - 2)
    
    end do
    close(1)
    
    deallocate(array, potentiarray, matrix)
    
end program tise_on_steroids
