! code of Zain Ul Abideen
! Solution to Time Independent Schrodinger Equation in Fortran95
! UPGRADED VERSION.

program tise_on_steroids

    implicit none
    
    real :: startpoint, endpoint, iter, stepsize, m, h
    character(len = 11) :: path = "matrix.fout"
    integer :: num, i, o, k
    
    real, dimension(:), allocatable :: array, potentiarray
    
    real, dimension(:, :), allocatable :: arr
    
    startpoint = 1
    endpoint = 10
    
    ! 91 for a neat linear space of values.
    num = 91
    
    m = 1
    h = 1
    
    allocate(array(num))
    allocate(potentiarray(num))
    allocate(arr(num - 2, num - 2))
    
    call linspace(startpoint, endpoint, num, array, stepsize)
    
    do i = 1, num
        
        call v(array(i), iter)
        potentiarray(i) = iter
    
    end do
    
    !print *, array
    !print *
    !print *, potentiarray
    !print *
    
    do i = 1, num - 2
    
        arr(i, i) = 2 + 2 * m * stepsize * stepsize * potentiarray(i) / (h * h) 
        
        if (i < (num - 2)) then
            
            arr(i + 1, i) = -1
            arr(i, i + 1) = -1
            
        end if
    
    end do
    
    !print *, arr
    
    open(1, file = "xvals.fout")
    
        do i = 1, num
        
            write (1,*) array(i)
        
        end do
    
    close(1)
    
        
    open(1, file = "potvals.fout")
    
        do i = 1, num
        
            write (1,*) potentiarray(i)
        
        end do
    
    close(1)
    
    open(1, file = "matrix.vout")
    open(2, file = "matrix.fout")
    
    do k = 1, num - 2
        
        write (1,*) (arr(o, k), o = 1, num - 2)
        
        do o = 1, num - 2
        
            write (2,*) arr(o, k)
        
        end do
        
    end do
    
    close(1)
    close(2)
    
    call eig(num - 2, path)
    
    deallocate(array)
    deallocate(potentiarray)
    deallocate(arr)
    
end program tise_on_steroids

subroutine v(x, y)

    implicit none
    
    real, intent(in) :: x
    real, intent(out) :: y
    
    y = -1 / x

end subroutine

subroutine linspace(startpoint, endpoint, num, answer, step)
    implicit none
    
    real, intent(in) :: startpoint, endpoint
    integer, intent(in) :: num
    real, intent(out) :: answer(num)
    real, intent(out) :: step
    integer :: i
  
    step = (endpoint - startpoint) / (num - 1)
  
    do i = 1, num
        answer(i) = startpoint + (i - 1) * step
    end do
    
end subroutine

subroutine eig(k, passedpath)

    implicit none
    
    ! declaring the variables for eigenvalues and vectors
    real, dimension(:, :), allocatable :: A
    integer, intent(in) :: k
    integer :: m, o
    character(len = 1) :: jobz = 'V'
    character(len = 1) :: uplo = 'L'
    character(len = 11), intent(in) :: passedpath
    character(len = 11) :: path
    integer :: n
    integer :: lda
    real, dimension(:), allocatable :: w, r
    real, dimension(:), allocatable :: work
    integer :: lwork
    integer :: info
    
    n = k
    lda = k
    lwork = 10 * n
    
    path = passedpath
    
    allocate(A(k, k))
    allocate(w(k))
    allocate(r(k * k))
    allocate(work(lwork))
    
    ! reading the array to calculate the eigen stuff
    open(1, file = path)
    
    do m = 1, k * k
    
        read (1,*) r(m)
        
    end do
    
    ! converting the linear array to square matrix
    A = reshape(r, (/k, k/))
    
    close(1)
    
    ! eigenvalues and vectors lapack subroutine for symmetric matrices
    call dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
    
    ! transposition for output vectors as columns
    A = transpose(A)
    
    ! writing the eigenvalues and vectors to the final output
    open(1, file = "result.fout")
    
        write (1, *) w
        write (1, *)
        
        do m = 1, k
        
            write (1, *) (A(o, m), o = 1, k)
            
        end do
    
    close(1)
    
    deallocate(A)
    deallocate(w)
    deallocate(r)
    deallocate(work)

end subroutine
