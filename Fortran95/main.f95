! code of Zain Ul Abideen
! Solution to Time Independent Schrodinger Equation in Fortran95

program tise

    implicit none
    
    ! declaring variables and arrays
    real :: m, h, d
    integer :: i, a, n, k, o
    real, dimension(:), allocatable :: x, vi
    real, dimension(:,:), allocatable :: arr
    
    real :: iter
    
    m = 1
    h = 1
    
    ! a is the upper limit, change to whatever
    a = 10
    
    ! d is the step size, dont mess with it
    d = 0.1
    
    ! n is the total number of points, f(d)
    n = a * 10 - 9
    
    ! allocating the arrays and matrices
    allocate(x(n))
    allocate(vi(n))
    allocate(arr(n - 2, n - 2))
    
    ! populating the position and potential arrays
    do i = 10, a * 10
        
        x(i - 9) = real(i) * d
        call v(x(i - 9), iter)
        vi(i - 9) = iter
    
    end do
    
    ! generating the matrix
    do i = 1, n - 2
    
        arr(i, i) = 2 + 2 * m * d * d * vi(i) / (h * h) 
        
        if (i < (n - 2)) then
            
            arr(i + 1, i) = -1
            arr(i, i + 1) = -1
            
        end if
    
    end do
    
    ! printing the position values for later plots    
    open(1, file = "xvals.fout")
    
        do i = 1, n
        
            write (1,*) x(i)
        
        end do
    
    close(1)
    
    open(1, file = "matrix.vout")
    open(2, file = "matrix.fout")
    
    ! printing the generated matrix in two formats
    ! .vout for visualization and .fout for later input
    do k = 1, n - 2
        
        ! first rows will be filled, then columns
        write (1,*) (arr(o, k), o = 1, n - 2)
        
        do o = 1, n - 2
        
            write (2,*) arr(o, k)
        
        end do
        
    end do
    
    close(1)
    close(2)
    
    ! calling the symmetric eigenvalues and vectors subroutine
    call eig(n - 2)
    
    deallocate(x)
    deallocate(vi)
    deallocate(arr)

end program tise

subroutine v(x, y)

    implicit none
    
    real, intent(in) :: x
    real, intent(out) :: y
    
    ! change the potential here to whatever
    y = x * x

end subroutine

subroutine eig(k)

    implicit none
    
    ! declaring the variables for eigenvalues and vectors
    real, dimension(:, :), allocatable :: A
    integer :: m, o, k
    character(len = 1) :: jobz = 'V'
    character(len = 1) :: uplo = 'L'
    integer :: n
    integer :: lda
    real, dimension(:), allocatable :: w, r
    real, dimension(:), allocatable :: work
    integer :: lwork
    integer :: info
    
    n = k
    lda = k
    lwork = 10 * n
    
    allocate(A(k, k))
    allocate(w(k))
    allocate(r(k * k))
    allocate(work(lwork))
    
    ! reading the array to calculate the eigen stuff
    open(1, file = "matrix.fout")
    
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
