module permutations

  use types, only: dp

  implicit none

  interface permuteColumns
     module procedure permuteColumns_i, permuteColumns_dp
  end interface permuteColumns

contains

  function permutation(n) result (p)
    !!! Copied and adapted from http://ftp.aset.psu.edu/pub/ger/fortran/hdk/byterand.f90
    !===Generate a random permutation, P, of the first N integers.
    !   (Shuffling, equivalent to sampling WITHOUT REPLACEMENT).
    !   Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
    integer n,p(n), k,j,i,ipj,itemp,m
    real u(100)
    do i=1,n
      p(i)=i
    end do
    !---Generate up to 100 U(0,1) numbers at a time.
    do i=1,n,100
      m=min(n-i+1,100)
      call random_number(u)
      do j=1,m
        ipj=i+j-1
        k=int(u(j)*(n-ipj+1))+ipj
        itemp=p(ipj)
        p(ipj)=p(k)
        p(k)=itemp
      end do
    end do
  end function permutation

  function permuteColumns_dp(matrix, perm) result(permuted)
    real(dp), dimension(:, :), intent(in) :: matrix
    integer, dimension(size(matrix, 2)), intent(in) :: perm
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) :: permuted

    integer :: i

    !do concurrent (i = 1:size(perm))
    do i = 1, size(perm)
       permuted(:, i) = matrix(:, perm(i))
    end do
  end function permuteColumns_dp

  function permuteColumns_i(matrix, perm) result(permuted)
    integer, dimension(:, :), intent(in) :: matrix
    integer, dimension(size(matrix, 2)), intent(in) :: perm
    integer, dimension(size(matrix, 1), size(matrix, 2)) :: permuted

    integer :: i

    !do concurrent (i = 1:size(perm))
    do i = 1, size(perm)
       permuted(:, i) = matrix(:, perm(i))
    end do
  end function permuteColumns_i

end module permutations
