module pscore

  use types, only: dp

  implicit none

contains

  function computeP(chrom1, p1, chrom2, p2) result(p)
    use poisbinom

    integer, dimension(:, :), intent(in) :: chrom1
    real(dp), dimension(size(chrom1, 1), size(chrom1, 2)), intent(in) :: p1
    integer, dimension(:, :), intent(in) :: chrom2
    real(dp), dimension(size(chrom2, 1), size(chrom2, 2)), intent(in) :: p2

    real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: p

    integer :: i, j

    !$omp parallel do collapse(2)
    do i = 1, size(chrom1, 1)
       do j = 1, size(chrom2, 1)
          p(i, j) = cdf(p1(i, :) * p2(j, :), dot_product(chrom1(i, :), chrom2(j, :)))
       end do
    end do
    !$omp end parallel do
  end function computeP

end module pscore
