module zscore

  use types, only: dp
  use permutations, only: permutation, permuteColumns
  implicit none

contains

  !function computeZ(chrom1, exp1, chrom2, exp2) result(z)
  !  integer, dimension(:, :), intent(in) :: chrom1
  !  real(dp), dimension(size(chrom1, 1), size(chrom1, 2)), intent(in) :: exp1
  !  integer, dimension(:, :), intent(in) :: chrom2
  !  real(dp), dimension(size(chrom2, 1), size(chrom2, 2)), intent(in) :: exp2
  !
  !  real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: z
  !
  !  !integer, dimension(size(z, 1), size(z, 2)) :: observed
  !  !real(dp), dimension(size(z, 1), size(z, 2)) :: expected
  !  integer, dimension(size(chrom1, 1), size(chrom2, 1)) :: observed
  !  real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: expected
  !
  !  observed = matmul(chrom1, transpose(chrom2))
  !  expected = matmul(exp1, transpose(exp2))
  !
  !  z = (observed - expected) / sqrt(expected * (1.0d0 - expected / dble(size(chrom1, 2))))
  !end function computeZ

  function computeZ(chrom1, p1, chrom2, p2) result(z)
    integer, dimension(:, :), intent(in) :: chrom1
    real(dp), dimension(size(chrom1, 1), size(chrom1, 2)), intent(in) :: p1
    integer, dimension(:, :), intent(in) :: chrom2
    real(dp), dimension(size(chrom2, 1), size(chrom2, 2)), intent(in) :: p2

    real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: z

    integer, dimension(size(chrom1, 1), size(chrom2, 1)) :: observed
    real(dp), dimension(size(chrom1, 1), size(chrom2, 1)) :: expected, stdev

    expected = matmul(p1, transpose(p2))
    stdev = sqrt(expected - matmul(p1**2, transpose(p2)**2))
    observed = matmul(chrom1, transpose(chrom2))

    z = (observed - expected) / stdev
  end function computeZ

  function sampleFromNull(chrom1, exp1, chrom2, exp2, numSamples) result (nullSamples)
    integer, dimension(:, :), intent(in) :: chrom1
    real(dp), dimension(size(chrom1, 1), size(chrom1, 2)), intent(in) :: exp1
    integer, dimension(:, :), intent(in) :: chrom2
    real(dp), dimension(size(chrom2, 1), size(chrom2, 2)), intent(in) :: exp2
    integer, intent(in) :: numSamples

    real(dp), dimension(numSamples) :: nullSamples

    integer :: i
    integer, dimension(size(chrom1, 2)) :: perm
    integer, dimension(size(chrom2, 1), size(chrom2, 2)) :: permutedChrom2
    real(dp), dimension(size(chrom2, 1), size(chrom2, 2)) :: permutedExp2

    call setSeed()

    do i = 1, numSamples
       perm = permutation(size(perm))
       permutedChrom2 = permuteColumns(chrom2, perm)
       permutedExp2 = permuteColumns(exp2, perm)
       nullSamples(i) = maxval(computeZ(chrom1, exp1, permutedChrom2, permutedExp2))
    end do
  end function sampleFromNull

  subroutine setSeed()
    integer :: size, i
    integer, allocatable, dimension(:) :: seed

    call random_seed(size=size)
    allocate(seed(size))

    do i = 1, size
       seed(i) = 1234 * i
    end do

    call random_seed(put=seed)
    deallocate(seed)
  end subroutine setSeed

end module zscore
