subroutine calculate_pseudo_dist(x, dist_matrix, N, T)
  implicit none
  real(8), dimension(N, T), intent(in) :: x   ! Input matrix (N x T)
  real(8), dimension(N, N), intent(out) :: dist_matrix  ! Output distance matrix (N x N)
  integer, intent(in) :: N  ! Number of data points (rows)
  integer, intent(in) :: T  ! Number of features (columns)
  
  integer :: i, j, k, m
  real(8) :: sum, max_dist, mean
  real(8), dimension(T) :: x_i, x_j, x_k
  integer :: i_block, j_block, i_block_end, j_block_end
  integer :: block_size
  
  ! Set block size, you can experiment with this value to find the best one for your system
  block_size = 64  ! Change as necessary
  
  ! Initialize distance matrix to 0
  dist_matrix = 0.0
  
  ! Calculate pseudo distance matrix using OpenMP for parallelism
  !$OMP PARALLEL DO PRIVATE(i, j, k, m, sum, max_dist, mean, x_i, x_j, x_k, i_block, j_block) SHARED(x, dist_matrix, N, T) 
  do i = 1, N, block_size
    do j = 1, N, block_size
      ! Block loop for i, j
      do i_block = i, min(i+block_size-1, N)
        do j_block = j, min(j+block_size-1, N)
          if (i_block /= j_block) then
            max_dist = -1.0E10  ! Initialize a large negative value for max dist

            ! Load data for x(i) and x(j) into temporary arrays
            x_i = x(i_block, :)
            x_j = x(j_block, :)

            ! Exclude i_block and j_block from k selection and compute the maximum distance
            do k = 1, N
              if (k /= i_block .and. k /= j_block) then
                x_k = x(k, :)

                sum = 0.0
                ! Compute the dot product of (x(i, m) - x(j, m)) with x(k, m)
                do m = 1, T
                  sum = sum + abs((x_i(m) - x_j(m)) * x_k(m))
                end do

                mean = sum / T
                if (mean > max_dist) then
                  max_dist = mean
                end if
              end if
            end do
            dist_matrix(i_block, j_block) = max_dist
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine calculate_pseudo_dist
