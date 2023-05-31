module print
    implicit none
  
  contains


subroutine print_2d(arr, mes)

    implicit none
    integer, intent(in) :: arr(:, :)
    character(len=*) :: mes
    integer :: i, j

    ! Array is 2D, print as a matrix
    write(6,*) shape(arr)
    write(6,*) mes

    do i = 1, size(arr, 1)
       do j = 1, size(arr, 2)
          write(6, ' (F8.6, 3x)', advance='no')  arr(i,j)
       end do
       write(6,*)
    end do

 end subroutine print_2d

end module print