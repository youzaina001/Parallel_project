module variables

    use parameters

    implicit none
    
    real(PR) :: ERROR
    real(PR), dimension(:), allocatable :: u, u0, b, x, Up, Up_o, Down, Down_o
    integer :: iBeg, iEnd, me, Np, statinfo, i, max, iteration
    character(13) :: filename
    
end module variables
