module functions

  use MPI
  use parameters

  implicit none

contains

  function f(x,y,t,cas)

    integer, intent(in) :: cas
    real(PR), intent(in) :: x, y, t
    real(PR) :: f

    select case(cas)

    case(1)

      f = 2._PR * (y - y**2 + x - x**2)

    case(2)

      f = sin(x) + cos(y)

    case(3)

      f = exp(-(x-(Lx/(2._PR)))**2)*exp(-(y-(Ly/(2._PR)))**2)*cos(pi*t/(2._PR))
    
    case default

      print*, "Error : this case is not implemented"
    
    end select

  end function f

  function g(x,y,t,cas)

    integer, intent(in) :: cas
    real(PR), intent(in) :: x, y, t
    real(PR) :: g

    select case(cas)

    case(1)

      g = 0._PR
    
    case(2)
      
      g = sin(x) + cos(y)
    
    case(3)
      
      g = 0._PR
    
    case default

      print*, "Error : this case is not implemented"
    
    end select

  end function g

  function h(x,y,t,cas)

    integer, intent(in) :: cas
    real(PR), intent(in) :: x, y, t
    real(PR) :: h

    select case(cas)

    case(1)

      h = 0._PR
    
    case(2)
      
      h = sin(x) + cos(y)
    
    case(3)
      
      h = 1._PR
    
    case default

      print*, "Error : this case is not implemented"
    
    end select

  end function h

  function loc_dot_product(x,y,iBeg,iEnd) result(z)

    integer, intent(in) :: iBeg, iEnd
    real(PR), dimension(iBeg:iEnd), intent(in) :: x, y
    real(PR) :: z
    integer :: i

    ! Initialisation
    z = 0

    do i = iBeg, iEnd
      
      z = z + x(i)*y(i)
    
    end do

  end function loc_dot_product

  subroutine output(me,filename,iBeg,iEnd,dx,dy,u)

    real(PR), dimension(iBeg:iEnd), intent(in) :: u
    real(PR), intent(in) :: dx, dy
    integer, intent(in) :: me, iBeg, iEnd
    character(13), intent(out) :: filename
    character(3) :: num
    integer :: i0, i2, i3, i

    i0 = me/100
    i2 =(me-100*i0)/10
    i3 = me-100*i0-10*i2
    num = char(i0+48)//char(i2+48)//char(i3+48)
    filename = "sol"//num//".txt"

    open(2,file = filename)

    do i = iBeg, iEnd

      write(2,*) (modulo(i-1,Nx)+1)*dx, (((i-1)/(Nx))+1)*dy, u(i)
    
    end do

    close(2)

  end subroutine output

  function schwarz_error(x,iBeg,iEnd,me,Np,Up,Down) result(ERROR_glob)

    integer, intent(in) :: iBeg, iEnd, me, Np
    real(PR), dimension(iBeg:iEnd), intent(in) :: x
    real(PR), dimension(:), intent(in) :: Up, Down
    integer :: statinfo
    real(PR) :: ERROR_l, ERROR_r, ERROR_loc, ERROR_glob
    real(PR), dimension(:), allocatable    :: Up_t, Down_t

    allocate(Up_t(1:Nx),Down_t(1:Nx))

    if(ar .eq. 0) then

      Up_t = Up
      Down_t = Down

    else

      Up_t = Up(Nx+1:2*Nx)
      Down_t = Down(Nx+1:2*Nx)
    
    end if

    if (me == 0) then

      ERROR_l = 0._PR
      ERROR_r = maxval(abs(x(iEnd-Nx+1:iEnd)-Down_t))
    
    else if(me == (Np-1)) then

      ERROR_l = maxval(abs(x(iBeg:iBeg+Nx-1)-Up_t))
      ERROR_r = 0._PR
    
    else

      ERROR_l = maxval(abs(x(iBeg:iBeg+Nx-1)-Up_t))
      ERROR_r = maxval(abs(x(iEnd-Nx+1:iEnd)-Down_t))
    
    end if

    ERROR_loc = max(ERROR_l,ERROR_r)

    Call MPI_Allreduce(ERROR_loc,ERROR_glob,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,statinfo)

    deallocate(Up_t,Down_t)

  end function schwarz_error

  subroutine data_allocation(iBeg,iEnd,u0,u,b,Up,Down,Up_o,Down_o)

    integer :: iBeg, iEnd
    real(PR), dimension(:), allocatable :: u, u0, b, x, Up, Up_o, Down, Down_o

    allocate(u0(iBeg:iEnd))
    allocate(u(iBeg:iEnd))
    allocate(b(iBeg:iEnd))

    if(ar == 0) then

      allocate(Up(1:Nx)); allocate(Up_o(1:Nx))
      allocate(Down(1:Nx)); allocate(Down_o(1:Nx))
      
    else
    
      allocate(Up(1:3*Nx)); allocate(Up_o(1:3*Nx))
      allocate(Down(1:3*Nx)); allocate(Down_o(1:3*Nx))
      
    end if
    
  end subroutine data_allocation

end module functions
