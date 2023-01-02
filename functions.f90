module functions

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

end module functions
