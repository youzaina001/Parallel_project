module matrix

  use MPI
  use parameters
  use functions
  use charge

  implicit none

contains

  Function matrix_vector_product(alpha,beta,gamma,x,iBeg,iEnd,dx,dy,me,Np) Result(res)

    integer, intent(in) :: iBeg, iEnd, me, Np
    real(PR), intent(in) :: alpha, beta, gamma, dx, dy
    real(PR), dimension(iBeg:iEnd), intent(in) :: x
    real(PR), dimension(iBeg:iEnd) :: res
    integer :: i

    ! Initialisation du vecteur produit
    res = 0._PR

    do i = iBeg + 1, iEnd - 1

      res(i) = res(i) + alpha*x(i)

      if (modulo(i,Nx) == 1) then

        res(i) = res(i) + beta*x(i+1)
          
      else if (modulo(i,Nx) == 0) then
             
        res(i) = res(i) + beta*x(i-1)
          
      else

        res(i) = res(i) + beta*(x(i+1) + x(i-1))
          
      end if
    
    end do

    res(iBeg) = res(iBeg) + alpha*x(iBeg)
    if (modulo(iBeg,Nx) == 1) then

      res(iBeg) = res(iBeg) + beta*x(iBeg+1)
    
    end if

    res(iEnd) = res(iEnd) + alpha*x(iEnd)
    if (modulo(iEnd,Nx) == 0) then

      res(iEnd) = res(iEnd) + beta*x(iEnd-1)
    
    end if

    do i = 1, iEnd-iBeg-Nx+1

      res(i+iBeg-1) = res(i+iBeg-1) + gamma*x(i+iBeg+Nx-1)
    
    end do

    do i = iBeg-iEnd+2*Nx, Nx

      res(iEnd-Nx+i) = res(iEnd-Nx+i) + gamma*x(iEnd-2*Nx+i)
    
    end do

    if (ar == 0) then

      if (me /= 0) then

        do i = iBeg, iBeg+Nx-1
          
          res(i) = x(i)/(dx**2+dy**2)
        
        end do

      end if

      if (me /= (Np-1)) then

        do i = iEnd-Nx+1, iEnd

          res(i) = x(i)/(dx**2+dy**2)
        
        end do
      
      end if
    
    else

      if(me /= 0) then

        do i = iBeg, iBeg+Nx-1

          res(i) = (1._PR + 2._PR*D*dt*(1._PR/(dx**2)+1._PR/(dy**2)+br/(ar*dy)))*x(i)
        
        end do

        do i = iBeg+1, iBeg+Nx-2

          res(i) = res(i) - D*dt*(x(i-1)+x(i+1))/(dx**2)
        
        end do

        res(iBeg) = res(iBeg) - D*dt*x(iBeg+1)/(dx**2)
        res(iBeg+Nx-1) = res(iBeg+Nx-1) - D*dt*x(iBeg+Nx-2)/(dx**2)

        do i = iBeg+Nx, iBeg+2*Nx-1

          res(i-Nx) = res(i-Nx) - 2._PR*D*dt/(dy**2)*x(i)
        
        end do

      end if

      if(me /= (Np-1)) then

        do i = iEnd-Nx+1, iEnd

          res(i) = (1._PR+2._PR*D*dt*(1._PR/(dx**2)+1._PR/(dy**2)+br/(ar*dy)))*x(i)
        
        end do

        do i = iEnd-Nx+2, iEnd-1
          
          res(i) = res(i) - D*dt*(x(i-1)+x(i+1))/(dx**2)
        
        end do

        res(iEnd-Nx+1) = res(iEnd-Nx+1)-D*dt*x(iEnd-Nx+2)/(dx**2)
        res(iEnd) = res(iEnd)-D*dt*x(iEnd-1)/(dx**2)

        do i = iEnd-2*Nx+1, iEnd-Nx
          
          res(i+Nx) = res(i+Nx)-2._PR*D*dt/(dy**2)*x(i)
        
        end do

      end if

    end if

  end function matrix_vector_product

  subroutine RHS(x0,t,dx,dy,b,iBeg,iEnd,me,Np,Up,Down)

    integer, intent(in) :: iBeg, iEnd, me, Np
    real(PR), dimension(iBeg:iEnd), intent(in) :: x0
    real(PR), intent(in) :: t, dx, dy
    real(PR), dimension(iBeg:iEnd), intent(inout) :: b
    real(PR), dimension(:), intent(in) :: Up, Down
    integer :: i, j, i1, j1

    do i = iBeg, iEnd

      i1 = modulo(i-1,Nx) + 1
      j1 = ((i-1)/(Nx))+1

      b(i) = x0(i) + dt*f(i1*dx,j1*dy,t,BC_choice)

      if ((i1 == 1) .or. (i1 == Nx)) then

        b(i) = b(i) + D*dt*h(i1*dx,j1*dy,t,BC_choice)/(dx**2)
      
      end if

      if ((j1 == 1) .or. (j1 == Ny)) then

        b(i) = b(i) + D*dt*g(i1*dx,j1*dy,t,BC_choice)/(dy**2)
      
      end if

    end do

    if(ar == 0) then

      if (me /= 0) then

        do i = iBeg, iBeg+Nx-1

          b(i) = Up(i-iBeg+1)/(dx**2+dy**2)
        
        end do
      
      end if

      if (me /= (Np-1)) then

        do i = iEnd-Nx+1, iEnd

          b(i) = Down(i-iEnd+Nx)/(dx**2+dy**2)
        
        end do

      end if

    else

      if (Me /= 0) then

        do i = iBeg, iBeg+Nx-1

          b(i) = b(i) + D*dt*(2._PR*br*Up(i-iBeg+1+Nx)/(ar*dy)+(Up(i-iBeg+1)-Up(i-iBeg+1+2*Nx))/(dy)**2)
        
        end do

      end if

      if (me /= (Np-1)) then

        do i = iEnd-Nx+1, iEnd
          
          b(i) = b(i) + D*dt*(2._PR*br*Down(i-iEnd+2*Nx)/(ar*dy)+(Down(i-iEnd+3*Nx)-Down(i-iEnd+Nx))/(dy)**2)
        
        end do
      
      end if
    
    end if

  end subroutine RHS

end module matrix
