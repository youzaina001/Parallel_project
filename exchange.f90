module exchange

  use MPI
  use parameters

  implicit none

contains

  subroutine communications(x0,iBeg,iEnd,me,Np,Up,Down)

    integer, intent(in) :: iBeg, iEnd, me, Np
    real(PR), dimension(iBeg:iEnd), intent(in) :: x0
    real(PR), dimension(:), intent(inout) :: Up, Down
    integer :: statinfo
    integer, dimension(MPI_STATUS_size) :: status

    if(alpha_robin == 0) then

      if (me /= (Np-1)) then

        call MPI_Send(x0(iEnd-2*rcv*Nx+1:iEnd-2*rcv*Nx+Nx),Nx,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,statinfo)

      end if

      if (me /= 0) then

        call MPI_Recv(Up,Nx,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,status,statinfo)

      end If

      if (me /= 0) then

        call MPI_Send(x0(iBeg+2*rcv*Nx-Nx:iBeg+2*rcv*Nx-1),Nx,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= (Np-1)) then
        
        call MPI_Recv(Down,Nx,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,status,statinfo)
      
      end if

    else

      if (me /= (Np-1)) then
        
        call MPI_Send(x0(iEnd-2*rcv*Nx-Nx+1:iEnd-2*rcv+2*Nx),3*Nx,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= 0) then

        call MPI_Recv(Up,3*Nx,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,status,statinfo)
      
      end if

      if (me /= 0) then

        call MPI_Send(x0(iBeg+2*rcv*Nx-2*Nx:iBeg+2*rcv*Nx+Nx-1),3*Nx,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= (Np-1)) then
        
        call MPI_Recv(Down,3*Nx,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,status,statinfo)
      
      end if

    end if

  end subroutine communications

  function schwarz_error(x,iBeg,iEnd,me,Np,Up,Down) result(ERROR_glob)

    integer, intent(in) :: iBeg, iEnd, me, Np
    real(PR), dimension(iBeg:iEnd), intent(in) :: x
    real(PR), dimension(:), intent(in) :: Up, Down
    integer :: statinfo
    real(PR) :: ERROR_l, ERROR_r, ERROR_loc, ERROR_glob
    real(PR), dimension(:), allocatable    :: Up_t, Down_t

    allocate(Up_t(1:Nx),Down_t(1:Nx))

    if(alpha_robin .eq. 0) then

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

end module exchange
