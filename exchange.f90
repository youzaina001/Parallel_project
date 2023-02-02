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

    if(ar == 0) then

      if (me /= (Np-1)) then

        call MPI_Send(x0(iEnd-2*OVERLAP*Nx+1:iEnd-2*OVERLAP*Nx+Nx),Nx,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,statinfo)

      end if

      if (me /= 0) then

        call MPI_Recv(Up,Nx,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,status,statinfo)

      end If

      if (me /= 0) then

        call MPI_Send(x0(iBeg+2*OVERLAP*Nx-Nx:iBeg+2*OVERLAP*Nx-1),Nx,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= (Np-1)) then
        
        call MPI_Recv(Down,Nx,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,status,statinfo)
      
      end if

    else

      if (me /= (Np-1)) then
        
        call MPI_Send(x0(iEnd-2*OVERLAP*Nx-Nx+1:iEnd-2*OVERLAP+2*Nx),3*Nx,MPI_DOUBLE,me+1,1,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= 0) then

        call MPI_Recv(Up,3*Nx,MPI_DOUBLE,me-1,1,MPI_COMM_WORLD,status,statinfo)
      
      end if

      if (me /= 0) then

        call MPI_Send(x0(iBeg+2*OVERLAP*Nx-2*Nx:iBeg+2*OVERLAP*Nx+Nx-1),3*Nx,MPI_DOUBLE,me-1,0,MPI_COMM_WORLD,statinfo)
      
      end if

      if (me /= (Np-1)) then
        
        call MPI_Recv(Down,3*Nx,MPI_DOUBLE,me+1,0,MPI_COMM_WORLD,status,statinfo)
      
      end if

    end if

  end subroutine communications

end module exchange
