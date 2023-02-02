program main

  use parameters
  use variables
  use exchange
  use functions
  use charge
  use matrix
  use solver
  use MPI

  implicit none

  call MPI_Init(statinfo)

  call MPI_Comm_size(MPI_COMM_WORLD,Np,statinfo)
  call MPI_Comm_rank(MPI_COMM_WORLD,Me,statinfo)

  call charge_proc(Ny,Np,me,iBeg,iEnd)
  iBeg = (iBeg - 1) * Nx + 1
  iEnd = iEnd * Nx
  print*, iBeg, iEnd, Me

  ! Allocating data
  call data_allocation(iBeg,iEnd,u0,u,b,Up,Down,Up_o,Down_o)

  ! Initialisation
  u0(iBeg:iEnd) = 0._PR

  ! First loop : tima loop
  do i = 1, nmax

    iteration = 0
    ERROR = 1._PR
    Up_o = 1._PR
    Down_o = 1._PR

    do while ((error .ge. eps_scwz) .and. (iteration .le. nmax_scwz))

      call communications(u0,iBeg,iEnd,me,Np,Up,Down)

      call RHS(u0,(i+1)*dt,dx,dy,b,iBeg,iEnd,me,Np,Up,Down)

      call BiCGSTAB(alpha,beta,gamma,b,u,iBeg,iEnd,dx,dy,me,Np)

      ERROR = schwarz_error(u,iBeg,iEnd,me,Np,Up_o,Down_o)

      u0(iBeg:iEnd) = u(iBeg:iEnd)

      Up_o = Up
      Down_o = Down

      iteration = iteration + 1

    end do

    u0(iBeg:iEnd) = u(iBeg:iEnd)
  
  end do

  ! Writing the local solution
  Call output(me,filename,iBeg,iEnd,dx,dy,u)

  deallocate(u0,u,b,Up,Up_o,Down,Down_o)

  call MPI_Finalize(statinfo)

end program main
