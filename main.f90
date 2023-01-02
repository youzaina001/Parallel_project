program main

  use parameters
  use exchange
  use functions
  use charge
  use matrix
  use solver
  use MPI

  implicit none

  real(PR) :: dx, dy, alpha, beta, gamma, ERROR
  real(PR), dimension(:), allocatable :: u, u0, b, x, Up, Up_o, Down, Down_o
  integer :: iBeg, iEnd, me, Np, statinfo, i, max, iteration
  character(13) :: filename

  call MPI_Init(statinfo)

  call MPI_Comm_size(MPI_COMM_WORLD,Np,statinfo)
  call MPI_Comm_rank(MPI_COMM_WORLD,Me,statinfo)

  ! Reading parameters
  call load_parameters()

  call charge_proc(Ny,Np,me,iBeg,iEnd)

  ! Computing dx and dy
  dx = Lx/(Nx+1)
  dy = Ly/(Ny+1)

  ! Loading the matrix's coefficients
  Call mat_coeff(dx,dy,alpha,beta,gamma)

  allocate(u0(iBeg:iEnd))
  allocate(u(iBeg:iEnd))
  allocate(b(iBeg:iEnd))

  if(alpha_robin == 0) then

    allocate(Up(1:Nx)); allocate(Up_o(1:Nx))
    allocate(Down(1:Nx)); allocate(Down_o(1:Nx))
  
  else

    allocate(Up(1:3*Nx)); allocate(Up_o(1:3*Nx))
    allocate(Down(1:3*Nx)); allocate(Down_o(1:3*Nx))
  
  end if

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
