Module parameters

  implicit none

  integer, parameter :: PR = 8 ! Double précision
  real(PR), parameter :: pi = 4._PR*atan(1._PR)
  integer :: Nx, Ny, BC_choice, kmax, nmax, rcv, nmax_scwz
  real(PR) :: Lx, Ly, D, dt, eps, eps_scwz, alpha_robin, beta_robin

contains

  subroutine load_parameters()

    open(1, file = 'data.txt')

    read(1,*) Nx
    read(1,*) Ny
    read(1,*) Lx
    read(1,*) Ly
    read(1,*) D
    read(1,*) dt
    read(1,*) BC_choice
    read(1,*) kmax
    read(1,*) eps
    read(1,*) nmax
    read(1,*) rcv
    read(1,*) eps_scwz
    read(1,*) nmax_scwz
    read(1,*) alpha_robin
    read(1,*) beta_robin

    close(1)

  end Subroutine load_parameters

end module parameters
