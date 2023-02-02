Module parameters

  implicit none

  integer, parameter :: PR = 8 ! Double précision
  integer, parameter :: Nx = 100, Ny = 100, BC_choice = 1
  integer, parameter :: kmax = 1000, nmax = 10, OVERLAP = 2, nmax_scwz = 100
  integer, parameter :: ar = 0, br = 1 
  real(PR), parameter :: pi = 4._PR*atan(1._PR)
  real(PR), parameter :: Lx = 1.0, Ly = 1.0, D = 1.0, dt = 1.0, eps = 1E-6, eps_scwz = 1E-6
  real(PR), parameter :: dx = Lx/(Nx+1), dy = Ly/(Ny+1)
  real(PR), parameter :: alpha = 1._PR + 2._PR*D*dt*(1._PR/(dx**2)+1._PR/(dy**2))
  real(PR), parameter :: beta = -D*dt/(dx**2)
  real(PR), parameter :: gamma = -D*dt/(dy**2)

end module parameters
