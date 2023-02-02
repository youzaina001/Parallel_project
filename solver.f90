module solver

  use MPI
  use parameters
  use matrix
  use charge

  implicit none

contains

  subroutine GC(alpha,beta,gamma,b,x,iBeg,iEnd,dx,dy,me,Np)

    integer, intent(In) :: iBeg, iEnd, me, Np
    real(PR), dimension(iBeg:iEnd), intent(in) :: b
    real(PR), dimension(iBeg:iEnd), intent(inout) :: x
    real(PR), intent(in) :: alpha, beta, gamma, dx, dy
    real(PR), dimension(iBeg:iEnd) :: r, p, p1, z
    real(PR) :: alpha_GC, beta_GC, gamma_GC
    integer :: k
    integer :: statinfo
    integer, dimension(MPI_STATUS_SIZE) :: status

    x = 1._PR
    r = b - matrix_vector_product(alpha,beta,gamma,x,iBeg,iEnd,dx,dy,me,Np)
    p = r
    beta_GC = sqrt(loc_dot_product(r,r,iBeg,iEnd))
    k = 0

    do while ((beta_GC > eps) .and. (k .le. kmax))

       z = matrix_vector_product(alpha,beta,gamma,p,iBeg,iEnd,dx,dy,me,Np)
       alpha_GC = loc_dot_product(r,r,iBeg,iEnd)/loc_dot_product(z,p,iBeg,iEnd)
       x = x + alpha_GC*p
       p1 = r
       r = r - alpha_GC*z
       gamma_GC = loc_dot_product(r,r,iBeg,iEnd)/loc_dot_product(p1,p1,iBeg,iEnd)
       p = r + gamma_GC*p
       beta_GC = sqrt(loc_dot_product(r,r,iBeg,iEnd))
       k = k + 1

    end do

  end subroutine GC

  subroutine BiCGSTAB(alpha,beta,gamma,b,x,iBeg,iEnd,dx,dy,me,Np)

      integer, intent(in) :: iBeg, iEnd, me, Np
      real(PR), dimension(iBeg:iEnd), intent(in) :: b
      real(PR), dimension(iBeg:iEnd), intent(inout) :: x
      real(PR), intent(in) :: alpha, beta, gamma, dx, dy
      real(PR), dimension(iBeg:iEnd) :: r, p, v, rs, s, t
      real(PR) :: beta_STB, alpha_STB, gamma_STB, psr, somme, somme2
      real(PR) :: rho, omega, norm_r, norm_b, rho_prev
      integer :: k
      integer :: Statinfo
      integer, dimension(MPI_STATUS_SIZE) :: status

      k = 0
      x = 0._PR
      r = b - matrix_vector_product(alpha,beta,gamma,x,iBeg,iEnd,dx,dy,me,Np)
      rs = r
      rho = 1._PR
      alpha_STB = 1._PR
      omega = 1._PR
      v = 0._PR
      p = 0._PR

      norm_r = dsqrt(dot_product(r,r))                !
      norm_b = dsqrt(dot_product(b,b))

      do while((norm_r .gt. eps*norm_b) .and. (k <= kmax))

          rho_prev = rho
          rho = dot_product(rs,r)
          beta_STB = (rho/rho_prev)*(alpha_STB/omega)
          p = r + beta_STB*(p-omega*v)
          v = matrix_vector_product(alpha,beta,gamma,p,iBeg,iEnd,dx,dy,me,Np)
          alpha_STB = rho/dot_product(rs,v)
          s = r - alpha_STB*v
          t = matrix_vector_product(alpha,beta,gamma,s,iBeg,iEnd,dx,dy,me,Np)
          omega = dot_product(t,s)/dot_product(t,t)
          x = x + alpha_STB*p + omega*s
          r = s - omega*t
          norm_r = dsqrt(dot_product(r,r))
          norm_b = dsqrt(dot_product(b,b))

          k = k+1

      end do

    end subroutine BiCGSTAB

end module solver
