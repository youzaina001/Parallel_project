module charge

  use parameters

  implicit none

contains

  subroutine charge_proc(n,Np,me,iBeg,iEnd)

    integer, intent(in)  :: n, Np, me
    integer, intent(out) :: iBeg, iEnd
    integer              :: r, q

    q = n/Np
    r = modulo(n,Np)

    if (me < r) then

      iBeg = me * (q + 1) + 1
      iEnd = (1 + me) * (q + 1)

    else

      iBeg = 1 + r + me*q
      iEnd = iBeg + q - 1

    end if

    if (me .ne. 0) then

      iBeg = iBeg - OVERLAP

    end if

    if (me .ne. (Np-1)) then

      iEnd = iEnd + OVERLAP

    end if

  end subroutine charge_proc

end module charge
