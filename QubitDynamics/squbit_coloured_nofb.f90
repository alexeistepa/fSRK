!--------
!squbit_coloured_nofb
!--------
!For thesis. Symmetric qubit with coloured noise and no feedback (and no extra dephasing). 
!
!--------

subroutine qubit_drf(n, np, t, q, Pr, f)
    implicit none

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), dimension(n), intent(in) :: q
    real(8), intent(in) :: t
    real(8), dimension(np), intent(in) :: Pr
    real(8), dimension(n), intent(out) :: f
    
    real(8) :: y, z, c
    real(8) :: lam, delI, gam, eps_w
    real(8) :: f_y, f_z, f_c
    
    ! Co-ordinate labels
    y = q(1)
    z = q(2)
    c = q(3)
    
    ! Parameter labels
    lam = Pr(1)
    delI = Pr(2)
    gam = Pr(3)
    eps_w = Pr(4)
    
    ! Drift function components
    f_y = - lam*z - (1 - eps_w)*(delI**2)*y*(z**2) - (0.5*delI**2)*y*eps_w - delI*y*z*c
    f_z = lam*y + (1 - eps_w)*(delI**2)*(z - z**3) - delI*(z + 1)*(z - 1)*c
    f_c = - gam*c
    
    ! Output
    f = (/ f_y, f_z, f_c /)
    
end subroutine

subroutine qubit_dif(n, np, t, q, Pr,g)
    implicit none

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), dimension(n), intent(in) :: q
    real(8), intent(in) :: t
    real(8), dimension(np), intent(in) :: Pr
    real(8), dimension(n), intent(out) :: g
    
    real(8) :: y, z, c
    real(8) :: lam, delI, gam, eps_w
    real(8) :: g_y, g_z, g_c
    
    ! Co-ordinate labels
    y = q(1)
    z = q(2)
    c = q(3)
    
    ! Parameter labels
    lam = Pr(1)
    delI = Pr(2)
    gam = Pr(3)
    eps_w = Pr(4)
    
    ! Diffusion function components
    g_y = - delI*y*z*eps_w
    g_z = - delI*(z + 1)*(z - 1)*eps_w 
    g_c = gam*(1-eps_w)
    
    ! Output
    g = (/ g_y, g_z, g_c /)
    
end subroutine


!------------------------------------------------------------------
! fSRK
!------------------------------------------------------------------
!2nd Order Weak Stochastic Runge-Kutta algorithm for numerically
!approximating the solution of stochastic ordinary differential
!equations in Ito form, dx = f(x,t)*dt + g(x,t)*dW, where dW is a
!scalar Weiner noise process.Designed to be imported into python
!with the wrapper f2py.
!------------------------------------------------------------------
!Author: Alexei Stepanenko
!-------------------------
!----------
!References
!----------
!1. Breuer and Petruccione - The Theory of Open Quantum Systems pg 362
!eq (7.47) - eq. (7.49).
!
!2. Kloeden and Platen - Numerical solutions of stochastic differential
!equations pg 486 eq. (1.1). This book contains generalisations of this
!algorithm to multiple Weiner processes.
!
!3. Milstein and  Tretyakov - Stochastic Numerics for Mathematical
!Physics pg 104 eq. 2.20 .
!-----------------------------------------------------------------------



subroutine drift(n, np, t, r, P, f)
    implicit none

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), intent(in) :: t
    real(8), dimension(n), intent(in) :: r
    real(8), dimension(np), intent(in) :: P
    real(8), dimension(n), intent(out) :: f
    
    call qubit_drf(n, np, t, r, P, f) ! EDIT HERE
    
end subroutine 

subroutine diffus(n, np, t, r, P, g)
    implicit none

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), intent(in) :: t
    real(8), dimension(n), intent(in) :: r
    real(8), dimension(np), intent(in) :: P
    real(8), dimension(n), intent(out) :: g
    
    call qubit_dif(n, np, t, r, P, g) ! EDIT HERE
    
end subroutine 

subroutine fSRKstp(n, np, r0, t0, delt, nstp, P, rout, tout)
    implicit none
    
!Evolves an initial state forward delt in time using nstp steps of    
!the 2nd order weak Runge-Kutta

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), dimension(n), intent(in) :: r0 ! Initial state 
    real(8), intent(in) :: t0 ! Initial time 
    real(8), intent(in) :: delt ! Time increment
    integer, intent(in) :: nstp
    real(8), dimension(np), intent(in) :: P ! Parameter array
    real(8), dimension(n), intent(out) :: rout ! Final state
    real(8), intent(out) :: tout  ! Final time
       
    
    real(8), dimension(n) :: rn ! Current state
    real(8) :: tn ! Current time
    real(8), dimension(n) :: rn_b, rn_p, rn_m ! Eq. (7.48) and (7.49) (Ref. 1)
    real(8), dimension(n) :: term1, term2, term3 ! Terms of eq. (7.47) (Ref.1)
    real(8) :: xi, dW ! Current Weiner proccess and increment
    real(8) :: dt ! Integration time step
    real(8) :: sdt ! sqrt(dt)
    integer :: i ! Dummy
    real(8), dimension(n) :: f_rn_, g_rn_, f_rnb_, g_rnp_, g_rnm_ ! 'Evaluation holders' (for speed)
    real(8) :: normal
    
    dt = delt/dble(nstp)
    sdt = sqrt(dt)
    
    ! Main loop
    rn = r0
    tn = t0
    do i = 1 , nstp
        ! Weiner 
        xi = normal()
        dW = xi*sdt    
        
        ! Intermediate states
        call drift(n, np, tn, rn, P, f_rn_)
        call diffus(n, np, tn, rn, P, g_rn_)
        rn_b = rn + f_rn_*dt + g_rn_*dW
        rn_p = rn + f_rn_*dt + g_rn_*sdt
        rn_m = rn + f_rn_*dt - g_rn_*sdt
        
        ! Terms 
        call drift(n, np, tn + dt, rn_b, P, f_rnb_)
        call diffus(n, np, tn + dt, rn_p, P, g_rnp_)
        call diffus(n, np, tn + dt, rn_m, P, g_rnm_)
        term1 = 0.5_8*dt*(f_rnb_ + f_rn_ )
        term2 = 0.25_8*dW*(g_rnp_ + g_rnm_ + 2*g_rn_)
        term3 = 0.25_8*(g_rnp_ - g_rnm_)*(dW**2 - dt)/sdt
        
        ! Next state and time
        rn = rn + term1 + term2 + term3
        tn = tn + dt
        
    end do
    
    ! Final state and time
    rout = rn
    tout = tn
    
end subroutine fSRKstp

subroutine fSRK(n, np, r0, t0, tf, nout, nps, P, rvals, tvals) 
    implicit none

    integer, intent(in) :: n ! Dimension of state-space
    integer, intent(in) :: np ! Number of parameters
    real(8), dimension(n), intent(in) :: r0 ! Initial state 
    real(8), intent(in) :: t0 ! Initial time 
    real(8), intent(in) :: tf ! Final time
    integer, intent(in) :: nout ! Number of output states
    integer, intent(in) :: nps ! Number of integration steps between outputs
    real(8), dimension(np), intent(in) :: P ! Parameter array        
    real(8), dimension(n,nout), intent(out) :: rvals ! Returned time series
    real(8), dimension(nout), intent(out) :: tvals ! Returned series of time values corrensponding to above time-series
    
    real(8), dimension(n) :: rn ! Current state
    real(8) :: tn ! Current time 
    real(8) :: delt ! Time between outputs 
    integer :: i
    
    delt = (tf - t0)/dble(nout - 1) 
    rvals(:,1) = r0(:)
    tvals(1) = t0
    rn = r0
    tn = t0
    do i = 1, nout - 1
        call fSRKstp(n, np, rn, tn, delt, nps, P, rn, tn)
        rvals(:,i+1) = rn(:)
        tvals(i+1) = tn
    end do
    
end subroutine

function normal()
    implicit none
    
!Normally distributed random number generator (mean = 0,std = 1).
!Uses box Muller method.
    
    real(8) :: normal
    real(8) :: pi, u1, u2
    pi = 4.0*atan(1.0)
    
    call random_number(u1)
    call random_number(u2)
    normal = (-2.0*log(u1))**0.5*cos(2.0*pi*u2)
    
end function

!--------
!End fSRK
!--------

! f2py -c -m squbit_coloured_nofb squbit_coloured_nofb.f90
