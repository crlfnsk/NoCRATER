!! Ich möchte herausfinden welche Konzentration von einem Isotop heruaskommt,
!! wenn ein Meteoroid durch den PR-Effekt zur Erde spiralt und dabei 
!! durch kosmisch bestrahlt wird.
!!
!! 1.0 
!! 2.1. Für den Anfang Produktionsrate = const.
!! 2.2. Produktionsrate ist von 1/r^2 abhängig.
!! 2.3. GCR + SCR.

program integrator
    use constants
    use setup
    use isotope_io
    use integrator_mod
    implicit none

    TYPE(IC) :: example
    Type(Particle) :: MiMi
    TYPE(Isotope) :: Be_euler

    
!   N'(t) = Produktionsrate(t) - const. * N(t)
!   N(0) = 0 
!   Produktionsrate = const.
!   

!   Analytische Lösung: N(t) = P/const. * (1 - e^const./t )
!   
!   Euler Ansatz: 
!   yˆ_i(x_o + h) = y_i(x_o) + h · y′(x) |x=x_o = y_i(x_o) + h · f_i(x_o, y_o)
    
!   y_i+1 = y_i + h * y_i'(y0, x0)

    ! Initial Condition: 

    example%sunmasses = 1.0_dp
    example%resolution = 1.0e2_dp
    example%kmax = 10000_dp
    example%diameter = 350_dp
    example%composition = 2.5_dp
    example%B_const = 3*L_Sun/(8*pi * c**2 * example%diameter*1e-6 * example%composition*1000) ! const for elliptical orbits in m2/s
    example%sma = 2.5_dp
    example%thetadeg = 0.0_dp
    example%pre_conc = 0.0_dp

    ! Setup particle

    MiMi%time = 0
    MiMi%a = example%sma * AU ! in m
    MiMi%Period = two_pi*sqrt(MiMi%a**3/(G*example%sunmasses*M_Sun)) ! in s
    MiMi%dt = MiMi%Period/example%resolution
    MiMi%N_Rev = 0

    print*, "Period = ", MiMi%Period/(yr), "Ma"
!    REAL(dp)                :: delta_concGalactic

    Be_euler%prodGalactic = PGCRrBe
    Be_euler%concGalactic = example%pre_conc

    Call euler(example, MiMi, Be_euler)

end program integrator