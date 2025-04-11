module integrator_mod
    use constants
    use setup
    use isotope_io
    implicit none
    
contains

!Define here the function that is supposed to be analyzed as f. 
FUNCTION dfunc(dt, prodrate, RNlambda, N_0) RESULT(N_dash)
    IMPLICIT NONE
    REAL(KIND=8) :: dt, prodrate, RNlambda, N_dash, N_0
    
    N_dash = prodrate * dt - RNlambda * N_0

END FUNCTION dfunc

subroutine euler(theseIC,  MM, radionuclide)
    TYPE(IC), intent(in)            ::  theseIC
    TYPE(Particle), intent(inout)   :: MM
    TYPE(Isotope), intent(inout)    :: radionuclide

    character(len=49) :: outputname
    Real(dp)    :: epsilon = 1.0e-6_dp
    Real(qp)    :: da
    Integer(i8) :: kk
    
    Write(outputname, '("euler", E0.1, ".dat")') real(theseIC%resolution)

    kk = 0

    Call print_header(trim(outputname))
    Call print_values(MM%time, MM%a/AU, radionuclide%concGalactic)

    do
        MM%dt = MM%Period/theseIC%resolution ! in s
        !radionuclide%delta_concGalactic = radionuclide%prodGalactic * MM%dt

        ! Poynting-Robertson effect auf Einheiten achten
        da = -2 * theseIC%B_const/MM%a * MM%dt ! in m
        !print*, "da = ", da

        MM%time = MM%time + MM%dt
        MM%a = MM%a + da

        !radionuclide%delta_concGalactic = dfunc(MM%dt, radionuclide%prodGalactic, Belambda, radionuclide%concGalactic ) * MM%dt
        !radionuclide%concGalactic = radionuclide%concGalactic + radionuclide%delta_concGalactic
        radionuclide%concGalactic = radionuclide%concGalactic + (radionuclide%prodGalactic - radionuclide%concGalactic*Belambda) * MM%dt

        kk = kk + 1

        if (MM%a .LT. 1.0_dp) then
            print*, "Semi major axis already smaller than 1 au."
            Exit
        Else if (abs(MM%a/AU - 1.0_dp) .LT. epsilon) then
            Exit
        Else if (kk .GT. theseIC%kmax) then
            call print_values(MM%time/(1e6*yr), MM%a/AU, radionuclide%concGalactic)
            kk = 0
        Else
            ! do nothing
        End if
    end do

    call print_values(MM%time/(1e6*yr), MM%a/AU, radionuclide%concGalactic)
    Close(10)
    
end subroutine euler

! Hier kommt Runge Kutta hin
    
end module integrator_mod