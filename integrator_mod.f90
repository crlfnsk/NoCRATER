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
    Real(dp)    :: da
    Integer(i8) :: kk
    
    Write(outputname, '("euler", E0.1)') real(theseIC%resolution)

    Call print_header(trim(outputname))
    Call print_values(MM%time, MM%a, radionuclide%concGalactic)

    kk = 0

    do
        MM%dt = MM%Period/theseIC%resolution
        radionuclide%delta_concGalactic = radionuclide%prodGalactic * MM%dt
        da = -2 * theseIC%B_const/MM%a * MM%dt/(1e6*yr)

        MM%time = MM%time + MM%dt
        MM%a = MM%a + da
        radionuclide%delta_concGalactic = dfunc(MM%dt, radionuclide%prodGalactic, Belambda, radionuclide%concGalactic ) * MM%dt
        radionuclide%concGalactic = radionuclide%concGalactic + radionuclide%delta_concGalactic

        kk = kk + 1

        if (MM%a < 1.0_dp) then
            print*, "Semi major axis already smaller than 1 au."
            Exit
        Else if (abs(MM%a - 1.0_dp) < epsilon) then
            Exit
        Else if (modulo(kk, Int(theseIC%resolution/theseIC%kmax)) == 0) then
            call print_values(MM%time/(1e6*yr), MM%a, radionuclide%concGalactic)
        End if
    end do

    call print_values(MM%time/(1e6*yr), MM%a, radionuclide%concGalactic)
    Close(10)
    
end subroutine euler

! Hier kommt Runge Kutta hin
    
end module integrator_mod