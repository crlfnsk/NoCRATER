module isotope_io
    use constants
    use setup

    implicit none
    
contains
    
subroutine print_header(filename)
    implicit none
    character(len=274), intent(in)    :: filename 

    Open(10, File=trim(filename)//'.txt')
    Write(10, '("time [Ma]", 1X, "a [au]", 1X, "NBe [10^9/g]")' )
    Write(10, '("--------", 1X, "--------", 1X, "---------------")')

end subroutine print_header

subroutine print_values(time,  sma, concentration)
    REAl(dp), intent(in) :: time,  sma, concentration

    Write (10,'(F8.3, 1X, F8.3, 1X, F15.9)') time, sma, concentration*1.0e-9
    
end subroutine print_values

end module isotope_io