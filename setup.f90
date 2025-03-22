!!+setup.f90
!!
!! Bis nicht aufgesetzt. Wenn es wirklich skalierbar gemacht werden sollte, benötigte es noch sehr viele Änderungen. 
!!
!! module: 
!!       
!!-
module setup
    use Constants
    implicit none

    ! Initial conditions:
    TYPE :: IC
    !    character(len=2)    :: solarmodulation
        real(dp)            :: sunmasses
        integer(i8)         :: resolution
        integer(i8)         :: kmax
    !    character(len=2)    :: Nrev !Anzahl von Revolutions als Abbruchkriterium an?
    !    integer(i8)         :: NR ! nach wie vielen revolutionen?
        real(dp)            :: diameter
        real(dp)            :: composition
        Real(dp)            :: B_const
        real(dp)            :: sma
        real(dp)            :: ecc
        real(dp)            :: thetadeg
        REAL(dp)            :: pre_conc ! , dimension(:), ALLOCATABLE
    END TYPE IC

    TYPE :: Isotope
    !    REAL(dp)                :: ONEau_prodGalactic
    !    REAL(dp)                :: ONEau_prodSolar
        REAL(dp)                :: prodGalactic
    !    REAL(dp), dimension(3)  :: prodSolar
        REAL(dp)                :: delta_concGalactic
    !    REAL(dp), dimension(3)  :: Delta_concSolar
        REAL(dp)                :: concGalactic
    !    REAL(dp), dimension(3)  :: concSolar
    END TYPE Isotope

    TYPE :: Particle
        REAL(dp)                :: time
        REAL(dp)                :: a
        REAL(dp)                :: e
        REAL(dp)                :: r
    !    REAL(dp)                :: rp
    !    REAL(dp)                :: ra
        REAL(dp)                :: Period
        REAL(dp)                :: dt
        INTEGER(i8)             :: N_Rev
    END TYPE Particle

contains
    



end module setup

