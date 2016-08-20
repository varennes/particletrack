module parameters

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)

!!! REAL UNIT PARAMETERS !!!
real(b8), parameter ::  aReal =     10.0_b8   ! in microns
real(b8), parameter ::  bReal =      5.0_b8   ! in microns
real(b8), parameter ::  lReal =    50.0_b8   ! in microns
real(b8), parameter :: syReal =    50.0_b8   ! in microns
real(b8), parameter :: szReal =    50.0_b8   ! in microns
real(b8), parameter ::  dReal =    10.0_b8   ! in microns^2/s
real(b8), parameter ::  kReal =    100.0_b8   ! in 1/s

end module

!!! PARAMETER DESCRIPTIONS !!!
!
!  aReal = cell radius
!  bReal = particle jump distance
!  lReal = system length in gradient direction
! syReal = system length perpendicular to gradient
! szReal = system length perpendicular to gradient
!      s = syReal * szReal; the area of the surface where particles are produced
!  dReal = diffusion coefficient
!  kReal = particle production rate
!
!!! PARAMETER DESCRIPTIONS !!!
