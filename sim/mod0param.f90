module parameters

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)

!!! REAL UNIT PARAMETERS !!!
real(b8), parameter ::  rReal =     5.0_b8   ! in microns
real(b8), parameter ::  bReal =     0.5_b8   ! in microns
real(b8), parameter ::  lReal =    20.0_b8   ! in microns
real(b8), parameter :: syReal =    20.0_b8   ! in microns
real(b8), parameter :: szReal =    20.0_b8   ! in microns
real(b8), parameter ::  dReal =    20.0_b8   ! in microns^2/s
real(b8), parameter ::  kReal =    50.0_b8   ! in 1/s

real(b8), parameter ::  hReal = rReal/2.0_b8 ! in microns

end module

!!! PARAMETER DESCRIPTIONS !!!
!
!  All length parameters are in units of microns.
!  All time parameters are in units of seconds.
!
!  rReal = cell radius
!  bReal = particle jump distance
!  lReal = system length in gradient direction
! syReal = system length perpendicular to gradient
! szReal = system length perpendicular to gradient
!      s = syReal * szReal; the area of the surface where particles are produced
!  dReal = diffusion coefficient
!  kReal = particle production rate
!
!  hReal = cell height, only used for 2D configuration. Cells are discs with
!          hReal and radius rReal
!
!!! PARAMETER DESCRIPTIONS !!!
