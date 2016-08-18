module parameters

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)

!!! REAL UNIT PARAMETERS !!!
real(b8), parameter ::  aReal =    10.0_b8   ! in microns
real(b8), parameter ::  bReal =     1.0_b8   ! in microns
real(b8), parameter ::  lReal =   100.0_b8   ! in microns
real(b8), parameter :: syReal =    10.0_b8   ! in microns
real(b8), parameter :: szReal =    10.0_b8   ! in microns
real(b8), parameter ::  dReal =   100.0_b8   ! in microns^2/s
real(b8), parameter ::  kReal =    10.0_b8   ! in 1/s

contains

end module
