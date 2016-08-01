program diff
! test diffusion in a 3d lattice

implicit none

integer, parameter:: b8 = selected_real_kind(14)

integer :: i, j, N, Ntotal, nc

real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax
real(b8) :: r
real(b8) :: dr(3), rsim(3,2)
real(b8), allocatable :: prtclArray(:,:)

! set total possible number of particles in system
prtclTotal = 20
! initialize simulation space size
do i = 1, 3
    rsim(i,1) = 0.0_b8
    rsim(i,2) = 1.0_b8
end do
! initialize particle movement step size
dr(1) = (rsim(1,2) - rsim(1,1)) / 5.0_b8
dr(2) = (rsim(2,2) - rsim(2,1)) / 5.0_b8
dr(3) = (rsim(3,2) - rsim(3,1)) / 5.0_b8
! initialize track array size
allocate( prtclArray( prtclTotal, 4))

write(*,*) 'Ntotal =',Ntotal
write(*,*) '     N =',N
do i = 1, 3
    write(*,*) 'dr =', dr(1)
enddo
write(*,*)

call init_random_seed()

! initialize particle positions positions
prtclArray(:,:) = 0.0_b8
do i = 1, (prtclTotal/2)
    prtclArray(i,4) = 1.0_b8
    do k = 1, 3
        call random_number(r) ! in cpm code use ran1() function
        prtclArray(i,j) = r *(rsim(j,2) - rsim(j,1)) + rsim(j,1)
    enddo
    write(*,*) 'i=', i, prtclArray(i,:)
enddo

! update particle location and check boundary conditions
call prtclUpdate( prtclTotal, dr, rsim, prtclArray)
write(*,*) '  move  '
do i = 1, prtclTotal
    if ( prtclArray(i,4) == 1.0_b8 ) then
        write(*,*) 'i=', i, prtclArray(i,:)
    end if
enddo

! add flux of particles
call prtclFlux( prtclTotal, dr, rsim, prtclArray)


contains
    ! gradient is assumed to be in the x-direction (1)
    ! add flux of particles at x-boundaries
    subroutine prtclFlux( N, dr, rsim, prtclArray)
        implicit none
        integer,  intent(in)    :: N
        real(b8), intent(in)    :: dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        real(b8) :: a, d, g
        integer  :: i, j, k, nJ
        ! cross-sectional area of flux
        a = (rsim(1,1)-rsim(1,2)) * (rsim(3,1)-rsim(3,2))
        ! set diffusion and gradient
        d = 4.00_b8
        g = 0.50_b8
        ! calculate flux
        nJ = int( d*g*a)
        write(*,*) 'flux =', nJ
        do i = 1, nJ
            j = 1
            do while( prtclArray(j,4) == 1.0_b8 )
                j = j + 1
            enddo
            prtclArray(j,1) = rsim(1,1) + dr(1)
            do k = 2, 3
                call random_number(r) ! in cpm code use ran1() function
                prtclArray(j,k) = r *(rsim(k,2) - rsim(k,1)) + rsim(k,1)
            enddo
        enddo
    end subroutine prtclFlux


    ! move particles and check boundary conditions
    subroutine prtclUpdate( N, dr, rsim, prtclArray)
        implicit none
        integer,  intent(in)    :: N
        real(b8), intent(in)    :: dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i
        real(b8) :: r

        do i = 1, N
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! move the particle by distance dr(j) in random direction
                call random_number(r)
                if ( r < (1.0_b8/6.0_b8) ) then
                    dim = 1
                    prtclArray(1,i) = prtclArray(i,1) - dr(1)
                elseif( r < (2.0_b8/6.0_b8) )then
                    dim = 1
                    prtclArray(1,i) = prtclArray(i,1) + dr(1)
                elseif( r < (3.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(2,i) = prtclArray(i,2) - dr(2)
                elseif( r < (4.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(2,i) = prtclArray(i,2) + dr(2)
                elseif( r < (5.0_b8/6.0_b8) )then
                    dim = 3
                    prtclArray(3,i) = prtclArray(i,3) - dr(3)
                else
                    dim = 3
                    prtclArray(3,i) = prtclArray(i,3) + dr(1)
                endif

                ! check boundary conditions
                if ( dim /= 1 ) then
                    ! reflective boundaries perpendicular to gradient
                    do j = 2, 3
                        if( prtclArray(i,j) < rsim(j,1) )then
                            write(*,*) 'too small', prtclArray(i,j), i, j
                            prtclArray(i,j) = rsim(j,2) - (rsim(j,1) - prtclArray(i,j))
                            exit
                        elseif( prtclArray(i,j) > rsim(j,2) )then
                            write(*,*) 'too big', prtclArray(i,j), i, j
                            prtclArray(i,j) = (prtclArray(i,j) - rsim(j,2)) + rsim(j,1)
                            exit
                        endif
                    enddo
                else
                    ! absorbing boundary parallel to gradient
                    if( (prtclArray(i,dim)<rsim(dim,1)) .OR. (prtclArray(i,dim)>rsim(dim,2)) )then
                        write(*,*) 'absorbed', prtclArray(i,1:3), i
                        prtclArray(i,4) = 0.0_b8
                    endif
                endif
            endif
        enddo
    end subroutine prtclUpdate

end program


! initialize RANDOM_SEED
subroutine init_random_seed()
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
end subroutine init_random_seed


! returns random number between 0 - 1
! real(b8) function ran1()
!     implicit none
!     real(b8) x
!     call random_number(x) ! built in fortran 90 random number function
!     ran1=x
! end function ran1
