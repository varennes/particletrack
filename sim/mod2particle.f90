module particle
! module for subroutines/functions related to particle tracking

use sysconfig

! number of particles added to simulation space
integer,  parameter :: nJ = 10

contains

    ! add flux of particles at top x-boundary
    ! gradient is assumed to be in the x-direction (1)
    subroutine prtclFlux( q, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: q, dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        real(b8) :: r
        integer  :: j, k
        ! add flux to top boundary
        ! first check probability of event occuring
        if ( q < 1.0 ) then
            call random_number(r)
            if ( q < r ) then
                ! write(*,*) ' q =', q, 'r =', r, 'no particle addition'
                return
            end if
        end if
        ! write(*,*) ' q =', q, 'r =', r, 'particle addition'
        j = 1
        do while( prtclArray(j,4) == 1.0_b8 )
            j = j + 1
            if ( j > prtclTotal ) then
                exit
            end if
        enddo
        if ( j > prtclTotal ) then
            return
        end if
        prtclArray(j,4) = 1.0_b8
        call random_number(r)
        prtclArray(j,1) = rsim(1,2) - r*dr(1)
        do k = 2, 3
            call random_number(r)
            prtclArray(j,k) = r *(rsim(k,2) - rsim(k,1)) + rsim(k,1)
        enddo
    end subroutine prtclFlux


    ! move particles and check boundary conditions
    ! x-boundaries are absorbing
    subroutine prtclUpdate( p, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: p, dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, j, dim
        real(b8) :: r

        ! first check probability of event occuring
        if ( p < 1.0 ) then
            call random_number(r)
            if ( p < r ) then
                ! write(*,*) ' q =', q, 'r =', r, 'no particle addition'
                return
            end if
        end if

        do i = 1, prtclTotal
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! move the particle by distance dr(j) in random direction
                call random_number(r)
                if ( r < (1.0_b8/6.0_b8) ) then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                elseif( r < (2.0_b8/6.0_b8) )then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                elseif( r < (3.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                elseif( r < (4.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                elseif( r < (5.0_b8/6.0_b8) )then
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                else
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) + dr(1)
                endif

                ! check boundary conditions
                if ( dim /= 1 ) then
                    ! periodic boundaries perpendicular to gradient
                    do j = 2, 3
                        if( prtclArray(i,j) < rsim(j,1) )then
                            ! write(*,*) 'too small', prtclArray(i,j), i, j
                            prtclArray(i,j) = rsim(j,2) - (rsim(j,1) - prtclArray(i,j))
                            exit
                        elseif( prtclArray(i,j) > rsim(j,2) )then
                            ! write(*,*) 'too big', prtclArray(i,j), i, j
                            prtclArray(i,j) = (prtclArray(i,j) - rsim(j,2)) + rsim(j,1)
                            exit
                        endif
                    enddo
                else
                    ! absorbing boundary parallel to gradient
                    if( (prtclArray(i,dim)<rsim(dim,1)) .OR. (prtclArray(i,dim)>rsim(dim,2)) )then
                        prtclArray(i,4) = 0.0_b8
                    endif
                endif
            endif
        enddo
    end subroutine prtclUpdate


    ! move particles and check boundary conditions
    ! all boundaries are periodic
    subroutine prtclUpdateAllPeriodic( N, dr, rsim, prtclArray)
        implicit none
        integer,  intent(in)    :: N
        real(b8), intent(in)    :: dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, j
        real(b8) :: r

        do i = 1, N
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! move the particle by distance dr(j) in random direction
                call random_number(r)
                if ( r < (1.0_b8/6.0_b8) ) then
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                elseif( r < (2.0_b8/6.0_b8) )then
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                elseif( r < (3.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                elseif( r < (4.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                elseif( r < (5.0_b8/6.0_b8) )then
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                else
                    prtclArray(i,3) = prtclArray(i,3) + dr(1)
                endif

                ! check periodic boundary conditions
                do j = 1, 3
                    if( prtclArray(i,j) < rsim(j,1) )then
                        ! write(*,*) 'too small', prtclArray(i,j), i, j
                        prtclArray(i,j) = rsim(j,2) - (rsim(j,1) - prtclArray(i,j))
                        exit
                    elseif( prtclArray(i,j) > rsim(j,2) )then
                        ! write(*,*) 'too big', prtclArray(i,j), i, j
                        prtclArray(i,j) = (prtclArray(i,j) - rsim(j,2)) + rsim(j,1)
                        exit
                    endif
                enddo
            endif
        enddo
    end subroutine prtclUpdateAllPeriodic


    ! output particle location data
    subroutine wrtPrtclLocation( N, nt, prtclArray)
        implicit none
        integer,  intent(in) :: N, nt
        real(b8), intent(in) :: prtclArray(:,:)
        integer :: i
        open(unit=11, file='prtclLocation.dat', action='write', status='old')
        do i = 1, N
            if ( prtclArray(i,4) == 1.0_b8 ) then
                write(11,*) prtclArray(i,1:3), nt
            end if
        enddo
    end subroutine wrtPrtclLocation

end module
