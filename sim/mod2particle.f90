module particle
! module for subroutines/functions related to particle tracking

use parameters
use sysconfig

contains

    ! add flux of particles at top x-boundary
    ! gradient is assumed to be in the x-direction (1)
    subroutine prtclFlux( q, rsim, prtclArray, overflow)
        implicit none
        real(b8), intent(in)    :: q, rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer,  intent(inout) :: overflow
        real(b8) :: r
        integer  :: j, k
        ! first check probability of event occuring
        if ( q < 1.0 ) then
            call random_number(r)
            if ( q < r ) then
                return
            end if
        end if
        j = 1
        do while( prtclArray(j,4) == 1.0_b8 )
            j = j + 1
            if ( j > prtclTotal ) then
                if ( overflow == 0 ) then
                    write(*,*) 'FLUX UPDATE: hit max number of particles'
                    write(*,*) '             consider increasing "prtclTotal"'
                    overflow = 1
                end if
                exit
            end if
        enddo
        if ( j > prtclTotal ) then
            return
        end if
        ! add a particle at a random location in y-z plane at the x-boundary
        prtclArray(j,4) = 1.0_b8
        prtclArray(j,1) = rsim(1,2)
        do k = 2, 3
            call random_number(r)
            prtclArray(j,k) = r *(rsim(k,2) - rsim(k,1)) + rsim(k,1)
        enddo

    end subroutine prtclFlux


    ! move particles and check boundary conditions
    ! x=rsim(1,1) boundary absorbing; x=rsim(1,2) boundary reflective
    ! all other boundaries periodic
    subroutine prtclUpdate( p, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: p, rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, j, dim
        real(b8) :: r

        do i = 1, prtclTotal
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! with probability 'p', move a distance dr(j) in random direction
                call random_number(r)
                if ( r <= p ) then
                    ! check what direction the particle moves
                    if ( r <= p / 6.0_b8 ) then
                        prtclArray(i,1) = prtclArray(i,1) + bReal
                    elseif ( r <= p / 3.0_b8 ) then
                        prtclArray(i,1) = prtclArray(i,1) - bReal
                    elseif ( r <= p / 2.0_b8 ) then
                        prtclArray(i,2) = prtclArray(i,2) + bReal
                    elseif ( r <= 2.0_b8 * p / 3.0_b8 ) then
                        prtclArray(i,2) = prtclArray(i,2) - bReal
                    elseif ( r <= 5.0_b8 * p / 6.0_b8 ) then
                        prtclArray(i,3) = prtclArray(i,3) + bReal
                    elseif ( r <= p ) then
                        prtclArray(i,3) = prtclArray(i,3) - bReal
                    end if

                    ! check boundary conditions
                    ! periodic boundaries perpendicular to x
                    ! x = 0 is an absorbing boundary, x = L is reflective
                    do j = 1, 3
                        if ( prtclArray(i,j) > rsim(j,2) ) then
                            if ( j == 1 ) then
                                prtclArray(i,j) = rsim(j,2) - ( prtclArray(i,j) - rsim(j,2) )
                            else
                                prtclArray(i,j) = rsim(j,1) + ( prtclArray(i,j) - rsim(j,2) )
                            end if
                        elseif ( prtclArray(i,j) < rsim(j,1) ) then
                            if ( j == 1 ) then
                                prtclArray(i,4) = 0.0_b8
                            else
                                prtclArray(i,j) = rsim(j,2) - ( rsim(j,1) - prtclArray(i,j) )
                            end if
                        end if
                    enddo
                endif
            endif
        enddo
    end subroutine prtclUpdate


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
