module particle
! module for subroutines/functions related to particle tracking

use sysconfig

contains

    ! add flux of particles at top x-boundary
    ! gradient is assumed to be in the x-direction (1)
    subroutine prtclFlux( q, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: q, dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        real(b8) :: r
        integer  :: j, k
        ! add flux to top x boundary
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
                write(*,*) 'FLUX UPDATE: hit max number of particles'
                write(*,*) '             consider increasing "prtclTotal"'
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
        do k = 1, 3
            if ( prtclArray(j,k) == rsim(k,1) .OR. prtclArray(j,k) == rsim(k,2) ) then
                write(*,*) 'FLUX: particle', j, 'on boundary', prtclArray(j,1:3)
            end if
        enddo

    end subroutine prtclFlux


    ! move particles and check boundary conditions
    ! x=rsim(1,1) boundary absorbing; x=rsim(1,2) boundary reflective
    ! all other boundaries periodic
    subroutine prtclUpdate( p, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: p, dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, j, dim
        real(b8) :: r

        do i = 1, prtclTotal
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! with probability 'p', move a distance dr(j) in random direction
                call random_number(r)
                if ( r < (p*1.0_b8/6.0_b8) ) then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                elseif( r < (p*2.0_b8/6.0_b8) )then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                elseif( r < (p*3.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                elseif( r < (p*4.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                elseif( r < (p*5.0_b8/6.0_b8) )then
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                elseif( r <= p )then
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) + dr(3)
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
                elseif( dim == 1 )then
                    ! absorbing boundary at x = 0, reflective at x = rsim(1,2)
                    if( prtclArray(i,dim) < rsim(dim,1) )then
                        prtclArray(i,4) = 0.0_b8
                    elseif( prtclArray(i,dim) > rsim(dim,2) )then
                        prtclArray(i,dim) = rsim(dim,2) - (prtclArray(i,dim) - rsim(dim,2))
                    endif
                endif
            endif
        enddo
    end subroutine prtclUpdate


    ! move particles and check boundary conditions
    ! x=rsim(1,1) boundary absorbing; x=rsim(1,2) boundary reflective
    ! all other boundaries reflective
    subroutine prtclUpdateReflect( p, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: p, dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, j, dim
        real(b8) :: r

        do i = 1, prtclTotal
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! with probability 'p', move a distance dr(j) in random direction
                call random_number(r)
                if ( r < (p*1.0_b8/6.0_b8) ) then
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                    ! check for absorbing boundary
                    if ( prtclArray(i,1) <= rsim(1,1) ) then
                        prtclArray(i,4) = 0.0_b8
                    end if
                elseif( r < (p*2.0_b8/6.0_b8) )then
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                    ! check for reflective boundary
                    if ( prtclArray(i,1) > rsim(1,2) ) then
                        prtclArray(i,1) = prtclArray(i,1) - (prtclArray(i,1) - rsim(1,2))
                    end if
                elseif( r < (p*3.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                    ! check for reflective boundary
                    if ( prtclArray(i,2) < rsim(2,1) ) then
                        prtclArray(i,2) = prtclArray(i,2) + (rsim(2,1) - prtclArray(i,2))
                    end if
                elseif( r < (p*4.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                    ! check for reflective boundary
                    if ( prtclArray(i,2) > rsim(1,2) ) then
                        prtclArray(i,2) = prtclArray(i,2) - (prtclArray(i,2) - rsim(2,2))
                    end if
                elseif( r < (p*5.0_b8/6.0_b8) )then
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                    ! check for reflective boundary
                    if ( prtclArray(i,3) < rsim(3,1) ) then
                        prtclArray(i,3) = prtclArray(i,3) + (rsim(3,1) - prtclArray(i,3))
                    end if
                elseif( r <= p )then
                    prtclArray(i,3) = prtclArray(i,3) + dr(3)
                    ! check for reflective boundary
                    if ( prtclArray(i,3) > rsim(3,2) ) then
                        prtclArray(i,3) = prtclArray(i,3) - (prtclArray(i,3) - rsim(3,2))
                    end if
                endif

                ! do j = 1, 3
                !     if ( prtclArray(i,j) == rsim(j,1) .OR. prtclArray(i,j) == rsim(j,2) ) then
                !         write(*,*) 'particle', i, 'on boundary', prtclArray(i,1:4)
                !     end if
                ! enddo
            endif
        enddo
    end subroutine prtclUpdateReflect


    ! move particles and check boundary conditions
    ! all boundaries are periodic
    subroutine prtclUpdateAllPeriodic( p, dr, rsim, prtclArray)
        implicit none
        real(b8), intent(in)    :: p ,dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        real(b8) :: r
        integer  :: i, j

        do i = 1, prtclTotal
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! move the particle by distance dr(j) in random direction
                call random_number(r)
                if ( r < (p/6.0_b8) ) then
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                elseif( r < (p*2.0_b8/6.0_b8) )then
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                elseif( r < (p*3.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                elseif( r < (p*4.0_b8/6.0_b8) )then
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                elseif( r < (p*5.0_b8/6.0_b8) )then
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                elseif( r <= p )then
                    prtclArray(i,3) = prtclArray(i,3) + dr(3)
                endif

                ! check periodic boundary conditions
                do j = 1, 3
                    if( prtclArray(i,j) < rsim(j,1) )then
                        prtclArray(i,j) = rsim(j,2) - (rsim(j,1) - prtclArray(i,j))
                        exit
                    elseif( prtclArray(i,j) > rsim(j,2) )then
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
