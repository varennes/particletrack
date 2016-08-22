program diff

use sysconfig
use particle
use polarization

implicit none

integer :: i, j, nGeo, nt, ntItl, run, size
real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax
real(b8) :: r
real(b8) :: dr(3), rsim(3,2)
real(b8), allocatable :: prtclArray(:,:)

real(b8) :: meanCount, varCount, avg, var
integer,  allocatable :: edgeList(:)
real(b8), allocatable :: cellArray(:,:,:), cellPolar(:,:), concentration(:,:,:), runCx(:,:)
real(b8), allocatable :: timePolar(:,:), timeCount(:)

character(len=1024) :: filename

! open/create file
open(unit=10, file='prtclLocation.dat', action='write', status='replace')
write(10,*) ' '
close(10)

! initialize simulation space size
do i = 1, 3
    rsim(i,1) = 0.0_b8
    rsim(i,2) = 1.0_b8
end do
! initialize particle movement step size
dr(1) = (rsim(1,2) - rsim(1,1)) / 10.0_b8
dr(2) = (rsim(2,2) - rsim(2,1)) / 10.0_b8
dr(3) = (rsim(3,2) - rsim(3,1)) / 10.0_b8
! set time-steps needed for sytem to reach equilibrium
ntItl = 10 * int( (rsim(1,2)-rsim(1,1))**2 / dr(1)**2 )

! allocate memory
allocate( prtclArray( prtclTotal, 4))
allocate( cellArray( cellTotal, 3, 2))
allocate( cellPolar( cellTotal, 3))
allocate( timePolar( 3, ntTotal))
allocate( edgeList( cellTotal))
allocate( timeCount( ntTotal))
! initialize concentration array
size = 10
allocate( concentration( size, size, size))
allocate(runCx(runTotal,size))


call itl3DClusterNN( cellArray, rsim)
call wrtCellLocation( cellArray)
call init_random_seed()

do nGeo = 1, geoTotal
    write(*,"(A8,I3)") '  nGeo =', nGeo
    ! open output data file
    write (filename, "(A1,I0.3,A4)") 's', nGeo, '.200'
    open( 12, file=filename)

    ! initialize arrays
    concentration = 0.0_b8
    runCx = 0.0_b8
    ! set cell configuration
    cellArray(:,:,:) = 0.0_b8
    ! call itl2DCellCluster( cellTotal, cellArray, rsim)
    call itlCellCluster( cellTotal, cellArray, rsim)
    edgeList = 0
    call clusterEdgeList( cellTotal, cellArray, rsim, edgeList)

    do run = 1, runTotal
        write(*,*) ' run', run

        timeCount(:)     = 0.0_b8
        cellPolar(:,:)   = 0.0_b8
        timePolar(:,:)   = 0.0_b8

        ! initialize particle positions
        nt = 1
        prtclArray(:,:) = 0.0_b8
        do i = 1, (prtclTotal/10)
            prtclArray(i,4) = 1.0_b8
            do j = 1, 3
                call random_number(r) ! in cpm code use ran1() function
                prtclArray(i,j) = r *(rsim(j,2) - rsim(j,1)) + rsim(j,1)
            enddo
        enddo
        ! call wrtPrtclLocation( prtclTotal, nt, prtclArray)
        ! call concentrationUpdate( prtclTotal, prtclArray, concentration)

        ! let system reach equilibrium
        do nt = 1, ntItl
            call prtclUpdate( prtclTotal, dr, rsim, prtclArray)
            ! add flux of particles
            call prtclFlux( prtclTotal, dr, rsim, prtclArray)
        enddo

        ! gather statistics
        do nt = 1, ntTotal
            ! update particle location and check boundary conditions
            call prtclUpdate( prtclTotal, dr, rsim, prtclArray)
            ! add flux of particles
            call prtclFlux( prtclTotal, dr, rsim, prtclArray)

            ! LONG TIME: count the particles within a cell
            call cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
            ! call cellpolar2DMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
            ! call cellpolarECNonAdpt( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
            do j = 1, 3
                timePolar(j,nt) = sum(cellPolar(:,j))
            enddo
        enddo
        ! output particle locations
        ! call wrtPrtclLocation( prtclTotal, run, prtclArray)
        ! call concentrationUpdate( prtclTotal, prtclArray, size, concentration)
        ! ! call concentrationXprj( size, concentration, runCx(run,:))
        ! do i = 1, size
        !     runCx(run,i) = concentration( i, 2, 3)
        ! enddo

        ! calculate and output time averaged molecule count
        do j = 1, 3
            cellPolar(1,j) = sum( timePolar(j,:)) / float(ntTotal)
        enddo
        call wrtPlrTotal( nGeo, run, cellPolar)
        ! call wrtCellLocation( cellArray)
    enddo

    ! call wrtCellLocation( cellArray)

    ! write concentration x projection
    ! do i = 1, size
    !     avg = sum(runCx(:,i)) / float(runTotal)
    !     var = 0.0_b8
    !     do j = 1, runTotal
    !         var = var + (((runCx(j,i)-avg)**2.0)/float(runTotal))
    !     enddo
    !     write(300,"(E16.8)", advance="no") avg
    !     write(300,"(E16.8)", advance="no") var
    !     write(300,"(E12.4)", advance="no") float(i) * minval(rsim(:,2)) / float(size) - (minval(rsim(:,2)) / float(size) / 2.0_b8)
    !     write(300,*) ''
    ! enddo
    !
    ! do i = 1, size
    !     do j = 1, runTotal
    !         write(310,'(E17.8)', advance='no') runCx(j,i)
    !     enddo
    !     write(310,*) ''
    ! enddo
    close(12)
enddo

! write out simulation information / parameters
write(*,*)
write(*,*)
write(*,*) ' geo Total =', geoTotal
write(*,*) ' Run Total =', runTotal
write(*,*) 'Cell Total =', cellTotal
write(*,*) 'prtclTotal =', prtclTotal
write(*,*) 'prtcl intl =', prtclTotal/2
write(*,*) '#/timestep =', nJ ! the number of particles added per timestep
do i = 1, 3
    write(*,*) i, 'dr =', dr(i), 'rsim =', rsim(i,:)
enddo
write(*,*) 'ntTotal =', ntTotal, 'ntItl =', ntItl

deallocate( prtclArray)
deallocate( cellArray)
deallocate( cellPolar)
deallocate( timePolar)
deallocate( timeCount)

contains

    subroutine concentrationUpdate( prtclTotal, prtclArray, size, concentration)
        implicit none
        integer,  intent(in)  :: prtclTotal, size
        real(b8), intent(in)  :: prtclArray(:,:)
        real(b8), intent(inout) :: concentration(:,:,:)
        real(b8) :: dc, il, jl, kl, voxl
        integer  :: i, j, k, n
        concentration = 0.0_b8
        ! set length of voxels
        il = rsim(1,2) / float(size)
        jl = rsim(2,2) / float(size)
        kl = rsim(3,2) / float(size)
        voxl = il * jl * kl
        dc   = 1.0_b8

        do n = 1, prtclTotal
            if ( prtclArray(n,4) == 1.0_b8 ) then
                i = floor( prtclArray(n,1) / il) + 1
                j = floor( prtclArray(n,2) / jl) + 1
                k = floor( prtclArray(n,3) / kl) + 1
                concentration(i,j,k) = concentration(i,j,k) + (dc / voxl)
            end if
        enddo
    end subroutine concentrationUpdate


    ! calculate x projection of concentration by averaging over y, z dimensions
    subroutine concentrationXprj( size, concentration, xConcentration)
        implicit none
        integer, intent(in)   :: size
        real(b8), intent(in)  :: concentration(:,:,:)
        real(b8), intent(out) :: xConcentration(:)
        real(b8) :: mC(size,size)
        integer :: i, j ,k
        mC = 0.0_b8
        ! average over z dimension
        do i = 1, size
            do j = 1, size
                mC(i,j) = sum(concentration(i,j,:)) / float(size)
            enddo
        enddo
        ! average over y dimension
        do i = 1, size
            xConcentration(i) = sum(mC(i,:)) / float(size)
        enddo
    end subroutine concentrationXprj


    ! count the number of particels within cells volume
    subroutine cellCount( cellTotal, prtclTotal, cellArray, prtclArray, countArray)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        integer,  intent(out) :: countArray(:)
        integer  :: i, j
        real(b8) :: check
        do i = 1, cellTotal
            countArray(i) = 0
            do j = 1, prtclTotal
                if ( prtclArray(j,4) == 1.0_b8 ) then
                    check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                    if ( check < 0.0 ) then
                        ! write(*,*) 'passed x-check:', prtclArray(j,1:3)
                        check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                        if ( check < 0.0 ) then
                            ! write(*,*) 'passed y-check:', prtclArray(j,1:3)
                            check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                            if ( check < 0.0 ) then
                                ! write(*,*) 'passed z-check:', prtclArray(j,1:3)
                                countArray(i) = countArray(i) + 1
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    end subroutine cellCount

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
