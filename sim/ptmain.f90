program diff

use sysconfig
use particle
use polarization

implicit none

integer :: i, j, nGeo, nt, ntItl, ntTotal, run, cSize(3)
real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax
real(b8) :: dtReal, p, q, r
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
call getSysLengthScales( dr, rsim)
! set event probabilities and time-steps needed for sytem to reach equilibrium
call getProbTimeScale( ntItl, dtReal, p, q)
ntTotal = ntItl
write(*,*) 'particle track'
write(*,*) 'p =', p, 'q =', q, ' dtReal =', dtReal
write(*,*) 'lReal =', lReal, 'dReal =', dReal
write(*,*)
do i = 1, 3
    write(*,*) i, 'dr =', dr(i), 'rsim =', rsim(i,:)
enddo
write(*,*)
write(*,*) 'ntTotal =', ntTotal, 'ntItl =', ntItl
write(*,*)

! allocate memory
allocate( prtclArray( prtclTotal, 4))
allocate( cellArray( cellTotal, 3, 2))
allocate( cellPolar( cellTotal, 3))
allocate( timePolar( 3, ntTotal))
allocate( edgeList( cellTotal))
allocate( timeCount( ntTotal))
! initialize concentration array
write(*,*) '  cSize:'
do i = 1, 3
    cSize(i) = ceiling( (rsim(i,2) - rsim(i,1)) / (10.0*dr(i)))
    write(*,*) '  ', cSize(i)
enddo
allocate( concentration( cSize(1), cSize(2), cSize(3)))
allocate( runCx(runTotal,cSize(1)))

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
    call itlCellCluster( cellTotal, cellArray, rsim)
    ! call itl2DClusterNN( cellArray, rsim)
    edgeList = 0
    call clusterEdgeList( cellTotal, cellArray, rsim, edgeList)

    do run = 1, runTotal
        write(*,*) ' run', run

        timeCount(:)     = 0.0_b8
        cellPolar(:,:)   = 0.0_b8
        timePolar(:,:)   = 0.0_b8

        ! initialize particle positions
        prtclArray(:,:) = 0.0_b8

        ! let system reach equilibrium
        do nt = 1, ntItl
            ! call prtclUpdate( p, dr, rsim, prtclArray)
            call prtclUpdateReflect( p, dr, rsim, prtclArray)
            ! add flux of particles
            call prtclFlux( q, dr, rsim, prtclArray)
        enddo

        ! gather statistics
        do nt = 1, ntTotal
            ! update particle location and check boundary conditions
            ! call prtclUpdate( p, dr, rsim, prtclArray)
            call prtclUpdateReflect( p, dr, rsim, prtclArray)
            ! add flux of particles
            call prtclFlux( q, dr, rsim, prtclArray)

            ! call cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
            ! call cellpolar2DMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
            ! call cellpolarECNonAdpt( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
            ! do j = 1, 3
            !     timePolar(j,nt) = sum(cellPolar(:,j))
            ! enddo
        enddo
        call concentrationUpdate( prtclTotal, prtclArray, cSize, concentration)
        call concentrationXprj( cSize(1), concentration, runCx(run,:))
        ! do i = 1, size
        !     runCx(run,i) = concentration( i, 2, 3)
        ! enddo

        ! calculate and output time averaged molecule count
        ! do j = 1, 3
        !     cellPolar(1,j) = sum( timePolar(j,:)) / float(ntTotal)
        ! enddo
        ! call wrtPlrTotal( nGeo, run, cellPolar)
        ! call wrtCellLocation( cellArray)
    enddo

    ! call wrtCellLocation( cellArray)

    ! write concentration x projection
    call wrtConcentrationX( cSize, runCx, rsim)
    close(12)
enddo

! write out simulation information / parameters
write(*,*)
write(*,*) ' geo Total =', geoTotal
write(*,*) ' Run Total =', runTotal
write(*,*) 'Cell Total =', cellTotal
write(*,*) 'prtclTotal =', prtclTotal
do i = 1, 3
    write(*,*) i, 'dr =', dr(i), 'rsim =', rsim(i,:)
enddo
write(*,*) 'ntTotal =', ntTotal, 'ntItl =', ntItl
write(*,*)
write(*,*) 'Length and Time Conversions:'
write(*,*) '  1 sim.   length =', rReal / rCell, 'microns'
write(*,*) '  1 sim. timestep =', dtReal, 'seconds'
write(*,*)
write(*,*) 'Average gradient:'
write(*,*) kReal / (dReal * syReal * szReal), '/ microns^4'
write(*,*) kReal / (dReal * syReal * szReal) * (rReal / rCell)**4.0, '/(sim. length)^4'
write(*,*)


deallocate( prtclArray)
deallocate( cellArray)
deallocate( cellPolar)
deallocate( timePolar)
deallocate( timeCount)

contains

    subroutine concentrationUpdate( prtclTotal, prtclArray, cSize, concentration)
        implicit none
        integer,  intent(in)  :: prtclTotal, cSize(3)
        real(b8), intent(in)  :: prtclArray(:,:)
        real(b8), intent(inout) :: concentration(:,:,:)
        real(b8) :: dc, il, jl, kl, voxl
        integer  :: i, j, k, n, ic(3), lc(3)
        concentration = 0.0_b8
        ! set length of voxels
        do i = 1, 3
            lc(i) = rsim(i,2) / float(cSize(i))
        enddo
        voxl = lc(1) * lc(2) * lc(3)
        dc   = 1.0_b8

        do n = 1, prtclTotal
            if ( prtclArray(n,4) == 1.0_b8 ) then
                do i = 1, 3
                    if ( prtclArray(n,i) == rsim(i,1) ) then
                        ic(i) = 1
                    elseif( prtclArray(n,i) == rsim(i,2) )then
                        ic(i) = cSize(i)
                    else
                        ic(i) = floor( prtclArray(n,i) / lc(i)) + 1
                    end if
                    if ( ic(i) > cSize(i) ) then
                        write(*,*) 'prtcl', n, prtclArray(n,1:4)
                    end if
                enddo
                concentration(ic(1),ic(2),ic(3)) = concentration(ic(1),ic(2),ic(3)) + (dc / voxl)
            end if
        enddo
    end subroutine concentrationUpdate


    ! calculate x projection of concentration by averaging over y, z dimensions
    subroutine concentrationXprj( cSize, concentration, xConcentration)
        implicit none
        integer, intent(in)   :: cSize(3)
        real(b8), intent(in)  :: concentration(:,:,:)
        real(b8), intent(out) :: xConcentration(:)
        real(b8) :: mC(cSize(1),cSize(2))
        integer :: i, j ,k
        mC = 0.0_b8
        ! average over z dimension
        do i = 1, cSize(1)
            do j = 1, cSize(2)
                mC(i,j) = sum(concentration(i,j,:)) / float(cSize(3))
            enddo
        enddo
        ! average over y dimension
        do i = 1, cSize(1)
            xConcentration(i) = sum(mC(i,:)) / float(cSize(2))
        enddo
    end subroutine concentrationXprj


    ! write concentration x projection
    subroutine wrtConcentrationX( cSize, runCx, rsim)
        implicit none
        integer,  intent(in) :: cSize(3)
        real(b8), intent(in) :: runCx(:,:), rsim(3,2)
        integer  :: i
        real(b8) :: avg, var
        do i = 1, cSize(1)
            avg = sum(runCx(:,i)) / float(runTotal)
            var = 0.0_b8
            do j = 1, runTotal
                var = var + (((runCx(j,i)-avg)**2.0)/float(runTotal))
            enddo
            write(300,"(E16.8)", advance="no") avg
            write(300,"(E16.8)", advance="no") var
            write(300,"(E12.4)", advance="no") float(i)*minval(rsim(:,2))/float(cSize(1))-(minval(rsim(:,2))/float(cSize(1))/2.0_b8)
            write(300,*) ''
        enddo
    end subroutine wrtConcentrationX


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
