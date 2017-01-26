program diff

use sysconfig
use particle
use polarization

implicit none

integer :: i, j, nGeo, nt, ntItl, ntTotal, run, cSize(3), overflow = 0
real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax, tCPU0, tCPU1
real(b8) :: dtReal, p, q, r
real(b8) :: dr(3), rsim(3,2), clstrCOM(3)
real(b8), allocatable :: prtclArray(:,:), prtclItl(:,:)

real(b8) :: meanCount, varCount, avg, var
integer,  allocatable :: edgeList(:)
real(b8), allocatable :: cellArray(:,:,:), cellPolar(:,:), concentration(:,:,:), runCx(:,:)
real(b8), allocatable :: timePolar(:,:)
real(b8), allocatable :: cellCenter(:,:)

real(b8) :: sMW, sEC

character(len=1024) :: filename

call cpu_time(tCPU0)

! initialize simulation space size
call getSysLengthScales( dr, rsim)
! set event probabilities and time-steps needed for sytem to reach equilibrium
call getProbTimeScale( ntItl, dtReal, p, q)
ntTotal = ntItl

write(*,*) 'particle track'
write(*,*) 'ntItl =', ntItl, ', in seconds:', float(ntItl) * dtReal
write(*,*) 'ntTotal =', ntTotal, ', in seconds:', float(ntTotal) * dtReal
write(*,*) 'p =', p, 'q =', q, ' dtReal =', dtReal
write(*,*)
write(*,*) 'Average gradient:'
write(*,*) kReal / (dReal * syReal * szReal), '/ microns^4'
write(*,*) 'Average concentration:'
write(*,*) lReal * kReal / (dReal * syReal * szReal * 2.0_b8), '/ microns^3'
write(*,*)

! allocate memory
allocate( prtclArray( prtclTotal, 4))
allocate( prtclItl( prtclTotal, 4))
allocate( cellArray( cellTotal, 3, 2))
allocate( cellPolar( cellTotal, 3))
allocate( timePolar( 3, ntTotal))
allocate( edgeList( cellTotal))
allocate( cellCenter( cellTotal, 3) )
! initialize concentration array
do i = 1, 3
    cSize(i) = 2 * ceiling( (rsim(i,2) - rsim(i,1)) / (10.0*dr(i)))
enddo
write(*,*) '  Concentration grid size:', cSize(:)
write(*,*) '  '
allocate( concentration( cSize(1), cSize(2), cSize(3)))
allocate( runCx(runTotal,cSize(1)))

call init_random_seed()

! initialize particles and let system reach equilibrium
prtclArray(:,:) = 0.0_b8
prtclItl(:,:)   = 0.0_b8

! open output data file
write (filename, "(A3,I0.3,A4)") 'sij', cellTotal, '.dat'
open( 12, file=filename)

do nGeo = 1, geoTotal
    write(*,"(A8,I3)") '  nGeo =', nGeo
    ! set cell configuration
    cellArray(:,:,:) = 0.0_b8
    call itl1DChain( cellArray, rsim)
    ! call itl3DRandom( cellTotal, cellArray, rsim)
    ! call itl3DClusterNN( cellArray, rsim)
    ! call itl2DClusterNN( cellArray, rsim)
    ! call itl2DRandom( cellTotal, cellArray, rsim)
    call getCellCenter( cellArray, cellCenter)
    call clusterCenter( cellArray, clstrCOM)

    edgeList(:) = 0
    ! call clusterEdgeList( cellTotal, cellArray, rsim, edgeList)
    ! call EdgeList2D( cellArray, edgeList)
    call EdgeList1D( cellArray, edgeList)

    call wrtOutClusterSys( cellTotal, cellArray, rsim)

    call SijMW( cellCenter, clstrCOM, sMW)
    call SijEC( cellCenter, clstrCOM, edgeList, sEC)

    write(12,"(E16.8)", advance="no") sMW
    write(12,"(E16.8)", advance="no") sEC
    write(12,"(I7)", advance="no") nGeo
    write(12,*) ''

enddo

close(12)

! write out simulation information / parameters
call cpu_time(tCPU1)
write(*,*)
write(*,*) ' geo Total =', geoTotal
write(*,*) ' Run Total =', runTotal
write(*,*) 'Cell Total =', cellTotal
write(*,*) 'prtclTotal =', prtclTotal
write(*,*) '   ntTotal =', ntTotal, 'ntItl =', ntItl
write(*,*)
write(*,*) 'CPU Run Time (min)', (tCPU1 - tCPU0) / 60.0

deallocate( prtclArray)
deallocate( prtclItl)
deallocate( cellArray)
deallocate( cellPolar)
deallocate( timePolar)
deallocate( cellCenter)

contains

    subroutine concentrationUpdate( prtclTotal, prtclArray, cSize, concentration)
        implicit none
        integer,  intent(in)  :: prtclTotal, cSize(3)
        real(b8), intent(in)  :: prtclArray(:,:)
        real(b8), intent(inout) :: concentration(:,:,:)
        real(b8) :: dc, lc(3), voxl
        integer  :: i, j, k, n, ic(3)
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
                    if ( abs(ic(i)) > cSize(i) ) then
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
    subroutine wrtConcentrationX( cSize, runCx, rsim, N)
        implicit none
        integer,  intent(in) :: cSize(3), N
        real(b8), intent(in) :: runCx(:,:), rsim(3,2)
        integer  :: i
        real(b8) :: avg, var
        do i = 1, cSize(1)
            avg = sum(runCx(:,i)) / float(N)
            var = 0.0_b8
            do j = 1, N
                var = var + (((runCx(j,i)-avg)**2.0)/float(N))
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
