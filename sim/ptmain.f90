program diff
! test diffusion in a 3d lattice

use sysconfig
use particle

implicit none

integer :: i, j, nt, ntItl, ntTotal, run, runTotal, size
integer :: prtclTotal
real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax
real(b8) :: r
real(b8) :: dr(3), rsim(3,2)
real(b8), allocatable :: prtclArray(:,:)

integer :: cellTotal
real(b8) :: meanCount, varCount, avg, var
integer,  allocatable :: countArray(:), edgeList(:)
real(b8), allocatable :: cellArray(:,:,:), cellPolar(:,:), concentration(:,:,:), runCx(:,:)
real(b8), allocatable :: timePolar(:,:), timeCount(:)


! open/create file
open(unit=10, file='prtclLocation.dat', action='write', status='replace')
write(10,*) ' '
close(10)

! set total number of cell in system
cellTotal = 1
! set total number of runs
runTotal = 100
ntTotal  = 1
! initialize simulation space size
do i = 1, 3
    rsim(i,1) = 0.0_b8
    rsim(i,2) = 1.0_b8
end do
! initialize particle movement step size
dr(1) = (rsim(1,2) - rsim(1,1)) / 50.0_b8
dr(2) = (rsim(2,2) - rsim(2,1)) / 50.0_b8
dr(3) = (rsim(3,2) - rsim(3,1)) / 50.0_b8
! set time-steps needed for sytem to reach equilibrium
ntItl = 100 * int( (rsim(1,2)-rsim(1,1))**2 / dr(1)**2 )
ntItl = 10000
! set total possible number of particles in system
prtclTotal = 10000

! allocate memory
allocate( prtclArray( prtclTotal, 4))
allocate( cellArray( cellTotal, 3, 2))
allocate( cellPolar( cellTotal, 3))
allocate( timePolar( 3, ntTotal))
allocate( countArray( cellTotal))
allocate( edgeList( cellTotal))
allocate( timeCount( ntTotal))

! initialize concentration array
size = 10
allocate( concentration( size, size, size))
allocate(runCx(runTotal,size))
concentration = 0.0_b8
runCx = 0.0_b8

write(*,*) ' Run Total =', runTotal
write(*,*) 'Cell Total =', cellTotal
write(*,*) 'prtclTotal =', prtclTotal
write(*,*) '  Ninitial =', prtclTotal/2
! write(*,*) '#/timestep =', (prtclTotal/2) / 100 ! the number of particles added per timestep
write(*,*) '#/timestep =', 5 ! the number of particles added per timestep
do i = 1, 3
    write(*,*) 'dr =', dr(i), 'rsim =', rsim(i,:)
enddo
write(*,*) 'ntTotal =', ntTotal, 'ntItl =', ntItl
write(*,*)

call init_random_seed()


do run = 1, runTotal
    ! write(*,*) ' run', run

    timeCount(:)     = 0.0_b8
    countArray(:)    = 0
    cellPolar(:,:)   = 0.0_b8
    cellArray(:,:,:) = 0.0_b8
    ! initialize cell position
    ! call itlClusterSys( cellTotal, cellArray, rsim)
    call itlCellCluster( cellTotal, cellArray, rsim)
    ! call wrtOutClusterSys( cellTotal, cellArray, rsim)

    !!!  EC polarization  !!!
    !!! find which cells are on the cluster edge !!!
    ! edgeList = 0
    ! call clusterEdgeList( cellTotal, cellArray, rsim, edgeList)

    ! initialize particle positions
    nt = 1
    prtclArray(:,:) = 0.0_b8
    do i = 1, (prtclTotal/2)
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

        ! INSTANTANEOUS
        ! if ( nt == ntTotal+300 ) then
        !     call cellCount( cellTotal, prtclTotal, cellArray, prtclArray, countArray)
        !     write(100,*) countArray(1), run
        !     exit
        ! end if

        ! LONG TIME: count the particles within a cell
        ! call cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
        ! ! call cellpolarECNonAdpt( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
        ! do j = 1, 3
        !     timePolar(j,nt) = sum(cellPolar(:,j))
        ! enddo
    enddo
    ! output particle locations
    call wrtPrtclLocation( prtclTotal, run, prtclArray)
    call concentrationUpdate( prtclTotal, prtclArray, size, concentration)
    ! call concentrationXprj( size, concentration, runCx(run,:))
    do i = 1, size
        runCx(run,i) = concentration( i, 2, 3)
    enddo

    ! calculate and output time averaged molecule count
    ! do j = 1, 3
    !     cellPolar(1,j) = sum( timePolar(j,:)) / float(ntTotal)
    ! enddo
    ! ! write(200,*) cellPolar(1,:), run
    ! write(200,"(E16.8)", advance="no") cellPolar(1,1)
    ! write(200,"(E17.8)", advance="no") cellPolar(1,2)
    ! write(200,"(E17.8)", advance="no") cellPolar(1,3)
    ! write(200,"(I7)", advance="no")    run
    ! write(200,*) ''

enddo

! write concentration x projection
do i = 1, size
    avg = sum(runCx(:,i)) / float(runTotal)
    var = 0.0_b8
    do j = 1, runTotal
        var = var + (((runCx(j,i)-avg)**2.0)/float(runTotal))
    enddo
    write(300,"(E16.8)", advance="no") avg
    write(300,"(E16.8)", advance="no") var
    write(300,"(E12.4)", advance="no") float(i) * minval(rsim(:,2)) / float(size) - 0.050_b8
    write(300,*) ''
enddo

do i = 1, size
    do j = 1, runTotal
        write(310,'(E17.8)', advance='no') runCx(j,i)
    enddo
    write(310,*) ''
enddo

deallocate( prtclArray)
deallocate( cellArray)
deallocate( cellPolar)
deallocate( timePolar)
deallocate( countArray)
deallocate( timeCount)

contains


    subroutine concentrationUpdate( prtclTotal, prtclArray, size, concentration)
        implicit none
        integer,  intent(in)  :: prtclTotal, size
        real(b8), intent(in)  :: prtclArray(:,:)
        real(b8), intent(inout) :: concentration(:,:,:)
        real(b8) :: dc, voxl
        integer  :: i, j, k, n
        concentration = 0.0_b8
        ! set length of voxels
        voxl = minval(rsim(:,2)) / float(size)
        dc   = 1.0_b8

        do n = 1, prtclTotal
            if ( prtclArray(n,4) == 1.0_b8 ) then
                i = floor( prtclArray(n,1) / voxl) + 1
                j = floor( prtclArray(n,2) / voxl) + 1
                k = floor( prtclArray(n,3) / voxl) + 1
                concentration(i,j,k) = concentration(i,j,k) + dc
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


    ! calculate EC cell polarization
    ! individual cell polarization vectors are NOT adaptive
    subroutine cellpolarECNonAdpt( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal, edgeList(:)
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: clstrCOM(3), cellCOM(3,2), center(3), check, q(3)
        real(b8) :: nCell, nCOM
        integer :: i, j, k

        cellPolar(:,:) = 0.0_b8
        do i = 1, cellTotal
            nCell = 0.0_b8
            ! check if cell is on the edge of the cluster
            if ( edgeList(i) == 1 ) then
                q = 0.0_b8
                do j = 1, 3
                    center(j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
                    q(j) = center(j) - clstrCOM(j)
                enddo
                q = q / sqrt(dot_product(q,q))
                ! calculate the concentration in cell i
                do j = 1, prtclTotal
                    if ( prtclArray(j,4) == 1.0_b8 ) then
                        check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                        if ( check < 0.0 ) then
                            check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                            if ( check < 0.0 ) then
                                check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                                if ( check < 0.0 ) then
                                    nCell = nCell + 1.0_b8
                                endif
                            endif
                        endif
                    endif
                enddo
                ! write(*,*) '   cell', i, 'q =', q, 'nCell =', nCell
                cellPolar(i,:) = nCell * q
            end if
        enddo
    end subroutine cellpolarECNonAdpt


    ! calculate EC cell polarization
    subroutine cellpolarEC( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal, edgeList(:)
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: clstrCOM(3), cellCOM(3,2), center(3), check, q(3)
        real(b8) :: nCell, nCOM
        integer :: i, j, k

        cellPolar(:,:) = 0.0_b8
        ! calculate Cluster COM
        clstrCOM = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                clstrCOM(j) = clstrCOM(j) + cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
            enddo
        enddo
        clstrCOM(:) = clstrCOM(:) / float(cellTotal)
        do i = 1,3
            cellCOM(i,1)  = clstrCOM(i) - (cellArray(1,i,2) - cellArray(1,i,1)) / (2.0_b8)
            cellCOM(i,2)  = clstrCOM(i) + (cellArray(1,i,2) - cellArray(1,i,1)) / (2.0_b8)
        enddo
        write(*,*) 'cluster COM', clstrCOM
        write(*,*) '  cluster x',cellCOM(1,:), '  cluster y',cellCOM(2,:),'  cluster z',cellCOM(3,:)
        ! calculate concentration at COM
        nCOM = 0.0_b8
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                check = (prtclArray(i,1)-cellCOM(1,1))*(prtclArray(i,1)-cellCOM(1,2))
                if ( check < 0.0 ) then
                    check = (prtclArray(i,2)-cellCOM(2,1))*(prtclArray(i,2)-cellCOM(2,2))
                    if ( check < 0.0 ) then
                        check = (prtclArray(i,3)-cellCOM(3,1))*(prtclArray(i,3)-cellCOM(3,2))
                        if ( check < 0.0 ) then
                            nCOM = nCOM + 1.0_b8
                        end if
                    end if
                end if
            end if
        enddo

        cellPolar(:,:) = 0.0_b8
        do i = 1, cellTotal
            nCell = 0.0_b8
            ! check if cell is on the edge of the cluster
            if ( edgeList(i) == 1 ) then
                q = 0.0_b8
                do j = 1, 3
                    center(j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
                    q(j) = center(j) - clstrCOM(j)
                enddo
                q = q / sqrt(dot_product(q,q))
                ! calculate the concentration in cell i
                do j = 1, prtclTotal
                    if ( prtclArray(j,4) == 1.0_b8 ) then
                        check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                        if ( check < 0.0 ) then
                            check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                            if ( check < 0.0 ) then
                                check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                                if ( check < 0.0 ) then
                                    nCell = nCell + 1.0_b8
                                endif
                            endif
                        endif
                    endif
                enddo
                write(*,*) '   cell', i, 'q =', q, 'nCell =', nCell, 'nCOM =', nCOM
                cellPolar(i,:) = q * (nCell - nCOM)
            end if
        enddo
    end subroutine cellpolarEC


    ! calculate MW cell polarization
    subroutine cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: center(3), check, q(3), qtot(3)
        integer :: i, j, k
        cellPolar(:,:) = 0.0_b8
        qtot = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                center(j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
            enddo
            qtot = 0.0_b8
            do j = 1, prtclTotal
                q = 0.0_b8
                if ( prtclArray(j,4) == 1.0_b8 ) then
                    check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                    if ( check < 0.0 ) then
                        check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                        if ( check < 0.0 ) then
                            check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                            if ( check < 0.0 ) then
                                do k = 1, 3
                                    q(k) = prtclArray(j,k) - center(k)
                                enddo
                                q = q / sqrt(dot_product(q,q))
                                ! q = q
                            endif
                        endif
                    endif
                endif
                qtot = qtot + q
            enddo
            ! cellPolar(i,:) = qtot / sqrt(dot_product(qtot,qtot))
            cellPolar(i,:) = qtot
        enddo

    end subroutine cellpolarMW


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


    ! output total cluster polarization to fort.141
    subroutine wrtPlrTotal( nRun, cellTotal, cellPolar, nt)
        implicit none
        integer,  intent(in) :: cellTotal, nRun, nt
        real(b8), intent(in), dimension(:,:) :: cellPolar
        real(b8) :: px, py, pz
        integer  :: i

        px = 0.0_b8
        py = 0.0_b8
        pz = 0.0_b8
        do i = 1, cellTotal
            px = px + cellPolar(i,1)
            py = py + cellPolar(i,2)
            pz = pz + cellPolar(i,3)
        enddo
        write(100+nRun,"(E16.8)", advance="no") px
        write(100+nRun,"(E17.8)", advance="no") py
        write(100+nRun,"(E17.8)", advance="no") pz
        write(100+nRun,"(I10)", advance="no")   nt
        write(100+nRun,*) ''
    end subroutine wrtPlrTotal

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
