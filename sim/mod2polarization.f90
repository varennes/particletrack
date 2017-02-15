module polarization
! module for subroutines/functions related to cell polarization

use sysconfig

contains

    ! calculate EC cell polarization, individual polarization vectors are NOT adaptive
    subroutine polar3DECnonadpt( cellArray, clstrCOM, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: edgeList(:)
        real(b8), intent(in)  :: cellArray(:,:,:), clstrCOM(3), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: nCell(cellTotal), rCell(cellTotal,3), r(3)
        integer :: i, j, n
        integer :: exitCount

        cellPolar(:,:) = 0.0_b8
        if ( cellTotal == 1 ) then
            return
        end if

        rCell = 0.0_b8
        nCell = 0.0_b8
        do i = 1, cellTotal
            if ( edgeList(i) == 1 ) then
                do j = 1, 3
                    rCell(i,j) = ((cellArray(i,j,2) + cellArray(i,j,1)) / 2.0_b8) - clstrCOM(j)
                enddo
                rCell(i,:) = rCell(i,:) / sqrt( rCell(i,1)**2 + rCell(i,2)**2 + rCell(i,3)**2)
            end if
        enddo

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell on the edge
                do n = 1, cellTotal
                    if ( edgeList(n) == 0 ) then
                        cycle
                    end if
                    if ( r(1) > cellArray(n,1,1) .AND. r(1) <= cellArray(n,1,2) ) then
                        if ( r(2) > cellArray(n,2,1) .AND. r(2) <= cellArray(n,2,2) ) then
                            if ( r(3) > cellArray(n,3,1) .AND. r(3) <= cellArray(n,3,2) ) then
                                nCell(n) = nCell(n) + 1.0_b8
                            end if
                        end if
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo

        do i = 1, cellTotal
            if ( edgeList(i) == 0 ) then
                cycle
            end if
            cellPolar(i,:) = nCell(i) * rCell(i,:)
        enddo
    end subroutine polar3DECnonadpt


    ! MW cell polarization vector
    subroutine polar3DMWv1( cellArray, prtclArray, cellPolar )
        implicit none
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: center(3), r(3), rmag
        integer  :: i, j, n, cellCheck
        integer  :: exitCount, count

        center(:)      = 0.0_b8
        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        count = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell
                do n = 1, cellTotal
                    cellCheck = 0
                    if ( r(1) > cellArray(n,1,1) .AND. r(1) <= cellArray(n,1,2) ) then
                        if ( r(2) > cellArray(n,2,1) .AND. r(2) <= cellArray(n,2,2) ) then
                            if ( r(3) > cellArray(n,3,1) .AND. r(3) <= cellArray(n,3,2) ) then
                                cellCheck = 1
                                count = count + 1
                                ! write(*,*) i, prtclArray(i,1), prtclArray(i,2), prtclArray(i,3)
                            end if
                        end if
                    end if
                    ! add q contribution and exit cell loop
                    if ( cellCheck == 1 ) then
                        do j = 1, 3
                            r(j) = prtclArray(i,j) - ((cellArray(n,j,2) + cellArray(n,j,1))/2.0_b8)
                        enddo
                        rmag = sqrt( r(1)**2 + r(2)**2 + r(3)**2 )
                        do j = 1, 3
                            cellPolar(n,j) = cellPolar(n,j) + r(j) / rmag
                        enddo
                        exit
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
        ! write(*,*) 'v1 count =',count
    end subroutine polar3DMWv1


    subroutine polar3DMWv2( cellArray, prtclArray, cellPolar )
        implicit none
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: center(3), r(3), rmag
        integer  :: i, j, n, cellCheck
        integer  :: exitCount, count

        center(:)      = 0.0_b8
        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        count = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell
                do n = 1, cellTotal
                    cellCheck = 0
                    if ( r(1) > cellArray(n,1,1) .AND. r(1) < cellArray(n,1,2) ) then
                        if ( r(2) > cellArray(n,2,1) .AND. r(2) < cellArray(n,2,2) ) then
                            if ( r(3) > cellArray(n,3,1) .AND. r(3) < cellArray(n,3,2) ) then
                                cellCheck = 1
                                count = count + 1
                                ! write(*,*) i, prtclArray(i,1), prtclArray(i,2), prtclArray(i,3)
                            end if
                        end if
                    end if
                    ! add q contribution and exit cell loop
                    if ( cellCheck == 1 ) then
                        do j = 1, 3
                            r(j) = prtclArray(i,j) - ((cellArray(n,j,2) + cellArray(n,j,1))/2.0_b8)
                        enddo
                        rmag = sqrt( r(1)**2 + r(2)**2 + r(3)**2 )
                        do j = 1, 3
                            cellPolar(n,j) = cellPolar(n,j) + r(j) / rmag
                        enddo
                        exit
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
        ! write(*,*) 'v2 count =',count
    end subroutine polar3DMWv2


    subroutine polarSphereMW( cellCenter, prtclArray, cellPolar )
        implicit none
        real(b8), intent(in)  :: cellCenter(:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: r(3), dr
        integer  :: i, j, n, exitCount

        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do n = 1, cellTotal
                    do j = 1, 3
                        r(j) = prtclArray(i,j) - cellCenter(n,j)
                    enddo
                    dr = sqrt( sum(r*r))
                    if ( dr < rReal ) then
                        ! add q contribution and exit cell loop
                        do j = 1, 3
                            cellPolar(n,j) = cellPolar(n,j) + r(j) / dr
                        enddo
                        exit
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
    end subroutine polarSphereMW


    ! calculate EC cell polarization, individual polarization vectors are NOT adaptive
    subroutine polarSphereEC( cellCenter, clstrCOM, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: edgeList(:)
        real(b8), intent(in)  :: cellCenter(:,:), clstrCOM(3), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: nCell(cellTotal), rCell(cellTotal,3), r(3), dr
        integer :: i, j, n
        integer :: exitCount

        cellPolar(:,:) = 0.0_b8
        if ( cellTotal == 1 ) then
            return
        end if

        rCell = 0.0_b8
        nCell = 0.0_b8
        do i = 1, cellTotal
            if ( edgeList(i) == 1 ) then
                do j = 1, 3
                    rCell(i,j) = cellCenter(i,j) - clstrCOM(j)
                enddo
                rCell(i,:) = rCell(i,:) / sqrt( rCell(i,1)**2 + rCell(i,2)**2 + rCell(i,3)**2)
            end if
        enddo

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell on the edge
                do n = 1, cellTotal
                    if ( edgeList(n) == 0 ) then
                        cycle
                    end if
                    do j = 1, 3
                        r(j) = prtclArray(i,j) - cellCenter(n,j)
                    enddo
                    dr = sqrt( sum(r*r))
                    if ( dr < rReal ) then
                        nCell(n) = nCell(n) + 1.0_b8
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo

        do i = 1, cellTotal
            if ( edgeList(i) == 0 ) then
                cycle
            end if
            cellPolar(i,:) = nCell(i) * rCell(i,:)
        enddo
    end subroutine polarSphereEC


    ! calculate EC cell polarization, individual polarization vectors are NOT adaptive
    subroutine polarDiscEC( cellCenter, clstrCOM, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: edgeList(:)
        real(b8), intent(in)  :: cellCenter(:,:), clstrCOM(3), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: nCell(cellTotal), rCell(cellTotal,2), r(2), dr
        integer :: i, j, n
        integer :: exitCount

        cellPolar(:,:) = 0.0_b8
        if ( cellTotal == 1 ) then
            return
        end if
        rCell = 0.0_b8
        nCell = 0.0_b8
        do i = 1, cellTotal
            if ( edgeList(i) == 1 ) then
                do j = 1, 2
                    rCell(i,j) = cellCenter(i,j) - clstrCOM(j)
                enddo
                rCell(i,:) = rCell(i,:) / sqrt( rCell(i,1)**2 + rCell(i,2)**2 )
            end if
        enddo

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                if ( abs(cellCenter(1,3)-prtclArray(i,3)) < (hReal/2.0_b8) ) then
                    do j = 1, 2
                        r(j) = prtclArray(i,j)
                    enddo
                    ! check if particle is in a cell on the edge
                    do n = 1, cellTotal
                        if ( edgeList(n) == 0 ) then
                            cycle
                        end if
                        do j = 1, 2
                            r(j) = prtclArray(i,j) - cellCenter(n,j)
                        enddo
                        dr = sqrt( sum(r*r))
                        if ( dr < rReal ) then
                            nCell(n) = nCell(n) + 1.0_b8
                        end if
                    enddo
                end if
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo

        do i = 1, cellTotal
            if ( edgeList(i) == 0 ) then
                cycle
            end if
            cellPolar(i,1) = nCell(i) * rCell(i,1)
            cellPolar(i,2) = nCell(i) * rCell(i,2)
        enddo
    end subroutine polarDiscEC


    subroutine polarDiscMW( cellCenter, prtclArray, cellPolar )
        implicit none
        real(b8), intent(in)  :: cellCenter(:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: r(2), dr
        integer  :: i, j, n, exitCount

        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                if ( abs(cellCenter(1,3)-prtclArray(i,3)) < (hReal/2.0_b8) ) then
                    do n = 1, cellTotal
                        do j = 1, 2
                            r(j) = prtclArray(i,j) - cellCenter(n,j)
                        enddo
                        dr = sqrt( sum(r*r))
                        if ( dr < rReal ) then
                            ! add q contribution and exit cell loop
                            do j = 1, 2
                                cellPolar(n,j) = cellPolar(n,j) + (r(j) / dr)
                            enddo
                            exit
                        end if
                    enddo
                end if
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
    end subroutine polarDiscMW


    ! calculate 2D MW cell polarization
    subroutine polar2DMW( cellArray, prtclArray, cellPolar)
        implicit none
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: r(3), rmag
        integer  :: i, j, n, cellCheck
        integer  :: exitCount, count

        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        count = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell
                do n = 1, cellTotal
                    cellCheck = 0
                    if ( r(1) > cellArray(n,1,1) .AND. r(1) < cellArray(n,1,2) ) then
                        if ( r(2) > cellArray(n,2,1) .AND. r(2) < cellArray(n,2,2) ) then
                            if ( r(3) > cellArray(n,3,1) .AND. r(3) < cellArray(n,3,2) ) then
                                cellCheck = 1
                                count = count + 1
                            end if
                        end if
                    end if
                    ! add q contribution and exit cell loop
                    if ( cellCheck == 1 ) then
                        do j = 1, 2
                            r(j) = prtclArray(i,j) - ((cellArray(n,j,2) + cellArray(n,j,1))/2.0_b8)
                        enddo
                        rmag = sqrt( r(1)**2 + r(2)**2 )
                        do j = 1, 2
                            cellPolar(n,j) = cellPolar(n,j) + (r(j) / rmag)
                        enddo
                        exit
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
    end subroutine polar2DMW


    ! 2D EC cell polarization, cube cells, individual polarization vectors are NOT adaptive
    subroutine polar2DEC( cellCenter, clstrCOM, edgeList, prtclArray, cellArray, cellPolar)
        implicit none
        integer,  intent(in)  :: edgeList(:)
        real(b8), intent(in)  :: cellCenter(:,:), clstrCOM(3), cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: nCell(cellTotal), rCell(cellTotal,2), r(2), dr
        integer :: i, j, n
        integer :: exitCount

        cellPolar(:,:) = 0.0_b8
        if ( cellTotal == 1 ) then
            return
        end if
        rCell = 0.0_b8
        nCell = 0.0_b8
        do i = 1, cellTotal
            if ( edgeList(i) == 1 ) then
                do j = 1, 2
                    rCell(i,j) = cellCenter(i,j) - clstrCOM(j)
                enddo
                rCell(i,:) = rCell(i,:) / sqrt( rCell(i,1)**2 + rCell(i,2)**2 )
            end if
        enddo

        exitCount = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                if ( abs(cellCenter(1,3)-prtclArray(i,3)) < (hReal/2.0_b8) ) then
                    ! check if particle is in a cell on the edge
                    do n = 1, cellTotal
                        if ( edgeList(n) == 0 ) then
                            cycle
                        end if
                        do j = 1, 2
                            r(j) = prtclArray(i,j)
                        enddo
                        if ( r(1) > cellArray(n,1,1) .AND. r(1) <= cellArray(n,1,2) ) then
                            if ( r(2) > cellArray(n,2,1) .AND. r(2) <= cellArray(n,2,2) ) then
                                nCell(n) = nCell(n) + 1.0_b8
                            end if
                        end if
                    enddo
                end if
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo

        do i = 1, cellTotal
            if ( edgeList(i) == 0 ) then
                cycle
            end if
            cellPolar(i,1) = nCell(i) * rCell(i,1)
            cellPolar(i,2) = nCell(i) * rCell(i,2)
        enddo
    end subroutine polar2DEC


    ! gradient sensing test
    subroutine gradientTest( cellArray, cellCenter, prtclArray, cellPolar)
        implicit none
        real(b8), intent(in)  :: cellArray(:,:,:), cellCenter(:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: r(3), rmag
        integer  :: i, j, n, cellCheck
        integer  :: count1, count2, exitCount

        cellPolar(:,:) = 0.0_b8

        exitCount = 0
        count1 = 0
        count2 = 0
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                exitCount = 0
                do j = 1, 3
                    r(j) = prtclArray(i,j)
                enddo
                ! check if particle is in a cell
                do n = 1, cellTotal
                    cellCheck = 0
                    if ( r(1) > cellArray(n,1,1) .AND. r(1) < cellArray(n,1,2) ) then
                        if ( r(2) > cellArray(n,2,1) .AND. r(2) < cellArray(n,2,2) ) then
                            if ( r(3) > cellArray(n,3,1) .AND. r(3) < cellArray(n,3,2) ) then
                                ! gradient sensing #1: difference of particle location and cell center
                                cellPolar(n,1) = cellPolar(n,1) + r(1) - cellCenter(n,1)
                                ! gradient sensing #2: divide cell into 2 bins
                                if ( (r(1) - cellCenter(n,1)) > 0.0_b8 ) then
                                    count2 = count2 + 1
                                else
                                    count1 = count1 + 1
                                end if
                            end if
                        end if
                    end if
                enddo
            else
                exitCount = exitCount + 1
                if ( exitCount > (prtclTotal/10) ) then
                    exit
                end if
            end if
        enddo
        do n = 1, cellTotal
            cellPolar(n,2) = float(count2 - count1)
        enddo
    end subroutine gradientTest


    ! output total cluster polarization
    subroutine wrtPlrTime( run, ntTotal, timePolar)
        implicit none
        integer,  intent(in) :: run, ntTotal
        real(b8), intent(in), dimension(:,:) :: timePolar
        real(b8) :: mean(3)
        integer  :: i, j

        mean = 0.0_b8
        do i = 1, 3
            mean(i) = sum( timePolar(i,:)) / float(ntTotal)
        enddo

        do i = 1, 3
            write(12,"(E16.8)", advance="no") mean(i)
        enddo
        write(12,"(I7)", advance="no") run
        write(12,*) ''
    end subroutine wrtPlrTime


    ! output total cluster polarization for both EC & MW
    subroutine wrtPlrECMW( run, ntTotal, timePolarEC, timePolarMW)
        implicit none
        integer,  intent(in) :: run, ntTotal
        real(b8), intent(in), dimension(:,:) :: timePolarEC, timePolarMW
        real(b8) :: mean(3)
        integer  :: i, j

        mean = 0.0_b8
        do i = 1, 3
            mean(i) = sum( timePolarEC(i,:)) / float(ntTotal)
        enddo
        do i = 1, 3
            write(12,"(E16.8)", advance="no") mean(i)
        enddo
        write(12,"(I7)", advance="no") run
        write(12,*) ''

        mean = 0.0_b8
        do i = 1, 3
            mean(i) = sum( timePolarMW(i,:)) / float(ntTotal)
        enddo
        do i = 1, 3
            write(13,"(E16.8)", advance="no") mean(i)
        enddo
        write(13,"(I7)", advance="no") run
        write(13,*) ''

    end subroutine wrtPlrECMW


    ! output total cluster polarization
    subroutine wrtPlrTotal( run, cellPolar)
        implicit none
        integer,  intent(in) :: run
        real(b8), intent(in), dimension(:,:) :: cellPolar

        write(12,"(E16.8)", advance="no") cellPolar(1,1)
        write(12,"(E16.8)", advance="no") cellPolar(1,2)
        write(12,"(E16.8)", advance="no") cellPolar(1,3)
        write(12,"(I7)", advance="no")    run
        write(12,*) ''
    end subroutine wrtPlrTotal

end module
