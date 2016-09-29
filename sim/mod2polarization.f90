module polarization
! module for subroutines/functions related to cell polarization

use sysconfig

contains

    ! calculate EC cell polarization, individual polarization vectors are NOT adaptive
    subroutine cellpolarECNonAdpt( cellTotal, prtclTotal, cellArray, edgeList, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal, edgeList(:)
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: clstrCOM(3), cellCOM(3,2), center(3), check, q(3)
        real(b8) :: nCell, nCOM
        integer :: i, j, k

        cellPolar(:,:) = 0.0_b8
        if ( cellTotal == 1 ) then
            return
        end if

        ! calculate Cluster COM
        clstrCOM = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                clstrCOM(j) = clstrCOM(j) + cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
            enddo
        enddo
        clstrCOM(:) = clstrCOM(:) / float(cellTotal)

        do i = 1, cellTotal
            nCell = 0.0_b8
            ! check if cell is on the edge of the cluster
            if ( edgeList(i) == 1 ) then
                q = 0.0_b8
                do j = 1, 3
                    center(j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
                    q(j) = center(j) - clstrCOM(j)
                    if ( abs(q(j)) < (10.0**(-15)) ) then
                        q(j) = 0.0_b8
                    end if
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


    ! calculate 3D MW cell polarization
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
                                ! q = q / sqrt(dot_product(q,q))
                                q = q
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


    ! calculate 2D MW cell polarization
    subroutine cellpolar2DMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
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
                                do k = 1, 2
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
    end subroutine cellpolar2DMW


    ! output total cluster polarization
    subroutine wrtPlrTime( run, ntTotal, timePolar)
        implicit none
        integer,  intent(in) :: run, ntTotal
        real(b8), intent(in), dimension(:,:) :: timePolar
        real(b8) :: mean(3), var(3)
        integer  :: i, j

        mean = 0.0_b8
        var  = 0.0_b8
        do i = 1, 3
            mean(i) = sum( timePolar(i,:)) / float(ntTotal)

            do j = 1, ntTotal
                var(i) = var(i) + (((timePolar(i,j)-mean(i))**2) / float(ntTotal))
            enddo
        enddo

        do i = 1, 3
            write(12,"(E16.8)", advance="no") mean(i)
            write(13,"(E16.8)", advance="no") var(i)
        enddo
        write(12,"(I7)", advance="no") run
        write(13,"(I7)", advance="no") run
        write(12,*) ''
        write(13,*) ''
    end subroutine wrtPlrTime


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
