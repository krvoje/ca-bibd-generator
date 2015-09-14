program blaser
    implicit none

    real nasumice

    integer, allocatable::c(:)
    integer i,j,n, otklen, doklen
    double precision ZBQLU01,ZBQLUAB

    call posij()

    do i=1,100000
        j=int(nasumice()*8)
        if (j==8)print *,j
    enddo

end program

subroutine posij()
    external ZBQLINI

    call ZBQLINI(time())
end subroutine

real function nasumice()

    double precision ZBQLU01,ZBQLUAB
    integer i

    nasumice=ZBQLU01(i)

end function
