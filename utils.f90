module utils
contains
  subroutine seedRandomGenerator()
    use mtmod

    external ZBQLINI
    !!write (*,*) "DEBUG: seedRandomGenerator"
    call ZBQLINI(time())
    call sgrnd(time())
    call srand(time())
  end subroutine seedRandomGenerator

  real function generateRandomNumber()
    use mtmod

    double precision ZBQLU01,ZBQLUAB, ZBQLEXP
    integer i

    !call seedRandomGenerator()

    !!write (*,*) "DEBUG: generateRandomNumber"
    if(modulo(time(),2) == 0) then
       generateRandomNumber=grnd() !<- Ova implementacija Mersenne twistera je koma za ovu svrhu
    else
       generateRandomNumber=ZBQLU01(i)
       !generateRandomNumber=rand()
    endif
    if (generateRandomNumber==1) generateRandomNumber=0 ! Računamo da se nikad ne dobije 1, jer bi to potrgalo sve
    return

  end function generateRandomNumber
end module utils
