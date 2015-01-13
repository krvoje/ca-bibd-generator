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

  integer function randomInt(upper)
    integer upper
    randomInt=int(generateRandomNumber()*upper)
  end function randomInt
  
  real function generateRandomNumber()
    use mtmod

    double precision ZBQLU01,ZBQLUAB, ZBQLEXP
    integer i

    !call seedRandomGenerator()

    !!write (*,*) "DEBUG: generateRandomNumber"
    !if(modulo(time(),2) == 0) then
    !   generateRandomNumber=grnd() !<- Ova implementacija Mersenne twistera je koma za ovu svrhu
    !else
    !   generateRandomNumber=ZBQLU01(i)
       !generateRandomNumber=rand()
    !endif
    generateRandomNumber=grnd()
    if (generateRandomNumber==1) generateRandomNumber=0 ! Računamo da se nikad ne dobije 1, jer bi to potrgalo sve
    return
  end function generateRandomNumber

  integer function generateZeroOne() result(rnd)
    real generated

    generated=generateRandomNumber()*99
    if(generated>49) then
       rnd=1
    else
       rnd=0
    endif
    return
  end function generateZeroOne

  subroutine increment(var,val)
    integer var,val
    var = var + val
    return
  end subroutine increment

  subroutine decrement(var,val)
    integer var,val
    var = var - val
    return
  end subroutine decrement
end module utils
